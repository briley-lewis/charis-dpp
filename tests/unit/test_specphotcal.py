from __future__ import annotations

import numpy as np
import pytest
from charis_dpp.astro.units import flam_to, to_flam
from charis_dpp.astro.zeropoints import filter_curve, zeropoint
from charis_dpp.io.wavelengths import wavelengths_um
from charis_dpp.specphotcal.bbphot import broadband_photometry
from charis_dpp.specphotcal.photcal_calc import (
    parse_sptype,
    resample_to_resolution,
    stellar_flux_density,
)
from charis_dpp.specphotcal.specphot_cal import measure_spot_fluxes
from charis_dpp.specphotcal.synth_psf import gaussian_psf_cube
from charis_dpp.synthetic import psf_image


def test_units_roundtrip():
    wav = np.array([1.2, 1.65, 2.2])
    flam = np.array([1e-12, 5e-13, 2e-13])
    for u in ("mJy", "Jy", "W/m^2/um", "ergs/s/cm^2/Hz", "ergs/s/cm^2/A"):
        assert np.allclose(to_flam(wav, flam_to(wav, flam, u), u), flam)
    # 1 Jy at 1.65 um in F_lambda: 1e-23 * c / lambda^2 with lambda in Angstrom
    f = to_flam(np.array([1.65]), np.array([1.0]), "Jy")[0]
    assert f == pytest.approx(1e-23 * 2.99792458e18 / 16500.0**2, rel=1e-6)


@pytest.mark.parametrize(
    "s,num,lum",
    [("A5V", 25.0, "V"), ("K3.5III", 53.5, "III"), ("M2", 62.0, "V"), ("B5-7V", 16.0, "V")],
)
def test_parse_sptype(s, num, lum):
    assert parse_sptype(s) == (num, lum)


def test_filter_curves_normalised():
    for f in ("J", "H", "Ks"):
        w, t = filter_curve(f)
        assert np.all(np.diff(w) > 0) and 0.5 < t.max() <= 1.05
        assert abs(w[np.argmax(t > 0.5)] - zeropoint(f).lam_um) < 0.3


def test_resample_conserves_flat():
    w = np.linspace(1.0, 2.5, 5000)
    out = resample_to_resolution(w, np.ones_like(w), 30.0, np.linspace(1.2, 2.3, 22))
    assert np.allclose(out, 1.0, atol=1e-6)


def test_stellar_flux_reproduces_zero_point():
    """An A0V star with H = 0 must have the Vega H-band flux density (within the Pickles vs STIS spread)."""
    w = wavelengths_um("lowres")
    for lib in ("pickles", "kurucz"):
        f = stellar_flux_density("A0V", 0.0, "lowres", w, library=lib, units="mJy")
        assert abs(np.interp(1.65, w, f) / (1e3 * zeropoint("H").f0_jy) - 1) < 0.05, lib
    fj = stellar_flux_density("A0V", 0.0, "J", wavelengths_um("J"), library="pickles", units="Jy")
    assert abs(np.interp(1.25, wavelengths_um("J"), fj) / zeropoint("J").f0_jy - 1) < 0.06
    # a fainter star scales by 10^-0.4 m
    f5 = stellar_flux_density("A0V", 5.0, "lowres", w, library="pickles", units="mJy")
    f0 = stellar_flux_density("A0V", 0.0, "lowres", w, library="pickles", units="mJy")
    assert np.allclose(f5 / f0, 10**-2.0, rtol=1e-6)
    # unmatched type falls back to the nearest one
    fa6 = stellar_flux_density("A6V", 0.0, "lowres", w, library="pickles")
    assert np.all(np.isfinite(fa6))


def test_broadband_photometry_flat_spectrum():
    w = wavelengths_um("lowres")
    flux = np.full(w.size, 1e-13)  # flat F_lambda
    res = broadband_photometry(w, flux, 0.05 * flux, units="ergs/s/cm^2/A")
    assert set(res) == {"J", "H", "Ks"}
    for b, r in res.items():
        expect = flam_to(np.array([r.lam_eff_um]), np.array([1e-13]), "mJy")[0]
        assert abs(r.flux_mjy / expect - 1) < 0.02, b
        assert abs(r.lam_eff_um - zeropoint(b).lam_um) < 0.1
        assert r.err_mjy > 0 and np.isfinite(r.mag)


def test_gaussian_psf_cube_normalised():
    psf = gaussian_psf_cube(np.array([2.5, 3.5]), 21)
    assert np.allclose(psf.sum(axis=(1, 2)), 1.0)
    assert psf[1, 10, 10] < psf[0, 10, 10]


def test_spot_photometry_recovers_flux(rng):
    shape = (101, 101)
    nlam = 3
    fwhm = np.array([2.5, 3.0, 3.5])
    pos = np.array([[[30.3, 70.2], [70.6, 70.4], [70.1, 30.2], [30.4, 30.5]]] * nlam)
    cube = np.zeros((nlam,) + shape)
    for il in range(nlam):
        for q in range(4):
            cube[il] += psf_image(shape, pos[il, q, 0], pos[il, q, 1], fwhm[il], 100.0)
        cube[il] += rng.normal(0, 0.5, shape)
    ph = measure_spot_fluxes(cube, pos, 0.5 * fwhm)
    # flux inside r = fwhm/2 of a Gaussian: 1 - 2^(-1/2)... use the direct integral
    sig = fwhm / 2.35482
    expect = 100.0 * 2 * np.pi * sig**2 * (1 - np.exp(-0.5 * (0.5 * fwhm / sig) ** 2))
    assert np.allclose(ph.flux, expect, rtol=0.06)  # pixelated Gaussian vs analytic
    assert np.all(ph.err > 0) and np.all(ph.err < 0.1 * ph.flux)
    ph2 = measure_spot_fluxes(cube, pos, 0.5 * fwhm, spot_select="left")
    assert np.allclose(ph2.flux, expect, rtol=0.06)
