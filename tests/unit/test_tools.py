"""Analysis tools: model-grid and empirical fits, spectral indices, background track, image quality."""

from __future__ import annotations

import numpy as np
from astropy.io import fits
from charis_dpp.astro.units import flam_to
from charis_dpp.io.catalogs import StarInfo
from charis_dpp.io.wavelengths import wavelengths_um
from charis_dpp.planetspec.convert import resample_to_charis
from charis_dpp.tools.astrometry import background_track
from charis_dpp.tools.atmosfit import R_JUP_PC, fit_model_grid, good_channels
from charis_dpp.tools.diskmodel_resample import resample_disk_model
from charis_dpp.tools.empcomp import (
    fit_empirical_library,
    gravity_class,
    luminosity_from_k,
    spectral_type_number,
    spt_gravity_indices,
)
from charis_dpp.tools.imagequality import encircled_energy, plane_fit_annulus


def _blackbody_flam(wav_um, teff):
    h, c, k = 6.626e-27, 2.998e10, 1.381e-16
    lam = wav_um * 1e-4
    return 2 * h * c**2 / lam**5 / (np.exp(h * c / (lam * k * teff)) - 1) * 1e-8  # per Å (cgs-ish)


def test_fit_model_grid_recovers_radius(tmp_path, rng):
    wav = wavelengths_um("lowres")
    grid = tmp_path / "grid"
    grid.mkdir()
    mw = np.linspace(0.95, 5.5, 20000)
    names = [
        "lte015-4.5-0.0a+0.0.BT-Settl.txt",
        "lte020-4.0-0.0a+0.0.BT-Settl.txt",
        "lte012-5.0-0.0a+0.0.BT-Settl.txt",
    ]
    for n, t in zip(names, (1500, 2000, 1200), strict=True):
        np.savetxt(grid / n, np.column_stack([mw * 1e4, _blackbody_flam(mw, t)]), fmt="%.6e")
    # "observed": the 1500 K model, 1.4 R_Jup at 25 pc, 5% noise
    d, radius = 25.0, 1.4
    pred, _ = resample_to_charis(mw, _blackbody_flam(mw, 1500), wav)
    obs = flam_to(wav, pred * (radius * R_JUP_PC / d) ** 2, "mJy")
    err = 0.05 * obs
    obs = obs + rng.normal(0, 1, obs.size) * err
    sp = tmp_path / "spec.dat"
    np.savetxt(sp, np.column_stack([wav, obs, err]))
    res = fit_model_grid(sp, grid_dir=grid, distance_pc=d, outdir=tmp_path / "fit", plot=True)
    assert res[0].name.startswith("lte015") and res[0].teff == 1500
    assert abs(res[0].radius_rjup - radius) < 0.1
    assert 0.3 < res[0].reduced < 3.0
    assert (tmp_path / "fit" / "fitoutcome.txt").exists() and (
        tmp_path / "fit" / f"{res[0].name}.png"
    ).exists()
    assert good_channels(wav, "h").size < good_channels(wav).size < wav.size


def test_spectral_types_and_indices():
    assert spectral_type_number("M8.5") == 8.5
    assert spectral_type_number("L1 low") == 11.0
    assert spectral_type_number("<T3:") == 23.0
    assert spectral_type_number("K7") == -3.0
    assert np.isnan(spectral_type_number("??"))
    assert (
        gravity_class("L1 low") == 2 and gravity_class("M9 int") == 1 and gravity_class("L5") == 0
    )
    assert gravity_class("L2 very low") == 3
    wav = wavelengths_um("lowres")
    idx = spt_gravity_indices(wav, fnu=np.ones(wav.size))
    assert all(np.isfinite(v) for v in idx.values())
    assert 0.5 < idx["H2K"] < 1.5 and 0.5 < idx["H-cont"] < 1.5
    logl, mbol = luminosity_from_k(14.0, 40.0, spt_number=12.0)
    assert -5 < logl < -2 and 10 < mbol < 20


def test_fit_empirical_library(tmp_path, rng):
    wav = wavelengths_um("lowres")
    lib = tmp_path / "lib"
    lib.mkdir()
    mw = np.linspace(0.9, 2.6, 4000)
    stems = ["ObjA_NIR_SpeX_Test2020", "ObjB_NIR_SpeX_Test2020"]
    temps = (1600, 2400)
    for s, t in zip(stems, temps, strict=True):
        f = _blackbody_flam(mw, t)
        fits.PrimaryHDU(np.vstack([mw, f, 0.02 * f]).astype(np.float32)).writeto(lib / f"{s}.fits")
    cols = ["Obj A"] + ["x"] * 4 + ["L1 low"] + ["x"] * 30 + [stems[0], "note"]
    cols2 = ["Obj B"] + ["x"] * 4 + ["M6"] + ["x"] * 30 + [stems[1], "note"]
    (lib / "montreal_library_full.tsv").write_text("\t".join(cols) + "\n" + "\t".join(cols2) + "\n")
    pred, _ = resample_to_charis(mw, _blackbody_flam(mw, 1600), wav)
    obs = flam_to(wav, pred * 1e-20, "mJy")
    err = 0.05 * obs
    sp = tmp_path / "spec.dat"
    np.savetxt(sp, np.column_stack([wav, obs + rng.normal(0, 1, obs.size) * err, err]))
    res = fit_empirical_library(sp, library_dir=lib, outdir=tmp_path / "emp", plot=True)
    assert res[0].name == "Obj A" and res[0].spt == "L1 low" and res[0].gravity == 2
    assert res[0].spt_number == 11.0 and res[0].reduced < res[1].reduced
    assert (tmp_path / "emp" / "chisqdist.png").exists()


def test_background_track_follows_proper_motion(tmp_path):
    star = StarInfo("test", 150.0, 20.0, plx_mas=20.0, pmra_mas_yr=-300.0, pmdec_mas_yr=100.0)
    tr = background_track(star, "2018-01-01", nyears=1.0, plot=tmp_path / "pm.png")
    # a background object drifts opposite to the star's proper motion (± the parallax ellipse)
    assert abs(tr.d_east_arcsec[-1] - 0.300) < 0.03 and abs(tr.d_north_arcsec[-1] + 0.100) < 0.03
    assert (tmp_path / "pm.png").exists()
    pf = tmp_path / "pos.txt"
    pf.write_text("0.5 0.5 0.01 0.01 2018-01-01\n0.8 0.4 0.01 0.01 2019-01-01\n")
    tr2 = background_track(star, "2018-01-01", nyears=1.0, posfile=pf)
    assert tr2.predicted_at_epochs.shape == (2, 2) and tr2.chisq is not None and tr2.chisq < 20


def test_encircled_energy_and_plane_fit(tmp_path):
    yy, xx = np.mgrid[0:61, 0:61]
    sig = 2.0
    psf = np.exp(-0.5 * ((xx - 30) ** 2 + (yy - 30) ** 2) / sig**2)
    r = encircled_energy(psf[None], dpsf=31, plot=tmp_path / "ee.png")
    ee = r.ee[0]
    assert np.all(np.diff(ee) >= -1e-9) and abs(ee[-1] - 1.0) < 1e-6
    # analytic EE at r = 2 sigma: 1 - exp(-2)
    j = int(np.argmin(np.abs(r.radii_px - 2 * sig)))
    assert abs(ee[j] - (1 - np.exp(-2))) < 0.03
    im = 3.0 + 0.1 * xx - 0.05 * yy + 50 * np.exp(-0.5 * ((xx - 40) ** 2 + (yy - 30) ** 2) / 1.5**2)
    sub, plane = plane_fit_annulus(im, (40, 30), 3.0)
    assert abs(sub[30, 40] - 50.0) < 0.5 and abs(plane[30, 40] - (3 + 4 - 1.5)) < 0.05


def test_resample_disk_model():
    yy, xx = np.mgrid[0:121, 0:121]
    model = np.exp(-0.5 * ((xx - 60) ** 2 + (yy - 60) ** 2) / 8.0**2)
    out = resample_disk_model(model, 0.0095, 0.019, out_size=61)
    assert out.shape == (61, 61)
    prof = out[30, :]
    half = np.flatnonzero(prof > 0.5)
    assert abs((half[-1] - half[0]) - 2.355 * 4.0) < 2.0  # width halves on the coarser grid
