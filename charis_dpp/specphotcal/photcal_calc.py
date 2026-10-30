"""Stellar flux density on the CHARIS wavelength grid
(``charis_photometric_calibration_calculation.pro``).

A library spectrum (Pickles, Castelli–Kurucz, or an empirical file) is
reddened (``av``), smoothed to the CHARIS spectral resolution, colour-
corrected against Vega in the chosen 2MASS/MKO band and scaled to the star's
magnitude, then converted to the requested units.

Fixes vs. IDL: the smoothing kernel width is derived from the resolving
power properly (#82: ``lambda[0]+dloglam`` mixed µm with a log step); Vega is
smoothed on its own grid (#83); unit labels are right (#84); the Kurucz column
follows the README's log g instead of a hard-wired ``g40``; an unmatched
spectral type falls back to the nearest tabulated one with a warning instead
of crashing (#80).
"""

from __future__ import annotations

import logging
import re
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter1d

from charis_dpp.astro.extinction import redden_av
from charis_dpp.astro.units import canonical_units, flam_to
from charis_dpp.astro.zeropoints import filter_curve
from charis_dpp.config.paths import models_dir
from charis_dpp.io.tables import Spectrum, read_kurucz, read_pickles_fits, read_spectrum, read_vega

log = logging.getLogger(__name__)

__all__ = ["stellar_flux_density", "resample_to_resolution", "parse_sptype", "band_filter_name"]

_SPT_RE = re.compile(
    r"^\s*([OBAFGKMLTY])\s*(\d+(?:\.\d+)?)?\s*[-/]?\s*(\d+(?:\.\d+)?)?\s*(I{1,3}V?|IV|V|VI)?", re.I
)
_CLASS_ORDER = "OBAFGKMLTY"


def parse_sptype(spt: str) -> tuple[float, str] | None:
    """``'A5V'`` → (numeric type 25.0, 'V').  Numeric type = 10·class index + subclass."""
    m = _SPT_RE.match(str(spt))
    if not m or not m.group(1):
        return None
    cls = m.group(1).upper()
    sub = float(m.group(2)) if m.group(2) else 5.0
    if m.group(3):
        sub = 0.5 * (sub + float(m.group(3)))
    lum = (m.group(4) or "V").upper()
    return 10.0 * _CLASS_ORDER.index(cls) + sub, lum


def band_filter_name(band: str) -> str:
    """Normalisation filter for a CHARIS band (H for broadband, as in IDL)."""
    return {"J": "J", "H": "H", "K": "Ks", "lowres": "H"}.get(band, "H")


@lru_cache(maxsize=1)
def _pickles_table() -> list[tuple[str, str, float]]:
    txt = (models_dir() / "starmodels" / "pickles" / "AA_README").read_text(errors="replace")
    rows = []
    for line in txt.splitlines():
        m = re.match(r"^\s*(pickles_uk_\d+)\s+(\S+)\s+([\d.]+)\s*$", line)
        if m:
            rows.append((m.group(1), m.group(2), float(m.group(3))))
    if not rows:
        raise RuntimeError("could not parse the Pickles AA_README table")
    return rows


@lru_cache(maxsize=1)
def _kurucz_table() -> list[tuple[str, float, float, str, str]]:
    txt = (models_dir() / "starmodels" / "kurucz" / "AA_README").read_text(errors="replace")
    rows = []
    for line in txt.splitlines():
        m = re.match(r"^\s*(\S+)\s+(\d+)\s+([+-]?[\d.]+)\s+(ckp00_\d+)\[(g\d+)\]\s*$", line)
        if m:
            rows.append((m.group(1), float(m.group(2)), float(m.group(3)), m.group(4), m.group(5)))
    if not rows:
        raise RuntimeError("could not parse the Kurucz AA_README table")
    return rows


def _match_sptype(spt: str, names: list[str]) -> str:
    """Exact (case-insensitive) match, else nearest numeric type with the same luminosity class."""
    for n in names:
        if n.upper() == spt.upper():
            return n
    want = parse_sptype(spt)
    if want is None:
        raise ValueError(f"cannot parse spectral type {spt!r}")
    best, bestd = None, np.inf
    for n in names:
        p = parse_sptype(n)
        if p is None:
            continue
        d = abs(p[0] - want[0]) + (0 if p[1] == want[1] else 5.0)
        if d < bestd:
            best, bestd = n, d
    if best is None:
        raise ValueError(f"no library entry resembles spectral type {spt!r}")
    log.warning("spectral type %r not in library; using nearest entry %r", spt, best)
    return best


def library_spectrum(
    sptype: str, library: str = "pickles", empspectrum: str | Path | None = None
) -> Spectrum:
    """Model spectrum (F_λ, wavelengths in µm) for ``sptype`` from the chosen library."""
    lib = str(library).lower()
    if lib in ("2", "pickles"):
        rows = _pickles_table()
        name = _match_sptype(sptype, [r[1] for r in rows])
        fname = next(r[0] for r in rows if r[1] == name)
        return read_pickles_fits(models_dir() / "starmodels" / "pickles" / f"{fname}.fits")
    if lib in ("1", "kurucz"):
        rows = _kurucz_table()
        name = _match_sptype(sptype, [r[0] for r in rows])
        _, _, logg, fname, gcol = next(r for r in rows if r[0] == name)
        sp = read_kurucz(models_dir() / "starmodels" / "kurucz" / f"{fname}.fits", logg)
        sp.flux = (
            sp.flux * 2.67 / (10 * 206265 * 100) ** 2
        )  # arbitrary (R/D)^2; removed by the mag scaling
        return sp
    if lib in ("3", "empirical", "emp"):
        if empspectrum is None:
            raise ValueError(
                "library='empirical' needs empspectrum=<file with wavelength[um], F_lambda>"
            )
        sp = read_spectrum(empspectrum)
        good = np.isfinite(sp.flux)
        return Spectrum(sp.wavelength_um[good], sp.flux[good], None, "flam", str(empspectrum))
    raise ValueError(f"unknown stellar library {library!r} (kurucz, pickles, empirical)")


def resample_to_resolution(
    wav_in: np.ndarray,
    flux_in: np.ndarray,
    resolving_power: float,
    wav_out: np.ndarray,
    *,
    oversample: int = 8,
) -> np.ndarray:
    """Smooth a spectrum to resolving power ``R`` (Gaussian, FWHM = λ/R) and sample it at ``wav_out``.

    Done on a uniform log-λ grid so a single kernel width applies everywhere.
    """
    wav_in = np.asarray(wav_in, dtype=float)
    flux_in = np.asarray(flux_in, dtype=float)
    good = np.isfinite(flux_in) & (wav_in > 0)
    wav_in, flux_in = wav_in[good], flux_in[good]
    o = np.argsort(wav_in)
    wav_in, flux_in = wav_in[o], flux_in[o]
    lo = max(wav_in[0], np.min(wav_out) * 0.8)
    hi = min(wav_in[-1], np.max(wav_out) * 1.25)
    dlog = (1.0 / resolving_power) / oversample
    grid = np.exp(np.arange(np.log(lo), np.log(hi), dlog))
    f = np.interp(grid, wav_in, flux_in)
    sigma_samples = (1.0 / resolving_power) / 2.35482 / dlog
    fs = gaussian_filter1d(f, sigma_samples, mode="nearest")
    return np.interp(wav_out, grid, fs)


def _band_integral(wav_um: np.ndarray, flux: np.ndarray, filt: str) -> float:
    fw, ft = filter_curve(filt)
    ft = np.where(ft > 0.10, ft, 0.0)  # IDL used only response > 10 %
    f = np.interp(fw, wav_um, flux)
    return float(np.trapezoid(ft * f, fw))


def stellar_flux_density(
    sptype: str,
    mag: float,
    band: str,
    wav_um: np.ndarray,
    *,
    library: str = "pickles",
    av: float = 0.0,
    empspectrum: str | Path | None = None,
    units: str = "mJy",
    resolving_power: float | None = None,
    dloglam: float | None = None,
    filter_name: str | None = None,
) -> np.ndarray:
    """Star flux density per CHARIS channel in ``units``.

    ``mag`` is the Vega magnitude in the normalisation filter (``filter_name``
    or the band's default: J, H, Ks; H for broadband).  ``resolving_power`` (or
    ``dloglam`` = 1/R) sets the smoothing; default R=30 for broadband, 100 for
    J/H/K.
    """
    wav = np.asarray(wav_um, dtype=float)
    units = canonical_units(units)
    if units == "contrast":
        return np.ones_like(wav)
    if resolving_power is None:
        resolving_power = 1.0 / dloglam if dloglam else (30.0 if band == "lowres" else 100.0)
    filt = filter_name or band_filter_name(band)
    model = library_spectrum(sptype, library, empspectrum)
    mflux = model.flux
    if av:
        mflux = redden_av(model.wavelength_um, mflux, float(av))
    fw, _ = filter_curve(filt)
    # smoothed model on the CHARIS grid and on the filter grid
    on_charis = resample_to_resolution(model.wavelength_um, mflux, resolving_power, wav)
    on_filter = resample_to_resolution(model.wavelength_um, mflux, resolving_power, fw)
    vega = read_vega(models_dir() / "starmodels" / "alpha_lyr_stis_005.fits")
    vega_on_filter = resample_to_resolution(vega.wavelength_um, vega.flux, resolving_power, fw)
    m_model = -2.5 * np.log10(
        _band_integral(fw, on_filter, filt) / _band_integral(fw, vega_on_filter, filt)
    )
    scale = 10 ** (-(mag - m_model) / 2.5)
    flam = on_charis * scale
    return flam_to(wav, flam, units)
