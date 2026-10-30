"""Broadband J/H/Ks photometry of a CHARIS spectrum (``charis_est_bbphot.pro``).

Photon-weighted integration ``∫λ T F dλ / ∫λ T dλ`` over the filter curves
smoothed to the CHARIS channel spacing.  All five flux units are accepted
(the IDL version only handled two, #116); errors are propagated in
quadrature (the IDL added them linearly, #118); zero points come from the
single instrument table (#117).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from charis_dpp.astro.units import canonical_units, flam_to, to_flam
from charis_dpp.astro.zeropoints import filter_curve, zeropoint
from charis_dpp.io.tables import read_ascii_columns
from charis_dpp.specphotcal.photcal_calc import resample_to_resolution

__all__ = ["broadband_photometry", "BandFlux"]


@dataclass
class BandFlux:
    band: str
    flux_mjy: float
    err_mjy: float
    lam_eff_um: float
    mag: float
    mag_err: float


def _bands_covered(wav: np.ndarray) -> list[str]:
    lo, hi = wav.min(), wav.max()
    out = []
    if lo < 1.17 and hi > 1.33:
        out.append("J")
    if lo < 1.52 and hi > 1.75:
        out.append("H")
    if lo < 2.05 and hi > 2.3:
        out.append("Ks")
    return out


def broadband_photometry(
    wav_um: np.ndarray,
    flux: np.ndarray,
    err: np.ndarray | None = None,
    *,
    units: str = "mJy",
    bands: list[str] | None = None,
) -> dict[str, BandFlux]:
    """Return per-band flux (mJy), error, effective wavelength and Vega magnitude."""
    wav = np.asarray(wav_um, dtype=float)
    units = canonical_units(units)
    flam = to_flam(wav, np.asarray(flux, dtype=float), units)
    eflam = (
        to_flam(wav, np.asarray(err, dtype=float), units)
        if err is not None
        else np.zeros_like(flam)
    )
    dl = np.median(np.diff(wav))
    R = float(np.median(wav) / dl)
    out: dict[str, BandFlux] = {}
    for b in bands or _bands_covered(wav):
        fw, ft = filter_curve(b)
        tf = resample_to_resolution(fw, ft, R, wav)
        good = tf > 0.10
        if good.sum() < 3:
            continue
        w, t, f, e = wav[good], tf[good], flam[good], eflam[good]
        denom = np.trapezoid(w * t, w)
        fl = np.trapezoid(w * t * f, w) / denom
        ef = np.sqrt(np.trapezoid((w * t * e) ** 2, w) / denom**2 * (w[-1] - w[0]) / max(len(w), 1))
        leff = np.trapezoid(w * t, w) / np.trapezoid(t, w)
        fl_mjy = float(flam_to(np.array([leff]), np.array([fl]), "mJy")[0])
        ef_mjy = float(flam_to(np.array([leff]), np.array([ef]), "mJy")[0])
        zp = zeropoint(b)
        mag = -2.5 * np.log10(fl_mjy / (1e3 * zp.f0_jy))
        mag_err = 1.0857 * ef_mjy / fl_mjy if fl_mjy > 0 else np.nan
        out[b] = BandFlux(b, fl_mjy, ef_mjy, float(leff), float(mag), float(mag_err))
    return out


def broadband_photometry_file(path, units: str = "mJy") -> dict[str, BandFlux]:
    arr = read_ascii_columns(path)
    err = arr[:, 2] if arr.shape[1] > 2 else None
    return broadband_photometry(arr[:, 0], arr[:, 1], err, units=units)
