"""Resample model/empirical spectra to the CHARIS channels
(``charis_convert_planetspec_to_lowres_calculation[2].pro``, ``charis_convert_planetphot_to_lowres_calculation.pro``).

Fixes vs. IDL: the output wavelength grid is the one passed in (#202); the
smoothing kernel is derived from the resolving power; the mJy conversion
is applied (the IDL version jumped over it, #203); 2-column files work (#201).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from charis_dpp.astro.units import flam_to
from charis_dpp.astro.zeropoints import filter_curve
from charis_dpp.io.tables import Spectrum, read_spectrum
from charis_dpp.specphotcal.photcal_calc import resample_to_resolution

__all__ = ["resample_to_charis", "load_model_spectrum", "synthetic_photometry"]


def load_model_spectrum(path: str | Path, wavelength_unit: str = "um") -> Spectrum:
    """Model spectrum (wavelength → µm) in whatever flux units the file has."""
    return read_spectrum(path, wavelength_unit=wavelength_unit)


def resample_to_charis(
    wav_in: np.ndarray,
    flux_in: np.ndarray,
    wav_out: np.ndarray,
    err_in: np.ndarray | None = None,
    *,
    resolving_power: float | None = None,
) -> tuple[np.ndarray, np.ndarray | None]:
    """Smooth to the CHARIS resolution and interpolate onto ``wav_out`` (µm).

    ``resolving_power`` defaults to the median ``λ/Δλ`` of ``wav_out``.
    Errors are propagated by resampling the variance.
    """
    wav_out = np.asarray(wav_out, dtype=float)
    if resolving_power is None:
        resolving_power = float(np.median(wav_out[1:] / np.diff(wav_out)))
    f = resample_to_resolution(wav_in, flux_in, resolving_power, wav_out)
    e = None
    if err_in is not None:
        v = resample_to_resolution(
            wav_in, np.asarray(err_in, dtype=float) ** 2, resolving_power, wav_out
        )
        # smoothing over ~N samples reduces the noise by sqrt(N)
        n_per = max(
            float(np.median(np.diff(wav_out)) / np.median(np.diff(np.asarray(wav_in)))), 1.0
        )
        e = np.sqrt(np.clip(v, 0, None) / n_per)
    return f, e


def synthetic_photometry(
    wav_um: np.ndarray, flux: np.ndarray, filter_name: str, *, units_in: str = "flam"
) -> tuple[float, float]:
    """Photon-weighted band flux in mJy and effective wavelength for a spectrum in F_λ (or mJy)."""
    fw, ft = filter_curve(filter_name)
    f = np.interp(fw, wav_um, flux)
    if units_in == "flam":
        f = flam_to(fw, f, "mJy")
    denom = np.trapezoid(fw * ft, fw)
    return float(np.trapezoid(fw * ft * f, fw) / denom), float(
        np.trapezoid(fw * ft, fw) / np.trapezoid(ft, fw)
    )
