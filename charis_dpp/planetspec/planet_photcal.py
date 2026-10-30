"""Scale a model planet spectrum to a contrast relative to the star
(``charis_planet_photometric_calibration_calculation.pro``).

The model (any F_λ scale) is smoothed to CHARIS resolution, converted to
F_ν [mJy] and scaled so that its band-averaged flux is ``10^contrast`` times
the star's band-averaged flux (``FSTAR_<l>`` from the calibrated cube), where
the band weighting is the filter response resampled to the CHARIS channels.

Fixes vs. IDL: no dependency on a dotted file name (#197), negative filter
responses are zeroed without the last-pixel bug (#198), unknown bands fall
back to H (#199).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from astropy.io.fits import Header

from charis_dpp.astro.units import flam_to
from charis_dpp.astro.zeropoints import filter_curve
from charis_dpp.planetspec.convert import load_model_spectrum, resample_to_charis
from charis_dpp.specphotcal.photcal_calc import band_filter_name

__all__ = ["planet_spectrum_mjy", "star_spectrum_from_header", "band_weights"]


def star_spectrum_from_header(h1: Header, nlam: int) -> np.ndarray:
    return np.array([float(h1.get(f"FSTAR_{l}", np.nan)) for l in range(nlam)])


def band_weights(wav_um: np.ndarray, band: str, filter_name: str | None = None) -> np.ndarray:
    fw, ft = filter_curve(filter_name or band_filter_name(band))
    R = float(np.median(wav_um[1:] / np.diff(wav_um)))
    w = resample_to_charis(fw, ft, wav_um, resolving_power=R)[0]
    return np.clip(w, 0.0, None)


def planet_spectrum_mjy(
    model: str | Path | np.ndarray,
    star_mjy: np.ndarray,
    wav_um: np.ndarray,
    band: str,
    log_contrast: float,
    *,
    filter_name: str | None = None,
    model_units: str = "flam",
) -> np.ndarray:
    """Planet F_ν [mJy] per channel at ``10^log_contrast`` of the star in the band.

    ``model`` is a file (Cruz/Montreal/Bonnefoy FITS or 2-column text, wavelength
    in µm) or an array ``(wav_um, flux)`` already on any grid.
    """
    wav_um = np.asarray(wav_um, dtype=float)
    if isinstance(model, np.ndarray):
        mw, mf = model[0], model[1]
    else:
        sp = load_model_spectrum(model)
        mw, mf = sp.wavelength_um, sp.flux
    f_on_grid, _ = resample_to_charis(mw, mf, wav_um)
    f_nu = flam_to(wav_um, f_on_grid, "mJy") if model_units == "flam" else f_on_grid
    w = band_weights(wav_um, band, filter_name)
    if w.sum() <= 0:
        raise ValueError("filter response does not overlap the CHARIS channels")
    model_band = np.trapezoid(f_nu * w, wav_um) / np.trapezoid(w, wav_um)
    star_band = np.trapezoid(np.asarray(star_mjy) * w, wav_um) / np.trapezoid(w, wav_um)
    scale = 10 ** float(log_contrast) * star_band / model_band
    return f_nu * scale
