"""Wavelength tables and channel sets.

Replaces ``cals/get_charis_wvlh.pro`` (which returned **nm**; everything here
is **µm**) and every hard-coded channel list in the IDL tree (``refslice=10``,
``[0:4]/[7:13]/[16:20]``, telluric ``[5,6,14,15]``…).  Channel sets are derived
from the wavelength grid and the band definitions in ``instrument.toml`` so
they are valid for the 22-channel broadband mode and the 15/20/17-channel
J/H/K modes alike.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np
from astropy.io.fits import Header

from charis_dpp.config.instrument import Instrument, instrument_for, load_instrument
from charis_dpp.config.paths import package_data


def band_from_header(header: Header) -> str:
    """Band from ``CAL_BAND`` (preferred) or ``FILTNAME``."""
    from charis_dpp.config.dataset import canonical_band

    for key in ("CAL_BAND", "FILTNAME"):
        if key in header:
            b = canonical_band(str(header[key]).strip())
            if b:
                return b
    raise KeyError("neither CAL_BAND nor FILTNAME identifies a CHARIS band in this header")


@lru_cache(maxsize=8)
def _table_um(band: str) -> np.ndarray:
    spec = load_instrument().bands[band]
    nm = np.loadtxt(package_data(spec["wavelengths_file"]), dtype=float)
    return np.atleast_1d(nm) * 1e-3


def wavelengths_um(band_or_header: str | Header, inst: Instrument | None = None) -> np.ndarray:
    """Channel wavelengths in µm for a band name or a cube header."""
    from charis_dpp.config.dataset import canonical_band

    if isinstance(band_or_header, Header):
        band = band_from_header(band_or_header)
    else:
        band = canonical_band(band_or_header)
        if band is None:
            raise ValueError(f"unknown CHARIS band {band_or_header!r}")
    return _table_um(band).copy()


def reference_channel(
    wav_um: np.ndarray, inst: Instrument | None = None, band: str = "lowres"
) -> int:
    """Channel nearest the band's geometric-mean wavelength (10 for broadband, like IDL's ``refslice``)."""
    spec = load_instrument().bands[band]
    ref = spec.get("reference_channel", "auto")
    if ref != "auto":
        return int(ref)
    gm = np.sqrt(wav_um[0] * wav_um[-1])
    return int(np.argmin(np.abs(wav_um - gm)))


def channels_in(wav_um: np.ndarray, lo_um: float, hi_um: float) -> np.ndarray:
    return np.flatnonzero((wav_um >= lo_um) & (wav_um <= hi_um))


def telluric_channels(
    wav_um: np.ndarray, inst: Instrument | None = None, band: str = "lowres"
) -> np.ndarray:
    """Channels inside the telluric water bands defined for the band."""
    spec = load_instrument().bands[band]
    idx: list[np.ndarray] = [channels_in(wav_um, lo, hi) for lo, hi in spec.get("telluric_um", [])]
    return np.unique(np.concatenate(idx)) if idx else np.zeros(0, dtype=int)


def band_channels(
    wav_um: np.ndarray, inst: Instrument | None = None, band: str = "lowres"
) -> dict[str, np.ndarray]:
    """``{"J": idx, "H": idx, "K": idx}`` sub-band channel sets for this grid."""
    spec = load_instrument().bands[band]
    return {k: channels_in(wav_um, lo, hi) for k, (lo, hi) in spec.get("sub_bands_um", {}).items()}


def fwhm_px(
    wav_um: np.ndarray, mjd: float | None = None, inst: Instrument | None = None
) -> np.ndarray:
    """λ/D in pixels per channel (what the IDL code called ``fwhm``)."""
    if inst is None:
        if mjd is None:
            raise ValueError("fwhm_px needs an Instrument or an MJD")
        inst = instrument_for(mjd)
    return np.asarray(inst.lambda_over_d_px(wav_um), dtype=float)
