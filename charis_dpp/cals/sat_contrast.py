"""Satellite-spot attenuation (``cals/get_charis_sat_contrast.pro``).

The spot/star flux ratio at wavelength λ for a modulation amplitude ``mod``
(nm) is ``10^(-0.4 m(mod)) (λ_ref/λ)²`` where ``m`` is the calibrated contrast
magnitude at ``λ_ref = 1.55 µm``.  The magnitude table changed at MJD 57995
(Aug 2017); both tables live in ``instrument.toml`` and are selected by date.
"""

from __future__ import annotations

import numpy as np

from charis_dpp.config.instrument import Instrument, instrument_for


def sat_contrast(
    modulation_nm: float,
    wavelengths_um: np.ndarray,
    mjd: float | None = None,
    *,
    inst: Instrument | None = None,
    manual: float | np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(attenuation, error)`` per channel: spot flux / star flux.

    ``manual`` gives the contrast in magnitudes directly (scalar, or one per
    channel — the IDL version silently used only the first element, BUGS.md
    #68/#122) and is applied without the wavelength scaling.
    """
    wav = np.asarray(wavelengths_um, dtype=float)
    if inst is None:
        if mjd is None:
            raise ValueError("sat_contrast needs an Instrument or an MJD")
        inst = instrument_for(mjd)
    spots = inst.sat_spots
    if manual is not None:
        mags = np.broadcast_to(np.asarray(manual, dtype=float), wav.shape)
        err_mag = float(spots.manual_err_mag)
        atten = 10 ** (-0.4 * mags)
        return atten, atten * 0.4 * np.log(10) * err_mag
    ref = float(spots.reference_wavelength_um)
    table = {float(k): float(v) for k, v in spots.contrast_mag.items()}
    errs = {float(k): float(v) for k, v in spots.contrast_err_mag.items()}
    mod = float(modulation_nm)
    if mod in table:
        mag = table[mod]
        emag = errs.get(mod, float(spots.manual_err_mag))
    else:
        # scale the 25 nm calibration as (mod/25)^2 in flux
        base = 25.0 if 25.0 in table else min(table)
        mag = table[base] - 2.5 * np.log10((mod / base) ** 2)
        emag = errs.get(base, float(spots.manual_err_mag))
    atten = 10 ** (-0.4 * mag) * (ref / wav) ** 2
    return atten, atten * 0.4 * np.log(10) * emag


def modulation_from_header(h0, h1=None, default: float | None = None) -> float:
    """``X_GRDAMP`` (µm in the header) → nm; falls back to ``default`` when absent/'unavailable'."""
    for h in (h1, h0):
        if h is None or "X_GRDAMP" not in h:
            continue
        v = h["X_GRDAMP"]
        try:
            return float(v) * 1e3
        except (TypeError, ValueError):
            continue
    if default is None:
        raise KeyError("X_GRDAMP not usable in headers and no modulation amplitude given")
    return float(default)
