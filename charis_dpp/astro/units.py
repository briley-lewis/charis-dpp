"""Flux-density unit conversions.

The IDL code carried six numbered unit codes (0 mJy, 1 Jy, 2 W/m²/µm,
3 erg/s/cm²/Å, 4 erg/s/cm²/Hz, 5 contrast) and mislabelled four of them in
its output (BUGS.md #84).  Here units are strings and every conversion goes
through F_λ in erg/s/cm²/Å.
"""

from __future__ import annotations

import numpy as np

C_UM_S = 2.99792458e14  # speed of light, µm/s

UNIT_NAMES = {
    "mJy": "mJy",
    "jy": "Jy",
    "w/m^2/um": "W/m^2/um",
    "ergs/s/cm^2/a": "ergs/s/cm^2/A",
    "flam": "ergs/s/cm^2/A",
    "ergs/s/cm^2/hz": "ergs/s/cm^2/Hz",
    "fnu": "ergs/s/cm^2/Hz",
    "contrast": "contrast",
}
IDL_UNIT_CODES = {
    0: "mJy",
    1: "Jy",
    2: "W/m^2/um",
    3: "ergs/s/cm^2/A",
    4: "ergs/s/cm^2/Hz",
    5: "contrast",
}


def canonical_units(u: str | int) -> str:
    if isinstance(u, int | np.integer):
        return IDL_UNIT_CODES[int(u)]
    key = str(u).strip()
    if key == "mJy":
        return "mJy"
    low = key.lower()
    if low in UNIT_NAMES:
        return UNIT_NAMES[low]
    if low == "mjy":
        return "mJy"
    raise ValueError(f"unknown flux units {u!r}")


def flam_to(wav_um: np.ndarray, flam: np.ndarray, units: str) -> np.ndarray:
    """erg/s/cm²/Å → ``units``."""
    units = canonical_units(units)
    wav = np.asarray(wav_um, dtype=float)
    flam = np.asarray(flam, dtype=float)
    if units == "ergs/s/cm^2/A":
        return flam
    if units == "W/m^2/um":
        return flam * 1e4 * 1e4 * 1e-7  # Å→µm, cm²→m², erg→J
    fnu = flam * (wav**2 / C_UM_S) * 1e4  # erg/s/cm²/Hz  (λ² dλ/dν with λ in µm → Å)
    if units == "ergs/s/cm^2/Hz":
        return fnu
    if units == "Jy":
        return fnu * 1e23
    if units == "mJy":
        return fnu * 1e26
    if units == "contrast":
        return np.ones_like(flam)
    raise ValueError(units)


def to_flam(wav_um: np.ndarray, flux: np.ndarray, units: str) -> np.ndarray:
    """``units`` → erg/s/cm²/Å (inverse of :func:`flam_to`)."""
    units = canonical_units(units)
    wav = np.asarray(wav_um, dtype=float)
    flux = np.asarray(flux, dtype=float)
    if units == "ergs/s/cm^2/A":
        return flux
    if units == "W/m^2/um":
        return flux / (1e4 * 1e4 * 1e-7)
    if units == "mJy":
        fnu = flux * 1e-26
    elif units == "Jy":
        fnu = flux * 1e-23
    elif units == "ergs/s/cm^2/Hz":
        fnu = flux
    else:
        raise ValueError(f"cannot convert {units} to F_lambda")
    return fnu / ((wav**2 / C_UM_S) * 1e4)
