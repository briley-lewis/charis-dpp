"""Photometric zero points and filter curves (one source of truth, from ``instrument.toml``)."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np

from charis_dpp.config.instrument import load_instrument
from charis_dpp.config.paths import package_data

_ALIASES = {
    "KS": "Ks",
    "LP": "Lp",
    "L'": "Lp",
    "MS": "Ms",
    "BRA": "BrA",
    "[3.1]": "H2O",
    "[3.3]": "PAH",
    "[4.05]": "BrA",
}


def canonical_filter(name: str) -> str:
    n = str(name).strip()
    zp = load_instrument().photometry["zeropoints"]
    if n in zp:
        return n
    up = n.upper()
    if up in _ALIASES:
        return _ALIASES[up]
    for k in zp:
        if k.upper() == up:
            return k
    raise KeyError(f"unknown filter {name!r}; known: {sorted(zp)}")


@dataclass(frozen=True)
class ZeroPoint:
    name: str
    f0_jy: float
    lam_um: float
    dlam_um: float
    system: str


def zeropoint(name: str) -> ZeroPoint:
    k = canonical_filter(name)
    d = load_instrument().photometry["zeropoints"][k]
    return ZeroPoint(
        k, float(d["f0_jy"]), float(d["lam_um"]), float(d["dlam_um"]), str(d.get("system", ""))
    )


def bandwidth_um(name: str) -> float:
    return zeropoint(name).dlam_um


def mag_to_mjy(mag: float, name: str) -> float:
    return 1e3 * zeropoint(name).f0_jy * 10 ** (-0.4 * mag)


def mjy_to_mag(f_mjy: float, name: str) -> float:
    return -2.5 * np.log10(f_mjy / (1e3 * zeropoint(name).f0_jy))


@lru_cache(maxsize=32)
def filter_curve(name: str) -> tuple[np.ndarray, np.ndarray]:
    """``(wavelength_um, transmission_fraction)`` sorted by wavelength, non-negative."""
    k = canonical_filter(name)
    filters = load_instrument().photometry["filters"]
    if k not in filters:
        raise KeyError(f"no filter response curve for {k!r}")
    spec = filters[k]
    path = package_data(spec["file"])
    hdr = spec.get("header_lines", 0)
    rows = []
    for i, line in enumerate(path.read_text().splitlines()):
        s = line.strip()
        if not s:
            continue
        if hdr == "#" and s.startswith("#"):
            continue
        if isinstance(hdr, int) and i < hdr:
            continue
        parts = s.replace(",", " ").split()
        try:
            rows.append((float(parts[0]), float(parts[1])))
        except (ValueError, IndexError):
            continue
    arr = np.array(rows, dtype=float)
    wav, tr = arr[:, 0], arr[:, 1]
    unit = spec["wavelength_unit"]
    if unit == "nm":
        wav = wav * 1e-3
    elif unit == "cm-1":
        wav = 1e4 / wav
    elif unit == "angstrom":
        wav = wav * 1e-4
    if spec.get("transmission") == "percent":
        tr = tr / 100.0
    tr = np.clip(tr, 0.0, None)
    order = np.argsort(wav)
    wav, tr = wav[order], tr[order]
    # drop duplicate wavelengths (Hband.dat has repeated rows)
    keep = np.concatenate([[True], np.diff(wav) > 0])
    return wav[keep], tr[keep]
