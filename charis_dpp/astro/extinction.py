"""Cardelli, Clayton & Mathis (1989) extinction with O'Donnell (1994) optical update — port of astrolib ``ccm_unred``."""

from __future__ import annotations

import numpy as np


def ccm_a_lambda(wave_um: np.ndarray, r_v: float = 3.1) -> np.ndarray:
    """A(λ)/A(V) for wavelengths in µm."""
    x = 1.0 / np.asarray(wave_um, dtype=float)  # inverse microns
    a = np.zeros_like(x)
    b = np.zeros_like(x)

    ir = (x > 0.3) & (x < 1.1)
    a[ir] = 0.574 * x[ir] ** 1.61
    b[ir] = -0.527 * x[ir] ** 1.61

    opt = (x >= 1.1) & (x < 3.3)
    y = x[opt] - 1.82
    c1 = [1.0, 0.104, -0.609, 0.701, 1.137, -1.718, -0.827, 1.647, -0.505]
    c2 = [0.0, 1.952, 2.908, -3.989, -7.985, 11.102, 5.491, -10.805, 3.347]
    a[opt] = np.polyval(c1[::-1], y)
    b[opt] = np.polyval(c2[::-1], y)

    uv = (x >= 3.3) & (x < 8.0)
    xu = x[uv]
    fa = np.zeros_like(xu)
    fb = np.zeros_like(xu)
    hi = xu > 5.9
    yy = xu[hi] - 5.9
    fa[hi] = -0.04473 * yy**2 - 0.009779 * yy**3
    fb[hi] = 0.2130 * yy**2 + 0.1207 * yy**3
    a[uv] = 1.752 - 0.316 * xu - 0.104 / ((xu - 4.67) ** 2 + 0.341) + fa
    b[uv] = -3.090 + 1.825 * xu + 1.206 / ((xu - 4.62) ** 2 + 0.263) + fb

    fuv = (x >= 8.0) & (x <= 11.0)
    yf = x[fuv] - 8.0
    a[fuv] = np.polyval([-0.070, 0.137, -0.628, -1.073][::-1], yf)
    b[fuv] = np.polyval([0.374, -0.420, 4.257, 13.670][::-1], yf)

    valid = (x > 0.3) & (x <= 11.0)
    out = a + b / r_v
    out[~valid] = 0.0
    return out


def ccm_unred(wave_um: np.ndarray, flux: np.ndarray, ebv: float, r_v: float = 3.1) -> np.ndarray:
    """Deredden ``flux`` for E(B-V) = ``ebv`` (negative ``ebv`` reddens)."""
    a_v = r_v * ebv
    a_lambda = a_v * ccm_a_lambda(wave_um, r_v)
    return np.asarray(flux, dtype=float) * 10 ** (0.4 * a_lambda)


def redden_av(wave_um: np.ndarray, flux: np.ndarray, a_v: float, r_v: float = 3.1) -> np.ndarray:
    """Apply A_V magnitudes of extinction (what ``av=`` means in ``specphot_cal``)."""
    return ccm_unred(wave_um, flux, -a_v / r_v, r_v)
