"""Robust statistics and frame combination.

Verbatim ports of the IDL Astronomy Library routines the pipeline relies on
(``resistant_mean``, ``robust_sigma``, ``robust_poly_fit``) so results match
the IDL reduction to round-off, plus ``combine_frames`` replacing
``charis_combinefits.pro`` / ``medfitsme``.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "idl_median",
    "resistant_mean",
    "robust_sigma",
    "robust_poly_fit",
    "combine_frames",
    "mad_sigma",
]

_SIGCOEFF = np.array([-0.15405, 0.90723, -0.23584, 0.020142])  # astrolib polynomial


def idl_median(x: np.ndarray, axis: int | None = None, even: bool = True) -> np.ndarray | float:
    """NaN-ignoring median.  ``even=False`` mimics IDL ``median()`` without ``/EVEN``
    (for an even count it returns the upper of the two middle values)."""
    x = np.asarray(x, dtype=float)
    if even:
        return np.nanmedian(x, axis=axis)
    if axis is None:
        v = np.sort(x[np.isfinite(x)])
        return np.nan if v.size == 0 else float(v[v.size // 2])
    x = np.moveaxis(x, axis, -1)
    n = np.isfinite(x).sum(axis=-1)
    s = np.sort(np.where(np.isfinite(x), x, np.inf), axis=-1)
    idx = np.clip(n // 2, 0, s.shape[-1] - 1)
    out = np.take_along_axis(s, idx[..., None], axis=-1)[..., 0]
    return np.where(n > 0, out, np.nan)


def mad_sigma(x: np.ndarray, axis: int | None = None) -> np.ndarray | float:
    """Median absolute deviation scaled to σ (÷0.6745)."""
    x = np.asarray(x, dtype=float)
    med = np.nanmedian(x, axis=axis, keepdims=axis is not None)
    return np.nanmedian(np.abs(x - med), axis=axis) / 0.6745


def _poly(x: float, c: np.ndarray) -> float:
    return float(np.polyval(c[::-1], x))


def resistant_mean(
    y: np.ndarray,
    cut: float = 3.0,
    axis: int | None = None,
    *,
    return_all: bool = False,
):
    """Astrolib ``resistant_mean``: outlier-resistant mean (and σ of the mean).

    Returns ``mean`` (or ``(mean, sigma_mean, num_rejected)`` when
    ``return_all``).  NaNs are ignored.  With ``axis`` the statistic is taken
    along that axis of ``y``.
    """
    y = np.asarray(y, dtype=float)
    if axis is None:
        y = y[np.isfinite(y)].ravel()
        if y.size == 0:
            return (np.nan, np.nan, 0) if return_all else np.nan
        y = y[None, :]
        res = _resistant_mean_axis(y, cut)
        m, s, n = (r[0] for r in res)
        return (float(m), float(s), int(n)) if return_all else float(m)
    y = np.moveaxis(y, axis, -1)
    shp = y.shape[:-1]
    res = _resistant_mean_axis(y.reshape(-1, y.shape[-1]), cut)
    m, s, n = (r.reshape(shp) for r in res)
    return (m, s, n) if return_all else m


def _resistant_mean_axis(y2: np.ndarray, cut: float):
    """Vectorised over rows of a 2-D array; NaN-aware."""
    finite = np.isfinite(y2)
    npts = finite.sum(axis=1)
    ymed = np.nanmedian(y2, axis=1)
    absdev = np.abs(y2 - ymed[:, None])
    mad = np.nanmedian(absdev, axis=1) / 0.6745
    small = mad < 1e-24
    if small.any():
        mad = np.where(small, np.nanmean(absdev, axis=1) / 0.8, mad)
    good = finite & (absdev <= (cut * mad)[:, None])
    ngood = good.sum(axis=1)
    mean = np.where(
        ngood > 0, np.nansum(np.where(good, y2, 0.0), axis=1) / np.maximum(ngood, 1), np.nan
    )
    dev2 = np.where(good, (y2 - mean[:, None]) ** 2, 0.0)
    sigma = np.sqrt(dev2.sum(axis=1) / np.maximum(ngood, 1))
    sc = max(cut, 1.0)
    if sc <= 4.5:
        sigma = sigma / _poly(sc, _SIGCOEFF)
    good = finite & (absdev <= (cut * sigma)[:, None])
    ngood = good.sum(axis=1)
    mean = np.where(
        ngood > 0, np.nansum(np.where(good, y2, 0.0), axis=1) / np.maximum(ngood, 1), np.nan
    )
    dev2 = np.where(good, (y2 - mean[:, None]) ** 2, 0.0)
    sigma = np.sqrt(dev2.sum(axis=1) / np.maximum(ngood, 1))
    if sc <= 4.5:
        sigma = sigma / _poly(sc, _SIGCOEFF)
    sigma = sigma / np.sqrt(np.maximum(ngood - 1, 1))
    return mean, sigma, npts - ngood


def robust_sigma(y: np.ndarray, zero: bool = False) -> float:
    """Astrolib ``robust_sigma`` (Tukey biweight); NaNs ignored; returns -1 if <3 points."""
    y = np.asarray(y, dtype=float)
    y = y[np.isfinite(y)]
    eps = 1e-20
    if y.size == 0:
        return np.nan
    y0 = 0.0 if zero else np.median(y)
    mad = np.median(np.abs(y - y0)) / 0.6745
    if mad < eps:
        mad = np.mean(np.abs(y - y0)) / 0.8
    if mad < eps:
        return 0.0
    u = (y - y0) / (6.0 * mad)
    uu = u * u
    q = uu <= 1.0
    if q.sum() < 3:
        return -1.0
    n = y.size
    num = np.sum((y[q] - y0) ** 2 * (1 - uu[q]) ** 4)
    den1 = np.sum((1 - uu[q]) * (1 - 5 * uu[q]))
    sig2 = n * num / (den1 * (den1 - 1))
    return float(np.sqrt(sig2)) if sig2 > 0 else 0.0


def robust_poly_fit(x: np.ndarray, y: np.ndarray, degree: int, numit: int = 25) -> np.ndarray:
    """Astrolib ``robust_poly_fit``: iteratively re-weighted (biweight) polynomial fit.

    Returns coefficients lowest order first (IDL convention).
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    n = x.size
    if n <= degree + 1:
        raise ValueError("robust_poly_fit: too few points")
    coeff = np.polyfit(x, y, degree)[::-1]
    yfit = np.polyval(coeff[::-1], x)
    sig = robust_sigma(y - yfit)
    if sig <= 0:
        return coeff
    w = np.ones(n)
    for _ in range(numit):
        resid = y - yfit
        u = resid / (6.0 * max(sig, 1e-30))
        uu = u * u
        w = np.where(uu < 1, (1 - uu) ** 2, 0.0)
        if (w > 0).sum() <= degree + 1:
            break
        coeff = np.polyfit(x, y, degree, w=np.sqrt(w))[::-1]
        yfit_new = np.polyval(coeff[::-1], x)
        new_sig = robust_sigma(y - yfit_new)
        if new_sig <= 0:
            break
        conv = abs(new_sig - sig) / max(sig, 1e-30)
        sig = new_sig
        yfit = yfit_new
        if conv < 0.03:
            break
    return coeff


def combine_frames(
    stack: np.ndarray,
    method: str = "median",
    axis: int = 0,
    *,
    cut: float = 3.0,
    var: np.ndarray | None = None,
) -> np.ndarray:
    """Combine frames along ``axis``: ``median`` (NaN-aware, ``/EVEN``), ``mean``
    (resistant mean, 3σ), ``ivar`` (inverse-variance weighted, needs ``var``)."""
    stack = np.asarray(stack, dtype=float)
    if method == "median":
        return np.nanmedian(stack, axis=axis)
    if method in ("mean", "resistant", "meanadd"):
        return resistant_mean(stack, cut, axis=axis)
    if method == "ivar":
        if var is None:
            raise ValueError("ivar combine needs var")
        good = np.isfinite(stack) & np.isfinite(var) & (var > 0)
        w = np.where(good, 1.0 / np.where(good, var, 1.0), 0.0)
        num = np.sum(np.where(good, stack, 0.0) * w, axis=axis)
        den = np.sum(w, axis=axis)
        with np.errstate(invalid="ignore", divide="ignore"):
            out = num / den
        out[den == 0] = np.nan
        return out
    raise ValueError(f"unknown combine method {method!r}")
