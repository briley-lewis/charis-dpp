"""Robust statistics and frame combination.

Verbatim ports of the IDL Astronomy Library routines the pipeline relies on
(``resistant_mean``, ``robust_sigma``, ``robust_poly_fit``) so results match
the IDL reduction to round-off, plus ``combine_frames`` replacing
``charis_combinefits.pro`` / ``medfitsme``.
"""

from __future__ import annotations

import logging

import numpy as np

log = logging.getLogger(__name__)

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


def rob_checkfit(
    y: np.ndarray, yfit: np.ndarray, eps: float = 1e-20, del_: float = 5e-7, bfac: float = 6.0
):
    """Astrolib ``rob_checkfit``: ``(istat, sig, ngood, w)``; ``istat == 0`` means the fit is already good."""
    dev = y - yfit
    sig = robust_sigma(dev, zero=True)
    if sig < eps:
        return 0, sig, y.size, None
    if del_ > 0:
        q = np.abs(yfit) > eps
        fracdev = 0.0 if q.sum() < 3 else float(np.median(np.abs(dev[q] / yfit[q])))
        if fracdev < del_:
            return 0, sig, y.size, None
    b = np.abs(dev) / (bfac * sig)
    nbad = int((b > 1.0).sum())
    b = np.minimum(b, 1.0)
    w = 1.0 - b**2
    w = w / w.sum()
    return 1, sig, y.size - nbad, w


def robust_poly_fit(x: np.ndarray, y: np.ndarray, degree: int, numit: int = 25) -> np.ndarray:
    """Astrolib ``robust_poly_fit`` (the ``NPTS < 10000`` least-squares branch), ported line by line.

    Centres ``x``/``y`` on their means (``degree < 6``), starts from an unweighted
    fit, then iterates ``poly_fit(measure_errors = 1/W^2)`` (weights ``W^2`` in the
    least-squares sense, ``W = 1 - (|dev| / 6 sigma)^2`` from ``rob_checkfit``),
    dropping zero-weight points, until the robust sigma changes by less than
    ``max(0.03 sqrt(0.5 / (n - 1)), 5e-7)`` or ``numit`` iterations.  Returns the
    coefficients lowest order first, un-centred like the IDL.  The earlier port
    used biweight weights ``(1 - u^2)^2`` and a different stopping rule, which
    moved the fitted PSF centres by up to 0.01 px (BUGS.md #333).
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    npts = x.size
    minpts = degree + 1
    if npts <= degree + 1:
        raise ValueError("robust_poly_fit: too few points")
    degmax = 6
    if degree < degmax:
        x0, y0 = x.mean(), y.mean()
        u, v = x - x0, y - y0
    else:
        x0 = y0 = 0.0
        u, v = x, y

    def poly_fit(uu, vv, w=None):
        c = np.polyfit(uu, vv, degree, w=w)[::-1]
        return c, np.polyval(c[::-1], uu)

    badfit = False
    cc, yfit = poly_fit(u, v)
    istat, sig, ngood, w = rob_checkfit(v, yfit)
    if istat != 0:
        if ngood < minpts:
            log.warning("robust_poly_fit: no fit possible")
            return np.zeros(degree + 1)
        close_enough = max(0.03 * np.sqrt(0.5 / (npts - 1)), 5e-7)
        diff = 1e10
        sig_1 = min(100.0 * sig, 1e20)
        nit = 0
        while diff > close_enough and nit < numit:
            nit += 1
            sig_2 = sig_1
            sig_1 = sig
            g = w > 0
            if g.sum() < w.size:
                u, v, w = u[g], v[g], w[g]
            cc, yfit = poly_fit(u, v, w=w**2)
            istat, sig, ngood, w = rob_checkfit(v, yfit)
            if istat == 0:
                break
            if ngood < minpts:
                log.warning("robust_poly_fit: questionable fit")
                badfit = True
                break
            diff = min(abs(sig_1 - sig) / sig, abs(sig_2 - sig) / sig)
    cc = np.asarray(cc, dtype=float)
    if degree < degmax and degree >= 1:
        # un-centre: p(x) = sum c_k (x - x0)^k + y0
        shifted = np.polynomial.Polynomial(cc)(np.polynomial.Polynomial([-x0, 1.0]))
        cc = np.asarray(shifted.coef, dtype=float)
        if cc.size < degree + 1:
            cc = np.concatenate([cc, np.zeros(degree + 1 - cc.size)])
        cc[0] += y0
    if badfit:
        cc = np.concatenate([cc, [0.0]])
    return cc


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
