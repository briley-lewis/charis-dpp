"""Centroiding.

* :func:`cntrd` — astrolib ``cntrd`` (derivative-weighted marginal centroid).
* :func:`gcntrd` — Gaussian centroid: marginal profiles fit with a fixed-σ
  Gaussian plus constant (same inputs/outputs as astrolib ``gcntrd``; failure
  is signalled with NaN instead of ``-1``).
* :func:`fit_gaussian2d` — elliptical 2-D Gaussian fit (replaces the missing
  ``calc_centroid_spots4`` used by ``register_cube /fitgauss``).
* :func:`center_of_symmetry` — ``charis_getcenter.pro`` (180° symmetry search).
* :func:`xcorr_shift` — ``charis_getrelshift.pro`` (cross-correlation shift).
"""

from __future__ import annotations

import numpy as np
from astropy.modeling import fitting, models
from scipy.optimize import least_squares

from charis_dpp.geometry.grids import image_center
from charis_dpp.image.interp import keys_cubic, rotate_image

__all__ = ["cntrd", "gcntrd", "fit_gaussian2d", "center_of_symmetry", "xcorr_shift"]

NAN2 = (np.nan, np.nan)


def _prep(im: np.ndarray, x: float, y: float, fwhm: float, keepcenter: bool):
    im = np.asarray(im, dtype=float)
    ny, nx = im.shape
    nhalf = max(int(0.637 * fwhm), 2)
    nbox = 2 * nhalf + 1
    ix, iy = int(round(x)), int(round(y))
    if not keepcenter:
        x0, x1 = max(ix - nhalf, 0), min(ix + nhalf, nx - 1)
        y0, y1 = max(iy - nhalf, 0), min(iy + nhalf, ny - 1)
        sub = im[y0 : y1 + 1, x0 : x1 + 1]
        if sub.size == 0 or not np.isfinite(sub).any():
            return None
        jy, jx = np.unravel_index(np.nanargmax(sub), sub.shape)
        ix, iy = x0 + int(jx), y0 + int(jy)
    if ix - nhalf < 0 or ix + nhalf > nx - 1 or iy - nhalf < 0 or iy + nhalf > ny - 1:
        return None
    d = im[iy - nhalf : iy + nhalf + 1, ix - nhalf : ix + nhalf + 1]
    if not np.isfinite(d).all():
        return None
    return d, ix, iy, nhalf, nbox


def cntrd(
    im: np.ndarray, x: float, y: float, fwhm: float, *, keepcenter: bool = False
) -> tuple[float, float]:
    """Astrolib ``cntrd``; returns ``(xcen, ycen)`` or ``(nan, nan)`` on failure."""
    p = _prep(im, x, y, fwhm, keepcenter)
    if p is None:
        return NAN2
    d, ix, iy, nhalf, nbox = p
    dd = np.arange(nbox - 1) + 0.5 - nhalf
    w = 1.0 - 0.5 * (np.abs(dd) - 0.5) / (nhalf - 0.5)
    sumc = w.sum()

    def axis_centroid(prof_deriv: np.ndarray) -> float | None:
        sumd = np.sum(w * prof_deriv)
        sumxd = np.sum(w * dd * prof_deriv)
        sumxsq = np.sum(w * dd**2)
        sumx = np.sum(w * dd)
        den = sumc * sumxd - sumd * sumx
        if den == 0 or sumxd >= 0:
            return None
        dx = sumxsq * sumd / den
        if abs(dx) > nhalf:
            return None
        return dx

    # x: derivative along x summed over y; y: derivative along y summed over x
    derx = (d[:, 1:] - d[:, :-1]).sum(axis=0)
    dery = (d[1:, :] - d[:-1, :]).sum(axis=1)
    dx = axis_centroid(derx)
    dy = axis_centroid(dery)
    if dx is None or dy is None:
        return NAN2
    return ix - dx, iy - dy


def gcntrd(
    im: np.ndarray, x: float, y: float, fwhm: float, *, keepcenter: bool = False
) -> tuple[float, float]:
    """Gaussian centroid of the marginal profiles (astrolib ``gcntrd`` semantics)."""
    p = _prep(im, x, y, fwhm, keepcenter)
    if p is None:
        return NAN2
    d, ix, iy, nhalf, nbox = p
    sig = fwhm / 2.35482
    pos = np.arange(nbox) - nhalf

    def fit1d(prof: np.ndarray) -> float | None:
        amp0 = prof.max() - prof.min()
        if amp0 <= 0:
            return None
        c0 = float(pos[np.argmax(prof)])

        def resid(p):
            a, c, b = p
            return a * np.exp(-0.5 * ((pos - c) / sig) ** 2) + b - prof

        res = least_squares(resid, [amp0, c0, prof.min()], method="lm")
        a, c, _ = res.x
        if not res.success or a <= 0 or abs(c) > nhalf:
            return None
        return float(c)

    cx = fit1d(d.sum(axis=0))
    cy = fit1d(d.sum(axis=1))
    if cx is None or cy is None:
        return NAN2
    return ix + cx, iy + cy


def fit_gaussian2d(
    im: np.ndarray,
    x: float,
    y: float,
    fwhm: float,
    box: int | None = None,
) -> dict[str, float]:
    """Elliptical 2-D Gaussian + constant fit around ``(x, y)``.

    Returns ``{"x", "y", "fwhm_x", "fwhm_y", "theta_deg", "amplitude", "background"}``;
    ``x``/``y`` are NaN on failure.
    """
    im = np.asarray(im, dtype=float)
    ny, nx = im.shape
    half = int(round(1.5 * fwhm)) if box is None else box // 2
    ix, iy = int(round(x)), int(round(y))
    x0, x1 = max(ix - half, 0), min(ix + half, nx - 1)
    y0, y1 = max(iy - half, 0), min(iy + half, ny - 1)
    sub = im[y0 : y1 + 1, x0 : x1 + 1]
    good = np.isfinite(sub)
    fail = {
        "x": np.nan,
        "y": np.nan,
        "fwhm_x": np.nan,
        "fwhm_y": np.nan,
        "theta_deg": np.nan,
        "amplitude": np.nan,
        "background": np.nan,
    }
    if good.sum() < 7:
        return fail
    yy, xx = np.mgrid[y0 : y1 + 1, x0 : x1 + 1]
    bkg = np.nanmedian(sub)
    amp = np.nanmax(sub) - bkg
    sig = fwhm / 2.35482
    g = models.Gaussian2D(
        amplitude=amp, x_mean=x, y_mean=y, x_stddev=sig, y_stddev=sig, theta=0.0
    ) + models.Const2D(bkg)
    g.x_stddev_0.bounds = (0.3, 5 * sig)
    g.y_stddev_0.bounds = (0.3, 5 * sig)
    fitter = fitting.LevMarLSQFitter()
    try:
        fit = fitter(g, xx[good], yy[good], sub[good], maxiter=200)
    except Exception:
        return fail
    xf, yf = float(fit.x_mean_0.value), float(fit.y_mean_0.value)
    if not (x0 <= xf <= x1 and y0 <= yf <= y1) or fit.amplitude_0.value <= 0:
        return fail
    return {
        "x": xf,
        "y": yf,
        "fwhm_x": 2.35482 * float(fit.x_stddev_0.value),
        "fwhm_y": 2.35482 * float(fit.y_stddev_0.value),
        "theta_deg": float(np.degrees(fit.theta_0.value)),
        "amplitude": float(fit.amplitude_0.value),
        "background": float(fit.amplitude_1.value),
    }


def _norm(v: np.ndarray) -> np.ndarray:
    v = v - np.nanmean(v)
    n = np.sqrt(np.nansum(v**2))
    return v / n if n > 0 else v


def center_of_symmetry(
    im: np.ndarray,
    x_guess: float,
    y_guess: float,
    roi: np.ndarray,
    *,
    smallsteps: bool = False,
) -> tuple[float, float]:
    """Star centre from 180° rotational symmetry (``charis_getcenter.pro``).

    ``roi`` is a boolean mask of pixels to compare.  The image and its 180°
    rotation are shifted in opposite directions over a grid of offsets; the
    offset minimising Σ|im − rot(im)| wins, refined down to 0.01 px.
    """
    im = np.asarray(im, dtype=float)
    imr = rotate_image(im, 180.0)
    xc, yc = image_center(im.shape)
    yroi, xroi = np.nonzero(roi)
    sx0 = round(xc - x_guess)
    sy0 = round(yc - y_guess)

    def score(sx: float, sy: float) -> float:
        a = _norm(keys_cubic(im, xroi - sx, yroi - sy, missing=0.0))
        b = _norm(keys_cubic(imr, xroi + sx, yroi + sy, missing=0.0))
        return -float(np.nansum(np.abs(a - b)))

    def search(sx_c: float, sy_c: float, ntry: int, step: float) -> tuple[float, float]:
        offs = (np.arange(ntry) - ntry // 2) * step
        best, bx, by = -np.inf, sx_c, sy_c
        for ox in offs:
            for oy in offs:
                s = score(sx_c + ox, sy_c + oy)
                if s > best:
                    best, bx, by = s, sx_c + ox, sy_c + oy
        return bx, by

    sx0, sy0 = search(sx0, sy0, 21, 0.5 if smallsteps else 2.0)
    for step in (0.5, 0.2, 0.1, 0.05, 0.01):
        sx0, sy0 = search(sx0, sy0, 11, step)
    return xc - sx0, yc - sy0


def xcorr_shift(
    im1: np.ndarray,
    im2: np.ndarray,
    roi: np.ndarray,
    sx_guess: float = 0.0,
    sy_guess: float = 0.0,
    *,
    subpix_only: bool = False,
    residual: bool = False,
) -> tuple[float, float]:
    """Shift that aligns ``im2`` onto ``im1`` within ``roi`` (``charis_getrelshift.pro``).

    Returns ``(x, y)`` in the IDL convention: the position in ``im2`` of the
    feature that is at the frame centre in ``im1`` — i.e. ``center - (sx, sy)``.
    Coarse integer search (unless ``subpix_only``) then refinement to 0.001 px.
    """
    im1 = np.asarray(im1, dtype=float)
    im2 = np.asarray(im2, dtype=float)
    xc, yc = image_center(im1.shape)
    yroi, xroi = np.nonzero(roi)
    ref = _norm(im1[yroi, xroi])
    sx0, sy0 = (
        (float(sx_guess), float(sy_guess)) if subpix_only else (round(sx_guess), round(sy_guess))
    )

    def score(sx: float, sy: float) -> float:
        t = _norm(keys_cubic(im2, xroi - sx, yroi - sy, missing=0.0))
        if residual:
            return -float(np.nansum(np.abs(ref - t)))
        return float(np.nansum(ref * t))

    def search(cx, cy, ntry, step):
        offs = (np.arange(ntry) - ntry // 2) * step
        best, bx, by = -np.inf, cx, cy
        for ox in offs:
            for oy in offs:
                s = score(cx + ox, cy + oy)
                if s > best:
                    best, bx, by = s, cx + ox, cy + oy
        return bx, by

    if not subpix_only:
        sx0, sy0 = search(sx0, sy0, 21, 1.0)
    for k in (1, 2, 3):
        sx0, sy0 = search(sx0, sy0, 21, 10.0**-k)
    return xc - sx0, yc - sy0
