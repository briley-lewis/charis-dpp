"""Sub-pixel interpolation, shifting, rotation and rescaling.

IDL's ``interpolate(..., cubic=-0.5)`` is the Keys cubic-convolution kernel
with ``a = -0.5``; scipy's ``map_coordinates(order=3)`` is a B-spline
interpolant and differs at the 1e-3 level, so the kernel is implemented here
(:func:`keys_cubic`).  :func:`shift_sub` and :func:`rotate_image` port
``subroutines/shift_sub.pro`` and ``subroutines/rotat.pro`` (including the
latter's NaN footprint logic and bilinear fallback).
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage

from charis_dpp.geometry.grids import image_center

__all__ = [
    "keys_cubic",
    "bilinear",
    "shift_sub",
    "shift_cube",
    "rotate_image",
    "rotate_cube",
    "rescale_radial",
    "subarr",
]


def _keys_weights(t: np.ndarray, a: float) -> np.ndarray:
    """Weights for offsets (-1, 0, 1, 2) relative to floor(x), t = frac(x), shape (..., 4)."""
    t = np.asarray(t, dtype=float)
    d = np.stack([t + 1.0, t, 1.0 - t, 2.0 - t], axis=-1)  # |x - x_k| for k = -1..2
    ad = np.abs(d)
    w_near = (a + 2.0) * ad**3 - (a + 3.0) * ad**2 + 1.0
    w_far = a * ad**3 - 5.0 * a * ad**2 + 8.0 * a * ad - 4.0 * a
    return np.where(ad < 1.0, w_near, np.where(ad < 2.0, w_far, 0.0))


def keys_cubic(
    im: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    a: float = -0.5,
    missing: float = np.nan,
) -> np.ndarray:
    """Sample ``im[y, x]`` at fractional ``(x, y)`` with the Keys cubic kernel.

    Points outside ``[0, nx-1] × [0, ny-1]`` get ``missing``; neighbours that
    fall outside the frame are clamped to the edge (IDL behaviour).  NaN pixels
    propagate to every output sample whose 4×4 stencil touches them.
    """
    im = np.asarray(im, dtype=float)
    ny, nx = im.shape
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    out = np.full(np.broadcast(x, y).shape, missing, dtype=float)
    inside = (x >= 0) & (x <= nx - 1) & (y >= 0) & (y <= ny - 1)
    if not inside.any():
        return out
    xs = np.broadcast_to(x, inside.shape)[inside]
    ys = np.broadcast_to(y, inside.shape)[inside]
    x0 = np.floor(xs).astype(int)
    y0 = np.floor(ys).astype(int)
    wx = _keys_weights(xs - x0, a)  # (n, 4)
    wy = _keys_weights(ys - y0, a)
    acc = np.zeros(xs.shape, dtype=float)
    for j in range(4):
        yy = np.clip(y0 + j - 1, 0, ny - 1)
        row = np.zeros(xs.shape, dtype=float)
        for i in range(4):
            xx = np.clip(x0 + i - 1, 0, nx - 1)
            row += wx[:, i] * im[yy, xx]
        acc += wy[:, j] * row
    out[inside] = acc
    return out


def bilinear(im: np.ndarray, x: np.ndarray, y: np.ndarray, missing: float = np.nan) -> np.ndarray:
    """Bilinear sampling with the same edge/missing conventions as :func:`keys_cubic`."""
    im = np.asarray(im, dtype=float)
    ny, nx = im.shape
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    out = np.full(np.broadcast(x, y).shape, missing, dtype=float)
    inside = (x >= 0) & (x <= nx - 1) & (y >= 0) & (y <= ny - 1)
    if not inside.any():
        return out
    xs = np.broadcast_to(x, inside.shape)[inside]
    ys = np.broadcast_to(y, inside.shape)[inside]
    x0 = np.clip(np.floor(xs).astype(int), 0, nx - 1)
    y0 = np.clip(np.floor(ys).astype(int), 0, ny - 1)
    x1 = np.clip(x0 + 1, 0, nx - 1)
    y1 = np.clip(y0 + 1, 0, ny - 1)
    tx = xs - x0
    ty = ys - y0
    v = (
        im[y0, x0] * (1 - tx) * (1 - ty)
        + im[y0, x1] * tx * (1 - ty)
        + im[y1, x0] * (1 - tx) * ty
        + im[y1, x1] * tx * ty
    )
    out[inside] = v
    return out


def shift_sub(im: np.ndarray, dx: float, dy: float, *, edge: str = "nan") -> np.ndarray:
    """Shift an image by ``(dx, dy)`` pixels (positive moves content to +x/+y).

    ``edge="nan"``: pixels shifted in from outside are NaN for both integer
    and fractional shifts.  ``edge="idl"`` reproduces ``shift_sub.pro``:
    integer shifts wrap around and fractional shifts clamp to the edge value.
    """
    im = np.asarray(im, dtype=float)
    ny, nx = im.shape
    if float(dx).is_integer() and float(dy).is_integer():
        ix, iy = int(dx), int(dy)
        if edge == "idl":
            return np.roll(np.roll(im, iy, axis=0), ix, axis=1)
        out = np.full_like(im, np.nan)
        ys, yd = (
            (slice(0, ny - iy), slice(iy, ny)) if iy >= 0 else (slice(-iy, ny), slice(0, ny + iy))
        )
        xs, xd = (
            (slice(0, nx - ix), slice(ix, nx)) if ix >= 0 else (slice(-ix, nx), slice(0, nx + ix))
        )
        out[yd, xd] = im[ys, xs]
        return out
    y, x = np.mgrid[0:ny, 0:nx].astype(float)
    xs = x - dx
    ys = y - dy
    if edge == "idl":
        xs = np.clip(xs, 0, nx - 1)
        ys = np.clip(ys, 0, ny - 1)
    return keys_cubic(im, xs, ys)


def shift_cube(cube: np.ndarray, dx: float, dy: float, **kw) -> np.ndarray:
    return np.stack([shift_sub(sl, dx, dy, **kw) for sl in cube])


def _rotation_maps(shape, deg, x0, y0, finite_mask):
    """Port of rotat.pro: forward-map finite pixels, dilate, inverse-map (xi, yi) for those pixels."""
    ny, nx = shape
    theta = np.deg2rad(deg)
    ct, st = np.cos(theta), np.sin(theta)
    yi, xi = np.nonzero(finite_mask)
    if xi.size == 0:
        return None
    xf = np.rint((x0 - ct * x0 - st * y0) + ct * xi + st * yi).astype(int)
    yf = np.rint((y0 + st * x0 - ct * y0) - st * xi + ct * yi).astype(int)
    ok = (xf >= 0) & (xf < nx) & (yf >= 0) & (yf < ny)
    mask = np.zeros((ny, nx), dtype=bool)
    mask[yf[ok], xf[ok]] = True
    mask = ndimage.binary_dilation(mask, structure=np.ones((3, 3), dtype=bool))
    yf, xf = np.nonzero(mask)
    xin = (x0 - ct * x0 + st * y0) + ct * xf - st * yf
    yin = (y0 - st * x0 - ct * y0) + st * xf + ct * yf
    return mask, xin, yin


def rotate_image(
    im: np.ndarray,
    deg: float,
    center: tuple[float, float] | None = None,
    *,
    missing: float = np.nan,
) -> np.ndarray:
    """Rotate ``im`` by ``deg`` about ``center`` (default frame centre), IDL ``rotat`` conventions.

    Positive ``deg`` rotates image content clockwise on a y-up display, exactly
    as ``rotat.pro`` does.  Only pixels that receive finite input are computed;
    where the cubic kernel yields NaN (stencil touching NaN), the bilinear
    value is used instead (as in ``rotat.pro``).
    """
    im = np.asarray(im, dtype=float)
    x0, y0 = image_center(im.shape) if center is None else center
    out = np.full_like(im, np.nan)
    maps = _rotation_maps(im.shape, deg, x0, y0, np.isfinite(im))
    if maps is None:
        return out
    mask, xin, yin = maps
    vals = keys_cubic(im, xin, yin, missing=missing)
    bad = ~np.isfinite(vals)
    if bad.any():
        vals[bad] = bilinear(im, xin[bad], yin[bad], missing=missing)
    out[mask] = vals
    return out


def rotate_cube(
    cube: np.ndarray,
    deg: float,
    center: tuple[float, float] | None = None,
    *,
    missing: float = np.nan,
) -> np.ndarray:
    """Rotate every slice of ``(nλ, ny, nx)`` with one footprint (fixes rotat_cube's ``yi mod dimy`` bug)."""
    cube = np.asarray(cube, dtype=float)
    x0, y0 = image_center(cube.shape[-2:]) if center is None else center
    out = np.full_like(cube, np.nan)
    maps = _rotation_maps(cube.shape[-2:], deg, x0, y0, np.isfinite(cube).any(axis=0))
    if maps is None:
        return out
    mask, xin, yin = maps
    for l in range(cube.shape[0]):
        vals = keys_cubic(cube[l], xin, yin, missing=missing)
        bad = ~np.isfinite(vals)
        if bad.any():
            vals[bad] = bilinear(cube[l], xin[bad], yin[bad], missing=missing)
        out[l][mask] = vals
    return out


def rescale_radial(
    im: np.ndarray,
    scale: float,
    center: tuple[float, float] | None = None,
    *,
    missing: float = np.nan,
) -> np.ndarray:
    """Magnify ``im`` about ``center`` by ``scale`` (>1 enlarges) with cubic interpolation.

    Used to align speckles between wavelengths: a slice at λ is rescaled by
    ``λ_ref / λ`` to match the reference slice.
    """
    im = np.asarray(im, dtype=float)
    ny, nx = im.shape
    x0, y0 = image_center(im.shape) if center is None else center
    y, x = np.mgrid[0:ny, 0:nx].astype(float)
    xs = x0 + (x - x0) / scale
    ys = y0 + (y - y0) / scale
    return keys_cubic(im, xs, ys, missing=missing)


def subarr(
    im: np.ndarray,
    size: int | tuple[int, int],
    center: tuple[float, float] | None = None,
    fill=np.nan,
) -> np.ndarray:
    """Extract (or pad to) a ``size`` box centred on ``center`` (integer centre), IDL ``subarr``."""
    im = np.asarray(im)
    ny, nx = im.shape[-2:]
    sy, sx = (size, size) if isinstance(size, int) else size
    xc, yc = (
        image_center(im.shape) if center is None else (int(round(center[0])), int(round(center[1])))
    )
    out_shape = im.shape[:-2] + (sy, sx)
    out = np.full(out_shape, fill, dtype=float)
    x0, y0 = xc - sx // 2, yc - sy // 2
    xs0, ys0 = max(0, x0), max(0, y0)
    xs1, ys1 = min(nx, x0 + sx), min(ny, y0 + sy)
    if xs1 <= xs0 or ys1 <= ys0:
        return out
    out[..., ys0 - y0 : ys1 - y0, xs0 - x0 : xs1 - x0] = im[..., ys0:ys1, xs0:xs1]
    return out
