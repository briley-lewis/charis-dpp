"""Coordinate grids.

One pixel-centre convention for the whole package: :func:`image_center`
returns ``(x, y) = (nx // 2, ny // 2)``.  For the odd-sized (201×201) CHARIS
cubes this equals both the IDL ``dim/2`` and ``(dim-1)/2`` conventions that
the original code mixed; for even sizes it differs from IDL ``dist_circle``
by half a pixel (documented in BUGS.md).

Arrays are indexed ``[y, x]``; all public functions take ``shape=(ny, nx)``
and return ``(x, y)`` ordered tuples or ``[y, x]`` arrays.
"""

from __future__ import annotations

import numpy as np

Shape = tuple[int, int]
XY = tuple[float, float]


def image_center(shape: Shape | np.ndarray) -> XY:
    """``(xc, yc)`` of a frame; accepts a shape tuple or an array."""
    if isinstance(shape, np.ndarray):
        shape = shape.shape[-2:]
    ny, nx = int(shape[-2]), int(shape[-1])
    return (nx // 2, ny // 2)


def xy_grids(shape: Shape, center: XY | None = None) -> tuple[np.ndarray, np.ndarray]:
    """``(dx, dy)`` pixel offsets from ``center`` (default frame centre), each ``[y, x]``."""
    ny, nx = int(shape[-2]), int(shape[-1])
    xc, yc = image_center(shape) if center is None else center
    y, x = np.mgrid[0:ny, 0:nx]
    return x - xc, y - yc


def radial_grid(shape: Shape, center: XY | None = None) -> np.ndarray:
    """Distance of each pixel from ``center`` (IDL ``dist_circle``)."""
    dx, dy = xy_grids(shape, center)
    return np.hypot(dx, dy)


def angle_grid(shape: Shape, center: XY | None = None, degrees: bool = False) -> np.ndarray:
    """Position angle of each pixel measured from +x towards +y (IDL ``angarr``), in (-π, π]."""
    dx, dy = xy_grids(shape, center)
    ang = np.arctan2(dy, dx)
    return np.degrees(ang) if degrees else ang


def polar_grids(shape: Shape, center: XY | None = None) -> tuple[np.ndarray, np.ndarray]:
    """``(r, theta)`` with theta in radians."""
    dx, dy = xy_grids(shape, center)
    return np.hypot(dx, dy), np.arctan2(dy, dx)
