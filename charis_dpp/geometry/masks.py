"""Boolean masks replacing the IDL index-array idioms (``where``, ``get_cind``,
``indices_*.dat`` scratch files, ``charis_generate_roi``, ``mkspider``).

All masks are ``[y, x]`` boolean arrays; ``center`` is ``(x, y)``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from charis_dpp.geometry.grids import angle_grid, image_center, radial_grid

__all__ = [
    "disk",
    "annulus",
    "sectors",
    "sector_mask",
    "fov_mask",
    "point_mask",
    "roi_rect_rotated",
    "generate_roi_mask",
    "spider_mask",
    "Zone",
    "loci_zones",
]

XY = tuple[float, float]


def disk(shape, r: float, center: XY | None = None) -> np.ndarray:
    return radial_grid(shape, center) <= r


def annulus(shape, r_in: float, r_out: float, center: XY | None = None) -> np.ndarray:
    """``r_in <= r < r_out`` (half-open, so adjacent annuli tile the plane)."""
    r = radial_grid(shape, center)
    return (r >= r_in) & (r < r_out)


def point_mask(shape, x: float, y: float, r: float) -> np.ndarray:
    return radial_grid(shape, (x, y)) <= r


def sectors(shape, n: int, center: XY | None = None, phase_deg: float = 0.0) -> np.ndarray:
    """Integer label (0..n-1) of the azimuthal sector containing each pixel."""
    th = (angle_grid(shape, center, degrees=True) - phase_deg) % 360.0
    return np.minimum((th / (360.0 / n)).astype(int), n - 1)


def sector_mask(
    shape,
    r_in: float,
    r_out: float,
    n: int,
    k: int,
    center: XY | None = None,
    phase_deg: float = 0.0,
) -> np.ndarray:
    return annulus(shape, r_in, r_out, center) & (sectors(shape, n, center, phase_deg) == k)


def fov_mask(
    shape,
    half_size: float,
    angle_deg: float,
    center: XY | None = None,
    half_size_y: float | None = None,
) -> np.ndarray:
    """Rectangular field of view of half-widths ``half_size`` (x) and
    ``half_size_y`` (y, default ``half_size``) rotated by ``angle_deg`` (PDI FOV).

    ``charis_generate_roi, roirange=[dx, dy, pa]`` takes *full* widths, so a
    caller translating an IDL ``roirange`` passes ``dx / 2`` and ``dy / 2``.
    ``angle_deg`` follows the IDL ``rotat`` convention used there: the mask
    ``rotat(square, pa)`` produces is the square rotated by ``-pa`` in the
    (x right, y up) frame, so the pixel coordinates are rotated by ``-pa``
    before the box test (checked against ``charis_generate_roi`` output on
    the 2023-07-30 PDI flats: with the sign the other way round the mask
    stuck out of the lenslet field and the flat normalisation was off by
    1-6 %, BUGS.md #328).
    """
    from charis_dpp.geometry.grids import xy_grids

    if half_size_y is None:
        half_size_y = half_size
    dx, dy = xy_grids(shape, center)
    t = -np.deg2rad(angle_deg)
    u = dx * np.cos(t) + dy * np.sin(t)
    v = -dx * np.sin(t) + dy * np.cos(t)
    return (np.abs(u) <= half_size) & (np.abs(v) <= half_size_y)


def generate_roi_mask(
    shape, dx: float, dy: float, angle_deg: float, offset: XY = (0.0, 0.0)
) -> np.ndarray:
    """``charis_generate_roi, roirange=[dx, dy, pa]``, geometry part, done the way the
    IDL does it: an axis-aligned ``dx × dy`` box of ones (full widths, IDL's truncating
    float subscripts) centred on ``(nx//2 + offx, ny//2 + offy)``, rotated with
    ``rotat`` (:func:`rotate_image`, cubic with bilinear fallback) and thresholded at
    ``> 0``.  The analytic :func:`fov_mask` is ~500 px smaller on a 201 px CHARIS
    frame because the interpolated edge counts as inside here (BUGS.md #334).
    """
    from charis_dpp.image.interp import rotate_image

    ny, nx = shape[-2], shape[-1]
    offx, offy = int(offset[0]), int(offset[1])
    cx, cy = nx // 2 + offx, ny // 2 + offy
    hx, hy = int(dx / 2.0), int(dy / 2.0)  # IDL float subscripts truncate
    box = np.zeros((ny, nx), dtype=float)
    box[max(cy - hy, 0) : cy + hy + 1, max(cx - hx, 0) : cx + hx + 1] = 1.0
    if angle_deg == 0:
        return box == 1.0
    r = rotate_image(box, float(angle_deg), (float(cx), float(cy)))
    return np.isfinite(r) & (r > 0)


def roi_rect_rotated(
    shape, dx: float, dy: float, angle_deg: float, offset: XY = (0.0, 0.0), center: XY | None = None
) -> np.ndarray:
    """Rectangle of size (dx, dy) centred at centre+offset, rotated by ``angle_deg`` (``charis_generate_roi`` ``roirange``)."""
    from charis_dpp.geometry.grids import xy_grids

    xc, yc = image_center(shape) if center is None else center
    gx, gy = xy_grids(shape, (xc + offset[0], yc + offset[1]))
    t = -np.deg2rad(angle_deg)  # rotat convention, as in fov_mask (BUGS.md #328)
    u = gx * np.cos(t) + gy * np.sin(t)
    v = -gx * np.sin(t) + gy * np.cos(t)
    return (np.abs(u) <= dx / 2.0) & (np.abs(v) <= dy / 2.0)


def spider_mask(
    shape, width: float, angles_deg, center: XY | None = None, *, fractional: bool = True
) -> np.ndarray:
    """Spider-vane coverage: ``1`` outside the vanes, ``0`` on them (``mkspider`` returns ``1 - spider``).

    Each vane is a half-line from the centre at the given position angle
    (degrees from +x towards +y) of full width ``width`` pixels.  With
    ``fractional`` the edge pixels get their covered-area fraction (sub-pixel
    integration), otherwise a boolean mask is returned.  Fixes the doubled
    horizontal vane in ``mkspider.pro`` (BUGS.md #266).
    """
    from charis_dpp.geometry.grids import xy_grids

    dx, dy = xy_grids(shape, center)
    cover = np.zeros(shape[-2:], dtype=float)
    n_sub = 8 if fractional else 1
    offs = (np.arange(n_sub) + 0.5) / n_sub - 0.5
    for ang in np.atleast_1d(angles_deg):
        t = np.deg2rad(float(ang))
        ux, uy = np.cos(t), np.sin(t)  # along the vane
        frac = np.zeros_like(cover)
        for ox in offs:
            for oy in offs:
                px, py = dx + ox, dy + oy
                along = px * ux + py * uy
                perp = -px * uy + py * ux
                frac += ((along >= 0) & (np.abs(perp) <= width / 2.0)).astype(float)
        cover += frac / (n_sub * n_sub)
    cover = np.clip(cover, 0.0, 1.0)
    return (cover == 0) if not fractional else 1.0 - cover


@dataclass(frozen=True)
class Zone:
    """One LOCI/KLIP subtraction zone and its optimisation region."""

    ring: int
    sector: int
    r_in: float
    r_out: float
    sub: np.ndarray  # bool mask of the subtraction region
    opt: np.ndarray  # bool mask of the optimisation region
    n_sectors: int


def loci_zones(
    shape,
    rmin: float,
    rmax: float,
    dr: float | np.ndarray,
    fwhm: float,
    na: float,
    geom: float,
    center: XY | None = None,
    *,
    pixmask: bool = False,
    opt_full_ring: bool = False,
) -> list[Zone]:
    """Lafrenière-style zone geometry (ported from ``charis_adialoci.pro``).

    Subtraction annuli of width ``dr`` from ``rmin`` to ``rmax``; each is split
    into ``nt`` sectors so the optimisation region has area ``na·π(fwhm/2)²``
    with radial/azimuthal aspect ``geom``.  The optimisation region starts at
    the annulus inner radius and extends ``dropt`` outwards; with ``pixmask``
    the subtraction region itself is excluded from it.  ``opt_full_ring``
    uses the whole annulus (SDI variant).
    """
    xc, yc = image_center(shape) if center is None else center
    r = radial_grid(shape, (xc, yc))
    th = angle_grid(shape, (xc, yc), degrees=True) % 360.0
    area = na * np.pi * (fwhm / 2.0) ** 2
    widths = np.atleast_1d(np.asarray(dr, dtype=float))
    zones: list[Zone] = []
    edges = _ring_edges(rmin, rmax, widths)
    for ring, (ri, ro) in enumerate(edges):
        dropt = np.sqrt(geom * area)
        r2opt = min(ri + dropt, r.max())
        dropt_eff = max(r2opt - ri, 1.0)
        nt = max(int(round((2 * np.pi * (ri + dropt_eff / 2.0) * dropt_eff) / area)), 1)
        for it in range(nt):
            lo = it * 360.0 / nt
            hi = (it + 1) * 360.0 / nt
            sec = (th >= lo) & (th < hi)
            sub = (r >= ri) & (r < ro) & sec
            if opt_full_ring:
                opt = (r >= ri) & (r < ro)
            else:
                opt = (r >= ri) & (r < r2opt) & sec
            if pixmask:
                opt = opt & ~sub
            if sub.any():
                zones.append(Zone(ring, it, ri, ro, sub, opt, nt))
    return zones


def _ring_edges(rmin: float, rmax: float, widths: np.ndarray) -> list[tuple[float, float]]:
    if widths.size == 1:
        w = float(widths[0])
        n = max(int(np.ceil((rmax - rmin) / w)), 1)
        return [(rmin + i * w, min(rmin + (i + 1) * w, rmax)) for i in range(n)]
    # variable widths: IDL semantics drsub=[w_inner, r_break, r_break2?, w_outer] are
    # ambiguous; accept an explicit list of ring widths from rmin outwards.
    edges = []
    r = rmin
    for w in widths:
        if r >= rmax:
            break
        edges.append((r, min(r + float(w), rmax)))
        r += float(w)
    while r < rmax:
        edges.append((r, min(r + float(widths[-1]), rmax)))
        r += float(widths[-1])
    return edges
