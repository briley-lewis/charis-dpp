"""Inject a model planet into a cube (``charis_insert_planet_into_cube.pro``).

The planet is a PSF stamp (normalised to unit flux inside ``R_AP<l>``) times
the model spectrum in the cube's flux units, smeared along the arc the
field rotates through during the exposure (``nsc`` sub-companions).

Position convention (one for the whole package): planet coordinates are
given in the **final north-up frame**; :func:`frame_position` converts them
to a given frame's detector coordinates using that frame's north-up
rotation (the inverse of the ``rotat`` forward map applied by ``northup``).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from astropy.io.fits import Header

from charis_dpp.astro.time_coords import parallactic_angle, ten
from charis_dpp.geometry.apertures import aperture_sum
from charis_dpp.image.interp import shift_sub

__all__ = ["Planet", "frame_position", "normalise_psf", "insert_planet", "paste_stamp"]


@dataclass
class Planet:
    x: float  # north-up frame, 0-based pixels
    y: float
    log_contrast: float = -4.0
    spectrum: np.ndarray | None = None  # per-channel flux in the cube's units (overrides contrast)


def frame_position(
    x: float, y: float, rot_deg: float, center: tuple[float, float]
) -> tuple[float, float]:
    """North-up position → detector position of a frame that needs rotation ``rot_deg`` to be north-up."""
    t = np.deg2rad(rot_deg)
    ct, st = np.cos(t), np.sin(t)
    ux, uy = x - center[0], y - center[1]
    return center[0] + ct * ux - st * uy, center[1] + st * ux + ct * uy


def normalise_psf(psf: np.ndarray, r_ap: np.ndarray) -> np.ndarray:
    """Unit flux inside radius ``r_ap[l]`` about the stamp centre, per channel."""
    out = np.empty_like(psf, dtype=float)
    c = psf.shape[-1] // 2
    for l in range(psf.shape[0]):
        s = aperture_sum(np.nan_to_num(psf[l]), c, c, float(r_ap[l]))
        out[l] = np.nan_to_num(psf[l]) / s if s > 0 else 0.0
    return out


def paste_stamp(im: np.ndarray, stamp: np.ndarray, x: float, y: float, scale: float = 1.0) -> None:
    """Add ``stamp`` (odd size, centred) at sub-pixel ``(x, y)`` into ``im`` in place, clipping at the edges."""
    d = stamp.shape[0]
    h = d // 2
    xe, ye = int(np.floor(x)), int(np.floor(y))
    s = shift_sub(stamp, x - xe, y - ye)
    s = np.nan_to_num(s) * scale
    ny, nx = im.shape
    x0, x1 = xe - h, xe + h + 1
    y0, y1 = ye - h, ye + h + 1
    sx0, sy0 = max(0, -x0), max(0, -y0)
    sx1, sy1 = d - max(0, x1 - nx), d - max(0, y1 - ny)
    if sx1 <= sx0 or sy1 <= sy0:
        return
    im[max(y0, 0) : min(y1, ny), max(x0, 0) : min(x1, nx)] += s[sy0:sy1, sx0:sx1]


def smear_offsets(h0: Header, dec_deg: float, lat_deg: float, nsc: int) -> np.ndarray:
    """Parallactic-angle offsets (deg) of ``nsc`` sub-companions across the exposure, centred on zero."""
    ha = ten(h0.get("HA", 0.0))
    tint = float(h0.get("EXP1TIME", 1.0)) * float(h0.get("COADDS", 1))
    hours = ha + (np.arange(nsc) / max(nsc - 1, 1)) * (tint / 3600.0)
    pa = np.unwrap(np.deg2rad(parallactic_angle(hours, dec_deg, lat_deg)))
    dpa = np.degrees(pa - pa[0])
    return dpa - 0.5 * (dpa.max() + dpa.min())


def insert_planet(
    cube: np.ndarray,
    x: float,
    y: float,
    spectrum: np.ndarray,
    psf: np.ndarray,
    *,
    center: tuple[float, float],
    negative: bool = False,
    h0: Header | None = None,
    dec_deg: float | None = None,
    lat_deg: float | None = None,
    nsc: int = 21,
) -> np.ndarray:
    """Return ``cube`` plus a planet at detector position ``(x, y)`` with per-channel flux ``spectrum``.

    ``psf`` must already be normalised (:func:`normalise_psf`).  With ``h0``,
    ``dec_deg`` and ``lat_deg`` the planet is smeared over the exposure.
    """
    out = np.array(cube, dtype=float, copy=True)
    spec = -np.asarray(spectrum, dtype=float) if negative else np.asarray(spectrum, dtype=float)
    dx, dy = x - center[0], y - center[1]
    rc = np.hypot(dx, dy)
    th0 = np.arctan2(dy, dx)
    if h0 is not None and dec_deg is not None and lat_deg is not None and nsc > 1:
        dpa = np.deg2rad(smear_offsets(h0, dec_deg, lat_deg, nsc))
    else:
        dpa = np.zeros(1)
    xs = center[0] + rc * np.cos(th0 + dpa)
    ys = center[1] + rc * np.sin(th0 + dpa)
    n = len(dpa)
    for l in range(out.shape[0]):
        for xi, yi in zip(xs, ys, strict=True):
            paste_stamp(out[l], psf[l], xi, yi, spec[l] / n)
    return out
