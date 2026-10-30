"""Frame classification (``charis_check_health.pro``): ``bad`` / ``science`` / ``sky``."""

from __future__ import annotations

import numpy as np

from charis_dpp.config.instrument import Instrument
from charis_dpp.geometry.grids import image_center, radial_grid

__all__ = ["classify_frame"]


def classify_frame(
    cube: np.ndarray,
    inst: Instrument | None = None,
    *,
    pdi: bool = False,
    center: tuple[float, float] | None = None,
    ratio: float | None = None,
) -> str:
    """Classify a cube as ``'bad'`` (NaN-riddled core), ``'science'`` (bright centre) or ``'sky'``.

    Radii come from the instrument table as fractions of the frame width. For
    PDI data the two Wollaston beams are tested about their own centres and
    only the left half of the frame (as defined by the split line) is used.
    """
    cube = np.asarray(cube, dtype=float)
    avg = np.nanmedian(cube, axis=0) if cube.ndim == 3 else cube
    ny, nx = avg.shape
    h = inst.values.health if inst is not None else None
    bad_r = (float(h.bad_radius_frac) if h else 0.30) * nx
    inner_r = (float(h.inner_radius_frac) if h else 0.10) * nx
    outer_r = (float(h.outer_radius_frac) if h else 0.30) * nx
    ratio = float(h.science_ratio) if (h and ratio is None) else (ratio or 3.0)
    max_bad = int(h.max_bad_pixels) if h else 5

    r = radial_grid(avg.shape, center)
    nbad = int(np.sum((r <= bad_r) & ~np.isfinite(avg)))
    if nbad > max_bad:
        return "bad"
    finite = np.isfinite(avg) & (avg != 0)
    if not pdi:
        inner = finite & (r < inner_r)
        outer = finite & (r > outer_r)
    else:
        xc, yc = image_center(avg.shape) if center is None else center
        offs = inst.pdi.beam_centre_offsets_px if inst is not None else ((-30, -16), (30, 16))
        ox, oy = offs[0]
        scale = nx / float(inst.detector.cube_size_px) if inst is not None else 1.0
        r2 = radial_grid(avg.shape, (xc + ox * scale, yc + oy * scale))
        yy, xx = np.mgrid[0:ny, 0:nx]
        sl = inst.pdi.split_line if inst is not None else {"slope": -1.887, "intercept_px": 286.8}
        half = yy < float(sl["slope"]) * xx + float(sl["intercept_px"]) * scale
        inner = finite & (r2 < inner_r) & half
        outer = finite & (r2 > outer_r) & half
    if inner.sum() == 0 or outer.sum() == 0:
        return "sky"
    return "science" if np.mean(avg[inner]) > ratio * np.mean(avg[outer]) else "sky"
