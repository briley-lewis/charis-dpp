"""Resample a disk model (or PSF) image to the CHARIS pixel scale
(``charis_diskmodel_resample.pro``, ``charis_resamplepsf.pro``).

Fixes vs. IDL: no syntax error (#262); models larger than the output are
cropped after resampling rather than compared against a fixed 201-px half
size (#263); the PSF helper does not write ``tmp.fits`` and documents that
``scalefac > 1`` magnifies (#261).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from charis_dpp.image.interp import rescale_radial, subarr
from charis_dpp.io.fits import read_image

__all__ = ["resample_disk_model", "resample_psf"]


def resample_disk_model(
    model: np.ndarray | str | Path,
    input_pixscale: float,
    output_pixscale: float,
    *,
    out_size: int | None = None,
) -> np.ndarray:
    """Rescale a north-up model from ``input_pixscale`` to ``output_pixscale`` (arcsec/px),
    conserving surface brightness, and crop/pad to ``out_size`` pixels (default: input size)."""
    im = read_image(model) if isinstance(model, str | Path) else np.asarray(model, dtype=float)
    if im.ndim == 3:
        return np.stack(
            [resample_disk_model(s, input_pixscale, output_pixscale, out_size=out_size) for s in im]
        )
    size = int(out_size or im.shape[0])
    big = max(size, im.shape[0], im.shape[1])
    work = subarr(im, big, fill=0.0) if big != im.shape[0] or big != im.shape[1] else im
    scale = input_pixscale / output_pixscale  # >1: the model gets bigger on the CHARIS grid
    out = rescale_radial(work, scale, missing=0.0)
    return subarr(out, size, fill=0.0)


def resample_psf(psf: np.ndarray, scalefac: float) -> np.ndarray:
    """Magnify (``scalefac > 1``) or shrink a PSF stamp about its centre; same size out."""
    psf = np.asarray(psf, dtype=float)
    if psf.ndim == 3:
        return np.stack([resample_psf(s, scalefac) for s in psf])
    out = rescale_radial(psf, scalefac, missing=0.0)
    return np.where(np.isfinite(out), out, 0.0)
