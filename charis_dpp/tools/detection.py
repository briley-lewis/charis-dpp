"""Point-source detection in a PSF-subtracted collapsed image (``charis_check_for_detection.pro``).

Fixes vs. IDL: works for any frame size (#240); when the cube carries no
LOCI/KLIP radius cards the search covers the whole frame instead of nothing
(#239); the de-duplication radius and centroid guard against failures.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from charis_dpp.config.dataset import DatasetConfig
from charis_dpp.geometry.apertures import aper_idl
from charis_dpp.geometry.grids import image_center, radial_grid
from charis_dpp.image.centroid import cntrd
from charis_dpp.image.snr import snr_map
from charis_dpp.io.fits import read_cube, write_image
from charis_dpp.planetspec.planet_photcal import star_spectrum_from_header
from charis_dpp.prelims.setupdir import Dirs

log = logging.getLogger(__name__)

__all__ = ["Detection", "find_detections"]


@dataclass
class Detection:
    x: float
    y: float
    snr: float
    flux: float
    log_contrast: float


def find_detections(
    cfg: DatasetConfig | str | Path,
    datacube: str | Path = "final.fits",
    *,
    threshold: float | None = None,
    excl: float | None = None,
    rangeconstant: float | None = None,
    outfile: str | Path | None = None,
    snrmap_out: str | Path | None = None,
) -> list[Detection]:
    cfg = cfg if isinstance(cfg, DatasetConfig) else DatasetConfig.load(cfg)
    res = cfg.resolve()
    dirs = Dirs(cfg.workdir)
    det = cfg.reduction.get("detection", {})
    excl = float(excl if excl is not None else det.get("excl_px", 2.0))
    rangeconstant = float(
        rangeconstant if rangeconstant is not None else det.get("rangeconstant", 1.5)
    )
    path = (
        Path(datacube)
        if Path(datacube).is_absolute() or Path(datacube).exists()
        else dirs.proc / datacube
    )
    col_path = path.with_name(path.stem + "_collapsed.fits")
    cube = read_cube(path)
    im = read_cube(col_path).data if col_path.exists() else np.nanmedian(cube.data, axis=0)
    h1 = cube.h1
    fwhm_col = float(np.median(res.fwhm_px))
    center = res.star_center or image_center(im.shape)
    rmin = float(h1.get("LOCI_RMI", h1.get("KLIP_RMI", 0.0)))
    rmax = float(h1.get("LOCI_RMA", h1.get("KLIP_RMA", min(im.shape) / 2.0)))
    if rmax <= 0:
        rmax = min(im.shape) / 2.0
    if threshold is None:
        threshold = (
            float(det.get("threshold", 5.0))
            if ("ZERO" in h1 or "MEANADD" in h1)
            else 2 * float(det.get("threshold", 5.0))
        )
    r = snr_map(im, fwhm_col, int(min(im.shape) / 2), center=center, zero=True, finite=True)
    snr = r.snrmap
    if snrmap_out:
        write_image(snrmap_out, snr)
    g = radial_grid(im.shape, center)
    cand = (
        np.isfinite(snr)
        & (snr >= threshold)
        & (g >= rmin + rangeconstant * fwhm_col)
        & (g <= rmax - rangeconstant * fwhm_col)
    )
    detections: list[Detection] = []
    if not cand.any():
        log.info("no pixels above SNR %.1f", threshold)
        if outfile:
            Path(outfile).write_text("")
        return detections
    detectmap = cand.copy()
    order = np.argsort(snr[cand])[::-1]
    ys, xs = np.nonzero(cand)
    ys, xs = ys[order], xs[order]
    for yi, xi in zip(ys, xs, strict=True):
        if not detectmap[yi, xi]:
            continue
        rad = excl + 1.5 * np.log10(max(snr[yi, xi] / 3.0, 1.0001)) * fwhm_col
        rr = radial_grid(im.shape, (xi, yi))
        kill = (rr > 1e-9) & (rr <= rad) & detectmap
        detectmap[kill] = False
    starflux = star_spectrum_from_header(h1, cube.nlam)
    med_star = float(np.nanmedian(starflux)) if np.isfinite(starflux).any() else np.nan
    for yi, xi in zip(*np.nonzero(detectmap), strict=True):
        cx, cy = cntrd(im, float(xi), float(yi), fwhm_col)
        if not np.isfinite(cx):
            cx, cy = float(xi), float(yi)
        flux, *_ = aper_idl(im, cx, cy, 0.5 * fwhm_col, setskyval=0.0)
        lc = float(np.log10(flux / med_star)) if np.isfinite(med_star) and flux > 0 else -9.0
        detections.append(Detection(cx, cy, float(snr[yi, xi]), flux, lc))
    detections.sort(key=lambda d: -d.snr)
    if outfile:
        np.savetxt(
            outfile,
            np.array([[d.x, d.y, d.log_contrast] for d in detections]),
            fmt="%.4f %.4f %.4f",
            header="x y log10(contrast)",
        )
    for d in detections:
        log.info(
            "detection at (%.1f, %.1f): SNR %.1f, log contrast %.2f",
            d.x,
            d.y,
            d.snr,
            d.log_contrast,
        )
    return detections
