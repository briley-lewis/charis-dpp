"""Synthetic PSF cube (``charis_makesynthpsf.pro``): flat disk or Gaussian per channel."""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
from astropy.io.fits import Header

from charis_dpp.config.dataset import DatasetConfig
from charis_dpp.geometry.grids import radial_grid
from charis_dpp.image.filters import idl_smooth
from charis_dpp.io.fits import write_cube, write_image
from charis_dpp.prelims.setupdir import Dirs
from charis_dpp.specphotcal.emp_psf import psf_header

log = logging.getLogger(__name__)

__all__ = ["make_synthetic_psf", "gaussian_psf_cube"]


def gaussian_psf_cube(fwhm_px: np.ndarray, size: int, rmax: float | None = None) -> np.ndarray:
    r = radial_grid((size, size))
    out = np.empty((len(fwhm_px), size, size))
    for il, f in enumerate(fwhm_px):
        sig = f / 2.35482
        g = np.exp(-0.5 * (r / sig) ** 2)
        if rmax is not None:
            g[r > rmax] = 0.0
        out[il] = g / g.sum()
    return out


def make_synthetic_psf(
    cfg: DatasetConfig | str | Path,
    *,
    kind: str = "gaussian",
    fwhm: float | np.ndarray | None = None,
    rmax: float | None = None,
    psfsize: int | None = None,
    smoothval: int = 2,
    outname: str | None = None,
) -> Path:
    cfg = cfg if isinstance(cfg, DatasetConfig) else DatasetConfig.load(cfg)
    res = cfg.resolve()
    dirs = Dirs(cfg.workdir)
    dirs.psfmodel.mkdir(parents=True, exist_ok=True)
    size = int(psfsize or cfg.reduction.get("psf", {}).get("dpsf_px", 21))
    if fwhm is None:
        fw = res.fwhm_px
    else:
        fw = np.broadcast_to(np.asarray(fwhm, dtype=float), (res.nlam,)).copy()
        if fw[0] < 1:  # arcsec given
            fw = fw / res.instrument.pixscale_arcsec
    if kind == "flat":
        r = radial_grid((size, size))
        col = (r <= (rmax if rmax is not None else 5.0)).astype(float)
        col = idl_smooth(col, smoothval)
        col /= col.sum()
        psf = np.repeat(col[None], res.nlam, axis=0)
    elif kind == "gaussian":
        psf = gaussian_psf_cube(fw, size, rmax)
        col = np.nanmedian(psf, axis=0)
    else:
        raise ValueError("kind must be 'flat' or 'gaussian'")
    name = outname or "synthcube_med.fits"
    out = dirs.psfmodel / name
    write_cube(out, psf, Header(), psf_header(res.wavelengths_um, {"PSFTYPE": kind}))
    write_image(dirs.psfmodel / (Path(name).stem + "col.fits"), col)
    log.info("synthetic %s PSF -> %s", kind, out)
    return out
