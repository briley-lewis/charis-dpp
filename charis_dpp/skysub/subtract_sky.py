"""Sky subtraction (``charis_subtract_sky.pro`` / ``charis_pdi_subtract_sky.pro``).

A master sky cube (median of the sky frames) is subtracted from every science
cube, optionally scaled by the ratio of the outer-region medians either in the
reddest channel (``scalesky``) or per channel (``scaleindchan``, only where
both medians exceed ``skylim``).  Output: ``reduc/prep/n####<suffix>_skysub.fits``.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np

from charis_dpp.config.dataset import DatasetConfig
from charis_dpp.geometry.grids import radial_grid
from charis_dpp.io.filelists import frame_filename
from charis_dpp.io.fits import read_cube, write_cube
from charis_dpp.prelims.setupdir import Dirs

log = logging.getLogger(__name__)

__all__ = ["subtract_sky", "master_sky"]


def master_sky(paths: list[Path]) -> np.ndarray:
    cubes = [read_cube(p).data for p in paths]
    stack = np.stack(cubes)
    return np.nanmedian(stack, axis=0) if len(cubes) > 1 else cubes[0]


def subtract_sky(
    cfg: DatasetConfig | str | Path,
    *,
    scalesky: bool = False,
    scaleindchan: bool = False,
    skyrad: float | None = None,
    skylim: float | None = None,
    ladder: bool = False,
    sky_files: list[Path] | None = None,
    prefix: str = "n",
    suffix: str | None = None,
    auto_suffix: bool = True,
    save_master: bool = False,
) -> list[Path]:
    cfg = cfg if isinstance(cfg, DatasetConfig) else DatasetConfig.load(cfg)
    dirs = Dirs(cfg.workdir)
    red = cfg.reduction.get("skysub", {})
    if suffix is None:
        suffix = "e"
        if auto_suffix and any(dirs.prep.glob(f"{prefix}*e_flat.fits")):
            suffix = "e_flat"
    nums = cfg.frames("ladder" if ladder else "science")
    if not nums:
        raise ValueError("no science frames listed in observation.frames")
    if sky_files is None:
        sky_nums = cfg.frames("sky")
        if not sky_nums:
            raise ValueError("no sky frames listed (observation.frames.sky) and no sky_files given")
        sky_files = [dirs.prep / frame_filename(n, suffix, prefix) for n in sky_nums]
    missing = [p for p in sky_files if not p.exists()]
    if missing:
        raise FileNotFoundError(f"sky frames not found: {missing}")
    sky = master_sky(sky_files)
    if save_master:
        write_cube(dirs.proc / "mastersky.fits", sky)
    nlam, ny, nx = sky.shape
    skyrad = float(
        skyrad if skyrad is not None else red.get("skyrad_px") or 0.6 * (min(nx, ny) // 2)
    )
    skyrad = min(skyrad, 0.95 * (min(nx, ny) // 2))
    skylim = float(skylim if skylim is not None else red.get("skylim", 25.0))
    r = radial_grid((ny, nx))
    outputs = []
    for n in nums:
        src = dirs.prep / frame_filename(n, suffix, prefix)
        c = read_cube(src)
        im = c.data
        if not (scalesky or scaleindchan):
            im = im - sky
        elif scalesky:
            targ, sk = im[-1], sky[-1]
            good = (r > skyrad) & (targ != 0) & np.isfinite(targ) & np.isfinite(sk)
            scale = np.median(targ[good]) / np.median(sk[good]) if good.any() else 1.0
            im = im - sky * scale
            c.h1["SKYSCALE"] = (float(scale), "sky scale factor (last channel)")
        else:
            im = im.copy()
            for il in range(nlam):
                targ, sk = im[il], sky[il]
                good = (r > skyrad) & (targ != 0) & np.isfinite(targ) & np.isfinite(sk)
                scale = 1.0
                if good.any():
                    mt, ms = np.median(targ[good]), np.median(sk[good])
                    if mt > skylim and ms > skylim:
                        scale = mt / ms
                im[il] = targ - sk * scale
                c.h1[f"SKYSC_{il}"] = (float(scale), "sky scale factor")
        c.h1["SKYSUB"] = (True, "master sky subtracted")
        c.h1["NSKY"] = (len(sky_files), "number of sky frames in master")
        out = dirs.prep / frame_filename(n, f"{suffix}_skysub", prefix)
        write_cube(out, im, c.h0, c.h1)
        outputs.append(out)
    log.info("sky-subtracted %d cubes using %d sky frames", len(outputs), len(sky_files))
    return outputs
