"""Spectrophotometric calibration of the two beams (``charis_pdi_specphot_cal.pro``):
the beam sum is calibrated with the generic routine and the resulting
per-channel scale factors and star-flux cards are applied to each beam
(``<n>leftreg_cal.fits`` / ``<n>rightreg_cal.fits``).

Fixes vs. IDL: FOV-mask guards (#220); the temporary files are removed with
``Path.unlink`` and the prefix is honoured; ``modamp`` and ``notimeratio`` are
not silently dropped with ``unsat`` (#221).
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np

from charis_dpp.config.dataset import DatasetConfig
from charis_dpp.geometry.masks import fov_mask
from charis_dpp.io.filelists import frame_filename
from charis_dpp.io.fits import read_cube, write_cube
from charis_dpp.prelims.setupdir import Dirs
from charis_dpp.specphotcal.specphot_cal import specphot_cal

log = logging.getLogger(__name__)

__all__ = ["pdi_specphot_cal"]

_COPY_KEYS = ("FLUXUNIT", "SKY_IN", "SKY_OUT")


def pdi_specphot_cal(
    cfg: DatasetConfig | str | Path,
    *,
    fov_mask_pars: tuple[float, float] | None = None,
    prefix: str = "n",
    left_suffix: str = "leftreg",
    right_suffix: str = "rightreg",
    **specphot_kwargs,
) -> list[Path]:
    cfg = cfg if isinstance(cfg, DatasetConfig) else DatasetConfig.load(cfg)
    dirs = Dirs(cfg.workdir)
    nums = cfg.frames("science")
    tmp = [dirs.reg / frame_filename(n, "_tmp", prefix) for n in nums]
    first = read_cube(dirs.reg / frame_filename(nums[0], left_suffix, prefix))
    outside = None
    if fov_mask_pars is not None:
        outside = ~fov_mask(first.shape, float(fov_mask_pars[0]), float(fov_mask_pars[1]))
    for n, t in zip(nums, tmp, strict=True):
        left = read_cube(dirs.reg / frame_filename(n, left_suffix, prefix))
        right = read_cube(dirs.reg / frame_filename(n, right_suffix, prefix))
        tot = left.data + right.data
        if outside is not None:
            tot[:, outside] = np.nan
        tot = np.where(np.isfinite(tot), tot, 0.0)
        write_cube(t, tot, left.h0, left.h1)
    try:
        specphot_cal(cfg, prefix=prefix, suffix="_tmp", **specphot_kwargs)
        out = []
        for n in nums:
            cal = read_cube(dirs.reg / frame_filename(n, "_tmp_cal", prefix))
            fscale = np.array([float(cal.h1[f"FSCALE{l}"]) for l in range(cal.nlam)])
            for side in (left_suffix, right_suffix):
                c = read_cube(dirs.reg / frame_filename(n, side, prefix))
                for k in _COPY_KEYS:
                    if k in cal.h1:
                        c.h1[k] = (cal.h1[k], cal.h1.comments[k])
                for l in range(cal.nlam):
                    for k in (f"FSTAR_{l}", f"R_AP{l}", f"FSCALE{l}", f"CERR{l}"):
                        if k in cal.h1:
                            c.h1[k] = (cal.h1[k], cal.h1.comments[k])
                p = dirs.reg / frame_filename(n, f"{side}_cal", prefix)
                write_cube(p, c.data * fscale[:, None, None], c.h0, c.h1)
                out.append(p)
    finally:
        for n, t in zip(nums, tmp, strict=True):
            t.unlink(missing_ok=True)
            (dirs.reg / frame_filename(n, "_tmp_cal", prefix)).unlink(missing_ok=True)
    log.info("pdi_specphot_cal: calibrated %d beam pairs", len(nums))
    return out
