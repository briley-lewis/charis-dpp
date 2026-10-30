"""Sky subtraction for PDI data (``charis_pdi_subtract_sky.pro``): the generic
sky subtraction on the flat-fielded cubes when they exist."""

from __future__ import annotations

from pathlib import Path

from charis_dpp.config.dataset import DatasetConfig
from charis_dpp.io.filelists import frame_filename
from charis_dpp.prelims.setupdir import Dirs
from charis_dpp.skysub.subtract_sky import subtract_sky

__all__ = ["pdi_subtract_sky"]


def pdi_subtract_sky(
    cfg: DatasetConfig | str | Path, *, prefix: str = "n", suffix: str | None = None, **kw
):
    cfg = cfg if isinstance(cfg, DatasetConfig) else DatasetConfig.load(cfg)
    dirs = Dirs(cfg.workdir)
    if suffix is None:
        n0 = cfg.frames("science")[0]
        suffix = "e_flat" if (dirs.prep / frame_filename(n0, "e_flat", prefix)).exists() else "e"
    return subtract_sky(cfg, prefix=prefix, suffix=suffix, **kw)
