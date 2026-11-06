"""PDI flat-fielding (``charis_pdi_flatfield.pro``).

A master flat is built from bright (minus dark) flat cubes, normalised per
channel by the median inside the field-of-view mask, and the prepared
science and sky cubes are divided by it (``<n>e_flat.fits``).

Fixes vs. IDL: the zero-pixel guard works (#215), the FOV mask comes from the
instrument table, ``suffix`` is honoured, a missing flat directory gives a
clear error and near-zero flat values are not divided by (#216).
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np

from charis_dpp.config.dataset import DatasetConfig
from charis_dpp.geometry.masks import generate_roi_mask
from charis_dpp.io.filelists import frame_filename
from charis_dpp.io.fits import read_cube, write_cube
from charis_dpp.prelims.setupdir import Dirs

log = logging.getLogger(__name__)

__all__ = ["pdi_flatfield", "make_master_flat"]


def make_master_flat(
    bright: list[Path], dark: list[Path] | None, fov: np.ndarray | None, *, floor: float = 0.05
) -> tuple[np.ndarray, object, object]:
    cubes = [read_cube(p) for p in bright]
    flat = np.nanmedian(np.stack([c.data for c in cubes]), axis=0)
    if dark:
        flat = flat - np.nanmedian(np.stack([read_cube(p).data for p in dark]), axis=0)
    for il in range(flat.shape[0]):
        sl = flat[il]
        norm = np.nanmedian(sl[fov] if fov is not None else sl)
        sl = sl / norm
        sl[~np.isfinite(sl) | (sl < floor)] = 1.0  # #215/#216: no division by ~0
        flat[il] = sl
    return flat, cubes[0].h0, cubes[0].h1


def pdi_flatfield(
    cfg: DatasetConfig | str | Path,
    *,
    brightflatdir: str | Path | None = None,
    darkflatdir: str | Path | None = None,
    master: str | Path = "master_pdiflat.fits",
    prefix: str = "n",
    suffix: str = "e",
    out_suffix: str | None = None,
) -> list[Path]:
    cfg = cfg if isinstance(cfg, DatasetConfig) else DatasetConfig.load(cfg)
    res = cfg.resolve()
    dirs = Dirs(cfg.workdir)
    pd = cfg.pdi
    master_path = Path(master) if Path(master).is_absolute() else dirs.root / master
    nums = cfg.frames("science") + cfg.frames("sky")
    if not nums:
        raise ValueError("no frames (observation.frames)")
    if brightflatdir is not None or not master_path.exists():
        bdir = Path(brightflatdir or pd.get("brightflatdir") or dirs.root / "brightflats")
        ddir = Path(darkflatdir or pd.get("darkflatdir") or dirs.root / "darkflats")
        bright = sorted(bdir.glob("*.fits"))
        if not bright:
            raise FileNotFoundError(f"no bright flat cubes in {bdir} (brightflatdir=)")
        dark = sorted(ddir.glob("*.fits")) if ddir.exists() else []
        first = read_cube(bright[0])
        fm = res.instrument.pdi.fov_mask
        scale = first.shape[1] / float(res.instrument.detector.cube_size_px)
        # charis_generate_roi, roirange=[125, 125, -26.9166] (full widths; #326, #328, #334)
        w = float(fm["width_px"]) * scale
        fov = generate_roi_mask(first.shape, w, w, float(fm["angle_deg"]))
        flat, h0, h1 = make_master_flat(bright, dark, fov)
        write_cube(master_path, flat, h0, h1)
        log.info(
            "master flat from %d bright / %d dark cubes -> %s", len(bright), len(dark), master_path
        )
    else:
        flat = read_cube(master_path).data
    out = []
    for n in nums:
        c = read_cube(dirs.prep / frame_filename(n, suffix, prefix))
        if c.data.shape != flat.shape:
            raise ValueError(f"flat shape {flat.shape} != cube shape {c.data.shape}")
        c.h1["FLATFILE"] = (master_path.name, "PDI master flat")
        p = dirs.prep / frame_filename(n, out_suffix or f"{suffix}_flat", prefix)
        write_cube(p, c.data / flat, c.h0, c.h1)
        out.append(p)
    return out
