"""Registration of the two Wollaston beams (``charis_pdi_register_cube.pro``, which
absorbed ``charis_pdi_split_pols.pro``; ``method='xcorr'`` absorbs
``python/nonstandard_pdi_reg.py``).

Each beam is registered to the frame centre with the generic registration
(satellite spots around the beam centre, guessed at ``±polcent_offset``),
giving ``<n>leftreg.fits`` / ``<n>rightreg.fits`` plus ``reduc_left.log`` and
``reduc_right.log`` (``reduc.log`` keeps the right-beam solution for the
total-intensity PSF subtraction).  Pixels outside the raw field of view and
pixels finite in only one beam are NaN; ``<n>reg.fits`` is the beam sum.

``method='xcorr'`` registers each beam by cross-correlation against the first
frame after centring that frame on its beam's centroid inside a box
(``[pdi.xcorr]`` config: ``left_box``/``right_box`` as ``[x0, x1, y0, y1]``
and optional ``left_polygon``/``right_polygon`` FOV masks).

Fixes vs. IDL: the file tests use counts (#217), no shell ``mv``/``cp``
(#219), the parity handling is explicit (#310: with ``parity = -1`` the
x-offset is mirrored and the beams swap, not the y-offset).
"""

from __future__ import annotations

import logging
import shutil
from pathlib import Path

import numpy as np

from charis_dpp.config.dataset import DatasetConfig
from charis_dpp.geometry.grids import image_center
from charis_dpp.image.centroid import fft_xcorr_shift
from charis_dpp.image.interp import shift_sub
from charis_dpp.io.filelists import frame_filename
from charis_dpp.io.fits import read_cube, write_cube
from charis_dpp.io.reduclog import ReducLog
from charis_dpp.prelims.setupdir import Dirs
from charis_dpp.reg.register_cube import register_cube

log = logging.getLogger(__name__)

__all__ = ["pdi_register_cube", "polcent_offset_px"]


def polcent_offset_px(cfg: DatasetConfig, res) -> np.ndarray:
    """Beam offset ``(x, y)`` of the right beam from the frame centre, scaled to the frame size."""
    v = cfg.pdi.get("polcent_offset_px", "instrument")
    inst = res.instrument
    off = np.array(inst.pdi.polcent_offset_px if v == "instrument" else v, dtype=float)
    shape = res.frame_shape
    if shape is not None:
        off = off * (float(shape[1]) / float(inst.detector.cube_size_px))
    if inst.parity == -1 and v == "instrument":
        off = off * np.array([-1.0, 1.0])
    return off


def _polygon_mask(shape, poly) -> np.ndarray:
    from matplotlib.path import Path as MplPath

    ny, nx = shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    return (
        MplPath(np.asarray(poly, dtype=float))
        .contains_points(np.column_stack([xx.ravel(), yy.ravel()]))
        .reshape(ny, nx)
    )


def _xcorr_register(cfg, res, dirs, nums, suffix, prefix, offset, side, reduclog_path):
    xc = cfg.pdi.get("xcorr", {})
    box = xc.get(f"{side}_box")
    poly = xc.get(f"{side}_polygon")
    first = read_cube(dirs.prep / frame_filename(nums[0], suffix, prefix))
    shape = first.shape
    cx, cy = image_center(shape)
    sign = -1.0 if side == "left" else 1.0
    gx, gy = cx + sign * offset[0], cy + sign * offset[1]
    # reference: first frame centred on the beam centroid inside the box (or ±offset)
    ref = np.empty_like(first.data)
    for il in range(first.nlam):
        im = np.where(np.isfinite(first.data[il]), first.data[il], 0.0)
        if box is not None:
            x0, x1, y0, y1 = (int(v) for v in box)
            m = np.zeros(shape, bool)
            m[y0:y1, x0:x1] = True
            w = np.where(m, im, 0.0)
            w = np.clip(w, 0, None)
            yy, xx = np.mgrid[0 : shape[0], 0 : shape[1]]
            bx, by = float((w * xx).sum() / w.sum()), float((w * yy).sum() / w.sum())
        else:
            bx, by = gx, gy
        ref[il] = shift_sub(im, cx - bx, cy - by)
    roi = np.ones(shape, bool)
    if poly is not None:
        roi = _polygon_mask(shape, poly)
    rl = ReducLog(path=reduclog_path)
    for n in nums:
        c = read_cube(dirs.prep / frame_filename(n, suffix, prefix))
        out = np.empty_like(c.data)
        cens = []
        for il in range(c.nlam):
            im = np.where(np.isfinite(c.data[il]), c.data[il], 0.0)
            dx, dy = fft_xcorr_shift(ref[il], im, roi)
            px, py = cx - dx, cy - dy
            out[il] = shift_sub(c.data[il], dx, dy)
            c.h1[f"SHIFT_{il}"] = (f"{dx:.4f} {dy:.4f}", "shift applied (x y)")
            c.h1[f"PSFC_{il}"] = (f"{px:.4f} {py:.4f}", "beam centre before shift (x y)")
            cens.append((px, py))
        mx, my = np.nanmedian(np.array(cens), axis=0)
        c.h1["PSFCENTX"], c.h1["PSFCENTY"] = float(mx), float(my)
        c.h1["REGMETH"] = ("xcorr", "registration method")
        write_cube(dirs.reg / frame_filename(n, f"{side}reg", prefix), out, c.h0, c.h1)
        rl.add(
            n,
            float(mx),
            float(my),
            -1.0,
            float(c.get("HA", 0.0)) if not isinstance(c.get("HA"), str) else 0.0,
            float(c.get("PA", c.get("PARANG", 0.0))),
        )
    rl.save()
    return rl


def pdi_register_cube(
    cfg: DatasetConfig | str | Path,
    *,
    method: str = "sats",
    polcent_offset: tuple[float, float] | None = None,
    astrogrid: float | None = None,
    nomask: bool = False,
    prefix: str = "n",
    suffix: str | None = None,
    ladder: bool = False,
    **register_kwargs,
) -> tuple[ReducLog, ReducLog]:
    cfg = cfg if isinstance(cfg, DatasetConfig) else DatasetConfig.load(cfg)
    res = cfg.resolve()
    dirs = Dirs(cfg.workdir)
    dirs.reg.mkdir(parents=True, exist_ok=True)
    nums = cfg.frames("ladder" if ladder else "science")
    if suffix is None:
        for cand in ("e_flat_skysub", "e_flat", "e_skysub", "e"):
            if (dirs.prep / frame_filename(nums[0], cand, prefix)).exists():
                suffix = cand
                break
        else:
            raise FileNotFoundError("no prepared cubes for the science frames")
    off = (
        np.asarray(polcent_offset, dtype=float)
        if polcent_offset is not None
        else polcent_offset_px(cfg, res)
    )
    astro = astrogrid if astrogrid is not None else cfg.observation.get("astrogrid_lambda_d")
    if astro is None:
        astro = float(res.instrument.pdi.astrogrid_lambda_d)
    left_log = dirs.root / ("reduc_left.log" if not ladder else "reduc_left_ladder.log")
    right_log = dirs.root / ("reduc_right.log" if not ladder else "reduc_right_ladder.log")
    if method == "xcorr":
        rl_l = _xcorr_register(cfg, res, dirs, nums, suffix, prefix, off, "left", left_log)
        rl_r = _xcorr_register(cfg, res, dirs, nums, suffix, prefix, off, "right", right_log)
    else:
        rl_l = register_cube(
            cfg,
            method=method,
            suffix=suffix,
            prefix=prefix,
            ladder=ladder,
            guessoffsets=(-off[0], -off[1]),
            astrogrid=astro,
            out_suffix="leftreg",
            reduclog_path=left_log,
            **register_kwargs,
        )
        rl_r = register_cube(
            cfg,
            method=method,
            suffix=suffix,
            prefix=prefix,
            ladder=ladder,
            guessoffsets=(off[0], off[1]),
            astrogrid=astro,
            out_suffix="rightreg",
            reduclog_path=right_log,
            **register_kwargs,
        )
    if not ladder:
        shutil.copyfile(right_log, dirs.reduclog)  # total-intensity PSF subtraction uses reduc.log
    # combine: mask outside the raw FOV, keep pixels finite in both beams, write the sum
    skysubbed = "skysub" in suffix
    for n in nums:
        c0 = read_cube(dirs.prep / frame_filename(n, suffix, prefix))
        im0 = c0.data.copy()
        if skysubbed:
            raw = dirs.prep / frame_filename(n, "e", prefix)
            fov0 = read_cube(raw).data if raw.exists() else im0
            outside = ~np.isfinite(fov0) | (fov0 == 0)
        else:
            outside = ~np.isfinite(im0) | (im0 == 0)
        if not nomask:
            im0[outside] = np.nan
        left = read_cube(dirs.reg / frame_filename(n, "leftreg", prefix))
        right = read_cube(dirs.reg / frame_filename(n, "rightreg", prefix))
        lreg = np.empty_like(im0)
        rreg = np.empty_like(im0)
        for il in range(c0.nlam):
            ldx, ldy = (float(v) for v in str(left.h1[f"SHIFT_{il}"]).split()[:2])
            rdx, rdy = (float(v) for v in str(right.h1[f"SHIFT_{il}"]).split()[:2])
            lreg[il] = shift_sub(im0[il], ldx, ldy)
            rreg[il] = shift_sub(im0[il], rdx, rdy)
        if not nomask:
            bad = ~(np.isfinite(lreg) & np.isfinite(rreg))
            lreg[bad] = np.nan
            rreg[bad] = np.nan
        tot = lreg + rreg
        write_cube(dirs.reg / frame_filename(n, "leftreg", prefix), lreg, left.h0, left.h1)
        write_cube(dirs.reg / frame_filename(n, "rightreg", prefix), rreg, right.h0, right.h1)
        write_cube(dirs.reg / frame_filename(n, "reg", prefix), tot, left.h0, left.h1)
    log.info("pdi_register_cube: %d frames, beams at ±(%.2f, %.2f) px", len(nums), off[0], off[1])
    return rl_l, rl_r
