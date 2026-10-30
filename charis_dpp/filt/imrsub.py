"""Spatial filtering before PSF subtraction (``filt/charis_imrsub.pro``).

Per slice: radial-profile subtraction (``prad``), azimuthal-profile
subtraction (``pang``), moving-box median high-pass (``medbox``), spider mask
(``mask``).  Optionally subtracts the sequence-median cube (``psfsub``) and/or
a simple SDI (``sdi``) before derotating; writes ``reduc/rsub/n####rsub[_fc].fits``
and the classical-ADI products ``reduc/proc/imcubefinal_{avg,derot}[_collapsed].fits``.

Fixes vs. IDL: the ``reg_cal`` suffix is chosen only if the files exist (#58),
PA/rsat are matched by frame number (#59), empty masks no longer corrupt the
last pixel (#60), ``nonorthup`` works (#61), outputs go to ``reduc/proc``.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np

from charis_dpp.config.dataset import DatasetConfig
from charis_dpp.geometry.grids import image_center
from charis_dpp.geometry.masks import spider_mask
from charis_dpp.image.filters import azimuthal_profile, idl_filter_image, radial_profile_2d
from charis_dpp.image.interp import rescale_radial, rotate_cube
from charis_dpp.image.wcs import northup
from charis_dpp.io.filelists import frame_filename
from charis_dpp.io.fits import read_cube, write_cube
from charis_dpp.io.reduclog import ReducLog
from charis_dpp.prelims.setupdir import Dirs

log = logging.getLogger(__name__)

__all__ = ["imrsub", "align_speckles"]


def align_speckles(cube: np.ndarray, wav: np.ndarray, ref: int, center=None) -> np.ndarray:
    """Rescale every channel so its speckle pattern matches channel ``ref`` (``charis_alignspeckle``)."""
    out = np.empty_like(cube)
    for il in range(cube.shape[0]):
        out[il] = cube[il] if il == ref else rescale_radial(cube[il], wav[ref] / wav[il], center)
    return out


def imrsub(
    cfg: DatasetConfig | str | Path,
    *,
    prad: bool = False,
    pang: bool = False,
    psfsub: bool = False,
    sdi: bool = False,
    medbox: int | None = None,
    mask: bool = False,
    rmax: float | None = None,
    fc: bool = False,
    northup_final: bool = True,
    prefix: str = "n",
    suffix: str | None = None,
) -> list[Path]:
    cfg = cfg if isinstance(cfg, DatasetConfig) else DatasetConfig.load(cfg)
    res = cfg.resolve()
    dirs = Dirs(cfg.workdir)
    dirs.rsub.mkdir(parents=True, exist_ok=True)
    nums = cfg.frames("science")
    if suffix is None:
        suffix = (
            "reg_cal" if (dirs.reg / frame_filename(nums[0], "reg_cal", prefix)).exists() else "reg"
        )
    in_suffix = suffix + ("_fc" if fc else "")
    out_suffix = "rsub" + ("_fc" if fc else "")
    rlog = ReducLog.load(dirs.reduclog)
    rlog.require(nums)
    pa = rlog.pa(nums)
    rsat = rlog.rsat(nums)
    first = read_cube(dirs.reg / frame_filename(nums[0], in_suffix, prefix))
    ny, nx = first.shape
    center = res.star_center or image_center((ny, nx))
    rmax = float(rmax if rmax is not None else min(nx, ny) // 2)
    spm = None
    if mask:
        sp = cfg.reduction.get("spider", {})
        spm = spider_mask(
            (ny, nx),
            float(sp.get("width_px", 3.0)),
            sp.get("angles_deg", [0, 90, 180, 270]),
            center,
            fractional=False,
        )
    outs = []
    stack = []
    headers = []
    for i, n in enumerate(nums):
        c = read_cube(dirs.reg / frame_filename(n, in_suffix, prefix))
        cube = c.data.copy()
        for il in range(c.nlam):
            im = cube[il]
            bad = (im == 0) | ~np.isfinite(im)
            im2 = np.where(bad, np.nan, im)
            if prad:
                pr = radial_profile_2d(im2, 2.0, max(rsat[i], 0.0), rmax, center)
                im = np.where(bad, im, im - pr)
            if pang:
                _, _, pmap = azimuthal_profile(
                    np.where(bad, np.nan, im), max(rsat[i], 1.0), rmax, center=center
                )
                im = im - pmap
            if medbox:
                imm = np.where((im != 0) & np.isfinite(im), im, np.nan)
                im = imm - idl_filter_image(imm, median=int(medbox))
            if spm is not None:
                im = np.where(spm, im, np.nan)
            cube[il] = im
        c.h1["RSUB"] = (1, "spatially filtered (imrsub)")
        c.h1["RSUBPRAD"] = int(prad)
        c.h1["RSUBMBOX"] = int(medbox or 0)
        out = dirs.rsub / frame_filename(n, out_suffix, prefix)
        write_cube(out, cube, c.h0, c.h1)
        outs.append(out)
        stack.append(cube.astype(np.float32))
        headers.append((c.h0, c.h1))
    stack = np.stack(stack)  # (nf, nλ, ny, nx)
    avg = np.nanmedian(stack, axis=0)
    dirs.proc.mkdir(parents=True, exist_ok=True)
    h0, h1 = headers[0]
    write_cube(dirs.proc / f"imcubefinal_avg{'_fc' if fc else ''}.fits", avg, h0, h1)
    write_cube(
        dirs.proc / f"imcubefinal_avg_collapsed{'_fc' if fc else ''}.fits",
        np.nanmedian(avg, axis=0),
        h0,
        h1,
    )

    ref_cube = avg if psfsub else None
    derot = []
    for i in range(len(nums)):
        use = stack[i].astype(float)
        if ref_cube is not None:
            use = use - ref_cube
        if sdi:
            use = _simple_sdi(use, res.wavelengths_um, center)
        derot.append(rotate_cube(use, pa[0] - pa[i], center))
    final = np.nanmedian(np.stack(derot), axis=0)
    h1f = headers[0][1].copy()
    if northup_final:
        final, _ = northup(final, h1f, res.instrument.north_offset_deg, center=center)
    tag = "_fc" if fc else ""
    write_cube(dirs.proc / f"imcubefinal_derot{tag}.fits", final, h0, h1f)
    write_cube(
        dirs.proc / f"imcubefinal_derot_collapsed{tag}.fits", np.nanmedian(final, axis=0), h0, h1f
    )
    log.info("imrsub: wrote %d filtered cubes to %s", len(outs), dirs.rsub)
    return outs


def _simple_sdi(cube: np.ndarray, wav: np.ndarray, center) -> np.ndarray:
    out = cube.copy()
    for il in range(cube.shape[0]):
        aligned = align_speckles(cube, wav, il, center)
        ref = np.nanmedian(aligned, axis=0)
        w = np.nansum(aligned[il] * ref) / np.nansum(ref * ref)
        out[il] = cube[il] - w * ref
    return out
