"""Empirical PSF from the satellite spots (``charis_makeemppsf.pro``).

For every channel and frame the four spots are cut out (after radial-profile
subtraction), re-centred to the stamp centre with sub-pixel shifts, background
subtracted, clipped of negative outliers and normalised; the median over
spots and frames is the PSF.  Output: ``psfmodel/psfcube_med.fits`` (with the
wavelengths in its header, unlike IDL — #98) and ``psf_medcol.fits``.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
from astropy.io.fits import Header

from charis_dpp.config.dataset import DatasetConfig
from charis_dpp.geometry.grids import radial_grid
from charis_dpp.image.filters import fixpix, idl_filter_image, radial_profile_2d
from charis_dpp.image.interp import shift_sub, subarr
from charis_dpp.io.filelists import frame_filename
from charis_dpp.io.fits import get_psf_center, get_sat_positions, read_cube, write_cube, write_image
from charis_dpp.prelims.setupdir import Dirs

log = logging.getLogger(__name__)

__all__ = ["make_empirical_psf", "psf_header"]


def psf_header(wav_um: np.ndarray, extra: dict | None = None) -> Header:
    h = Header()
    for l, w in enumerate(wav_um):
        h[f"WAVEL{l}"] = (float(w), "channel wavelength [um]")
    h["NLAM"] = int(len(wav_um))
    for k, v in (extra or {}).items():
        h[k] = v
    return h


def make_empirical_psf(
    cfg: DatasetConfig | str | Path,
    *,
    psfsize: int | None = None,
    fsize: int | None = None,
    unsat: bool = False,
    filt: bool = False,
    nozero: bool = False,
    nonorm: bool = False,
    ladder: bool = False,
    files: list[Path] | None = None,
    prefix: str = "n",
    suffix: str | None = None,
    outname: str | None = None,
) -> Path:
    cfg = cfg if isinstance(cfg, DatasetConfig) else DatasetConfig.load(cfg)
    res = cfg.resolve()
    dirs = Dirs(cfg.workdir)
    dirs.psfmodel.mkdir(parents=True, exist_ok=True)
    pc = cfg.reduction.get("psf", {})
    dpsf = int(psfsize or pc.get("dpsf_px", 21))
    fsize = int(fsize or pc.get("fsize_px", dpsf))
    if files is None:
        nums = cfg.frames("ladder" if ladder else "science")
        if suffix is None:
            suffix = (
                "reg_cal"
                if (dirs.reg / frame_filename(nums[0], "reg_cal", prefix)).exists()
                else "reg"
            )
        files = [dirs.reg / frame_filename(n, suffix, prefix) for n in nums]
    nlam = res.nlam
    fwhm = res.fwhm_px
    rstamp = radial_grid((dpsf, dpsf))
    cubes = [read_cube(p) for p in files]
    psf_per_frame = np.full((len(cubes), nlam, dpsf, dpsf), np.nan)
    for nf, c in enumerate(cubes):
        sats = get_sat_positions(c.h1, nlam) if not unsat else None
        xc, yc = get_psf_center(c.h1) or (c.shape[1] // 2, c.shape[0] // 2)
        for il in range(nlam):
            im = c.data[il].copy()
            if not unsat:
                bad = (im <= 0) | ~np.isfinite(im)
                im2 = np.where(bad, np.nan, im)
                pr = radial_profile_2d(im2, 1.0, 0.0, None)
                im[~bad] -= pr[~bad]
                if filt:
                    im = im - idl_filter_image(im, median=fsize)
                stamps = []
                for q in range(4):
                    pos = sats[il, q]
                    if not np.isfinite(pos).all() or np.any(pos > 1000):
                        continue
                    base = np.floor(pos)
                    s = subarr(im, dpsf, (int(base[0]), int(base[1])))
                    s = shift_sub(s, -(pos[0] - base[0]), -(pos[1] - base[1]))
                    if not nozero:
                        s = s - np.nanmedian(s[rstamp > 3 * fwhm[il]])
                    s[s < -0.01 * np.nanmax(s)] = 0.0
                    if not nonorm:
                        s = s / np.nansum(s)
                    stamps.append(s)
                if stamps:
                    psf_per_frame[nf, il] = np.nanmedian(np.stack(stamps), axis=0)
            else:
                s = subarr(im, dpsf, (int(round(xc)), int(round(yc))))
                if not nozero:
                    s = s - np.nanmedian(s[rstamp > 2 * fwhm[il]])
                psf_per_frame[nf, il] = s
    if len(cubes) > 1:
        psf = np.nanmedian(psf_per_frame, axis=0)
    else:
        psf = np.stack([fixpix(np.nan_to_num(sl), 3) for sl in psf_per_frame[0]])
    if not nonorm:
        psf = psf / np.nansum(psf, axis=(1, 2), keepdims=True)
    tag = "" if not nonorm else "_nonorm"
    name = outname or f"psfcube_med{tag}.fits"
    out = dirs.psfmodel / name
    hdr = psf_header(
        res.wavelengths_um, {"PSFTYPE": "empirical", "NFRAMES": len(cubes), "NORMED": not nonorm}
    )
    write_cube(out, psf, Header(), hdr)
    col = np.nanmedian(psf, axis=0)
    write_image(
        dirs.psfmodel
        / (
            Path(name).stem.replace("psfcube", "psf") + "col.fits"
            if outname
            else f"psf_medcol{tag}.fits"
        ),
        col,
    )
    log.info("empirical PSF from %d frames -> %s", len(cubes), out)
    return out
