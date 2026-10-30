"""Forward-model disk images through the PSF subtraction
(``charis_{aloci,klip}_disk_fwdmod.pro``).

For each model: scale it to the data inside the region of interest, convolve
with the intrinsic PSF per channel, rotate it into every frame, run the
recorded reduction in forward-model mode, re-scale the processed model to
the data by weighted least squares and compute the binned reduced χ².

Outputs under ``reduc/proc/``: ``model_in/<model>_inputcolsc.fits`` (scaled
input), ``model_psfsub/<model>_psfsubcol.fits`` (processed, scaled; header
``SFAC``, ``CSQ_COL``, ``NPTS`` and the model parameters), and with
``calcattenfact`` the per-channel cubes plus ``model_atten/<model>_relattencube.fits``;
``model_stat/fit_outcomes_<timestamp>.dat`` summarises the fits.

Fixes vs. IDL: all reduction types work because their keywords exist (#177);
the cube and header are read from the same directory (#178); ``snrlim`` is
honoured (#179); ``rsub`` is read from the right key (#180); ``maskrad``
defaults to 0 (#181); the Gaussian PSF uses the frame size (#182); header keys
are valid FITS keys on the right header (#183); the outcomes table widths are
sane (#184); zero-noise pixels are excluded (#185); the model is rotated once
per frame, or really smeared with ``nsc`` sub-steps (#186); one rotation
convention shared with the planet code (#187).
"""

from __future__ import annotations

import logging
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

from charis_dpp.config.dataset import DatasetConfig
from charis_dpp.filt.imrsub import imrsub
from charis_dpp.geometry.apertures import sumaper_im
from charis_dpp.geometry.grids import image_center, radial_grid
from charis_dpp.image.snr import snr_map
from charis_dpp.image.wcs import rotation_from_header
from charis_dpp.io.filelists import frame_filename
from charis_dpp.io.fits import read_cube, write_cube
from charis_dpp.partialsub.common import ReducedCube, build_psf, run_reduction
from charis_dpp.partialsub.diskfit import (
    MODEL_KEYS,
    DiskFitResult,
    DiskModel,
    chisq_map,
    convolve_model,
    fit_outcomes_header,
    fit_outcomes_line,
    load_disk_model,
    model_in_frame,
    roi_mask,
    scale_model_ls,
)
from charis_dpp.prelims.setupdir import Dirs
from charis_dpp.tools.chisq import binned_chisq

log = logging.getLogger(__name__)

__all__ = ["disk_fwdmod", "reference_noise"]


def reference_noise(
    refcol: np.ndarray, fwhm: float, rmax: float, center, *, maskrad: float = 0.0
) -> tuple[np.ndarray, np.ndarray]:
    """(aperture-summed collapsed image, noise map) with the inner ``maskrad`` filled by the maximum valid noise."""
    r = snr_map(refcol, fwhm, rmax, center=center, zero=True, finite=True)
    sig = r.noisemap.copy()
    if maskrad > 0:
        rad = radial_grid(refcol.shape, center)
        ok = (rad > maskrad) & np.isfinite(sig) & (sig > 0)
        if ok.any():
            sig[rad <= maskrad] = np.max(sig[ok])
    return r.convimage, sig


def _model_cards(h, params: dict[str, float], scale: float, chisq: float, nbins: int) -> None:
    h["SFAC"] = (scale, "model scale factor (weighted LS in the ROI)")
    h["CSQ_COL"] = (chisq, "reduced chi^2 of the binned collapsed image")
    h["NPTS"] = (nbins, "number of resolution-element bins in the ROI")
    for k in MODEL_KEYS:
        if k in params:
            h[k] = (params[k], "disk model parameter")


def disk_fwdmod(
    cfg: DatasetConfig | str | Path,
    models: list[str | Path],
    *,
    reducname: str = "final.fits",
    psf: str | Path = "empirical",
    noconv: bool = False,
    absflux: bool = False,
    synthscale: float | None = None,
    simplescale: bool = False,
    manscale: bool = False,
    roirange: tuple[float, float, float] | None = None,
    roioffset: tuple[float, float] = (0.0, 0.0),
    snrlim: float | None = None,
    lrmin: float | None = None,
    lrmax: float | None = None,
    maskrad: float = 0.0,
    calcattenfact: bool = False,
    reductype: int = 0,
    sdi_reducname: str | None = None,
    refpath: str | Path | None = None,
    nsc: int | None = None,
    prefix: str = "n",
) -> list[DiskFitResult]:
    cfg = cfg if isinstance(cfg, DatasetConfig) else DatasetConfig.load(cfg)
    res = cfg.resolve()
    dirs = Dirs(cfg.workdir)
    inst = res.instrument
    fm = cfg.reduction.get("fwdmod", {})
    nsc = int(nsc if nsc is not None else fm.get("nsc", 21))
    if not models:
        raise ValueError("give at least one disk model FITS file")
    rc = ReducedCube.load(dirs.proc / reducname, res.nlam)
    sdi = ReducedCube.load(dirs.proc / sdi_reducname, res.nlam) if sdi_reducname else None
    stem = Path(reducname).stem
    refcube = rc.cube.data
    colpath = dirs.proc / f"{stem}_collapsed.fits"
    refc = read_cube(colpath)
    refcol = refc.data[0] if refc.data.ndim == 3 else refc.data
    if reductype in (2, 3) and sdi is not None:
        sc = read_cube(dirs.proc / f"{Path(sdi_reducname).stem}_collapsed.fits")
        refcol = sc.data[0] if sc.data.ndim == 3 else sc.data
    shape = refcol.shape
    center = res.star_center or image_center(shape)
    fwhm = res.fwhm_px
    medfwhm = float(np.median(fwhm))
    lrmin = float(lrmin if lrmin is not None else fm.get("lrmin_px", 5.0))
    lrmax = float(lrmax if lrmax is not None else min(shape) / 2.0)
    zrmax = float(rc.params.get("rmax", lrmax))
    refcol_con, sigdata = reference_noise(refcol, medfwhm, zrmax, center, maskrad=maskrad)
    with np.errstate(divide="ignore", invalid="ignore"):
        snrmap = refcol_con / sigdata
    roi = roi_mask(
        shape,
        center=center,
        lrmin=lrmin,
        lrmax=lrmax,
        roirange=roirange,
        roioffset=roioffset,
        snrmap=snrmap if snrlim is not None else None,
        snrlim=snrlim,
    )
    cal = read_cube(dirs.reg / frame_filename(cfg.frames("science")[0], "reg_cal", prefix))
    r_ap = (
        rc.r_ap
        if np.isfinite(rc.r_ap).all()
        else np.array([float(cal.h1[f"R_AP{l}"]) for l in range(res.nlam)])
    )
    psf_cube = None if noconv else build_psf(cfg, res, r_ap, kind=psf)
    for d in ("model_in", "model_psfsub", "model_stat", "model_atten"):
        (dirs.proc / d).mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%S")
    outcomes = dirs.proc / "model_stat" / f"fit_outcomes_{stamp}.dat"
    outcomes.write_text(fit_outcomes_header() + "\n")
    nums = cfg.frames("science")
    dec = float(cfg.target.get("dec_deg", 0.0))
    results: list[DiskFitResult] = []
    for mpath in models:
        model: DiskModel = load_disk_model(mpath)
        if model.data.shape[-2:] != shape:
            raise ValueError(f"{mpath}: model shape {model.data.shape} != data {shape}")
        if model.is_cube and model.data.shape[0] != res.nlam:
            raise ValueError(f"{mpath}: model cube has {model.data.shape[0]} channels")
        # initial scaling to the (unprocessed) collapsed data
        mcol = np.nanmedian(model.data, axis=0) if model.is_cube else model.data
        guess = float(synthscale) if synthscale is not None else scale_model_ls(refcol, mcol, roi)
        mdata = model.data if absflux else model.data * guess
        # per-channel model, convolved with the intrinsic PSF
        chan = np.empty((res.nlam,) + shape)
        for il in range(res.nlam):
            m = mdata[il] if model.is_cube else mdata
            chan[il] = m if noconv else convolve_model(m, psf_cube[il])
        chan = np.where(np.isfinite(chan), chan, 0.0) + 1e-10
        # rotate into each frame and write the pure-model _fc cubes
        for n in nums:
            c = read_cube(dirs.reg / frame_filename(n, "reg_cal", prefix))
            rot = rotation_from_header(c.h1) - inst.north_offset_deg
            out = np.empty_like(c.data, dtype=float)
            for il in range(res.nlam):
                out[il] = model_in_frame(
                    chan[il],
                    rot,
                    center,
                    h0=c.h0,
                    dec_deg=dec,
                    lat_deg=inst.latitude_deg,
                    nsc=nsc,
                )
            c.h1["DISKMOD"] = (model.name, "forward-modelled disk model")
            write_cube(dirs.reg / frame_filename(n, "reg_cal_fc", prefix), out, c.h0, c.h1)
        suffix = None
        if rc.rsub > 0:
            imrsub(cfg, prad=True, fc=True, rmax=min(80.0, min(shape) / 2.0), prefix=prefix)
            suffix = "rsub"
        res_path = run_reduction(
            cfg,
            rc,
            reductype=reductype,
            sdi=sdi,
            refpath=refpath,
            fwdmod=True,
            outfile=f"disk_{model.name}.fits",
            suffix=suffix,
        )
        outc = read_cube(res_path)
        imcube = outc.data
        colp = res_path.with_name(res_path.stem + "_collapsed.fits")
        imcol = read_cube(colp).data
        imcol = imcol[0] if imcol.ndim == 3 else imcol
        imcol_con = sumaper_im(imcol, medfwhm * 0.5, lrmax, center)
        if manscale:
            scale = 1.0
        elif simplescale:
            scale = scale_model_ls(refcol_con, imcol_con, roi)
        else:
            scale = scale_model_ls(refcol_con, imcol_con, roi, sigdata)
        imcol_con = imcol_con * scale
        imcol_s = imcol * scale
        input_col = (np.nanmedian(chan, axis=0)) * scale
        cmap = chisq_map(refcol_con, imcol_con, sigdata, roi)
        chisq_red, nbins, cmap_bin = binned_chisq(cmap, medfwhm)
        h1col = refc.h1.copy()
        _model_cards(h1col, model.params, scale, chisq_red, nbins)
        p_sub = dirs.proc / "model_psfsub" / f"{model.name}_psfsubcol.fits"
        p_in = dirs.proc / "model_in" / f"{model.name}_inputcolsc.fits"
        write_cube(p_sub, imcol_s, refc.h0, h1col)
        write_cube(p_in, input_col, refc.h0, h1col)
        write_cube(dirs.proc / "model_stat" / f"{model.name}_chisqmap.fits", cmap, refc.h0, h1col)
        write_cube(
            dirs.proc / "model_stat" / f"{model.name}_chisqmap_bin.fits", cmap_bin, refc.h0, h1col
        )
        write_cube(
            dirs.proc / "model_stat" / f"{model.name}_diff.fits", refcol - imcol_s, refc.h0, h1col
        )
        r = DiskFitResult(model.name, scale, chisq_red, nbins, model.params, p_sub, p_in)
        if calcattenfact:
            sc_l = np.empty(res.nlam)
            cs_l = np.empty(res.nlam)
            cube_s = imcube.copy()
            spec_in = chan.copy()
            for il in range(res.nlam):
                rcon = sumaper_im(refcube[il], medfwhm / 2.0, lrmax, center)
                icon = sumaper_im(imcube[il], medfwhm / 2.0, lrmax, center)
                sc_l[il] = 1.0 if manscale else scale_model_ls(rcon, icon, roi)
                cube_s[il] = imcube[il] * sc_l[il]
                spec_in[il] = chan[il] * sc_l[il]
                d = refcube[il] - cube_s[il]
                sel = roi & np.isfinite(d)
                cs_l[il] = float(np.sum(d[sel] ** 2) / max(sel.sum() - 1, 1))
            h1cube = outc.h1.copy()
            _model_cards(h1cube, model.params, scale, chisq_red, nbins)
            for il in range(res.nlam):
                h1cube[f"SFAC_{il}"] = (sc_l[il], "per-channel model scale")
                h1cube[f"CSQ_{il}"] = (cs_l[il], "per-channel chi^2 in the ROI")
            p_cube = dirs.proc / "model_psfsub" / f"{model.name}_psfsubcube.fits"
            write_cube(p_cube, cube_s, outc.h0, h1cube)
            write_cube(
                dirs.proc / "model_in" / f"{model.name}_inputcubesc.fits", spec_in, outc.h0, h1cube
            )
            with np.errstate(divide="ignore", invalid="ignore"):
                relatten = (spec_in - cube_s) / spec_in
            relatten = np.where(np.abs(relatten) > 1.5, np.nan, relatten)
            p_att = dirs.proc / "model_atten" / f"{model.name}_relattencube.fits"
            write_cube(p_att, relatten, outc.h0, h1cube)
            r.psfsub_cube, r.relatten_cube = p_cube, p_att
            r.scale_per_channel, r.chisq_per_channel = sc_l, cs_l
        with outcomes.open("a") as f:
            f.write(fit_outcomes_line(model.name, model.params, scale, chisq_red) + "\n")
        log.info(
            "disk model %s: scale %.4g, binned reduced chi^2 %.3f over %d bins",
            model.name,
            scale,
            chisq_red,
            nbins,
        )
        results.append(r)
    return results
