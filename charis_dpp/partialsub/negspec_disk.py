"""Negative-disk injection (``charis_aloci_negspec_disk.pro``): subtract a
(forward-modelled, scaled) disk model cube from every frame and re-reduce;
a good model nulls the disk in ``reduc/proc/<outfile>``.

Fixes vs. IDL: the reduction is called with ``fc`` instead of the
non-existent ``np`` keyword (#188), the NaN and zero masks are separate (#189),
no per-frame ``psf.fits`` scratch file (#190), one rotation convention (#187).
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np

from charis_dpp.config.dataset import DatasetConfig
from charis_dpp.filt.imrsub import imrsub
from charis_dpp.geometry.grids import image_center
from charis_dpp.image.wcs import rotation_from_header
from charis_dpp.io.filelists import frame_filename
from charis_dpp.io.fits import read_cube, write_cube
from charis_dpp.partialsub.common import ReducedCube, run_reduction
from charis_dpp.partialsub.diskfit import load_disk_model, model_in_frame
from charis_dpp.prelims.setupdir import Dirs

log = logging.getLogger(__name__)

__all__ = ["negative_disk"]


def negative_disk(
    cfg: DatasetConfig | str | Path,
    diskmodel: str | Path,
    *,
    reducname: str = "final.fits",
    reductype: int = 0,
    sdi_reducname: str | None = None,
    refpath: str | Path | None = None,
    outfile: str = "res0_np.fits",
    nsc: int | None = None,
    prefix: str = "n",
) -> Path:
    cfg = cfg if isinstance(cfg, DatasetConfig) else DatasetConfig.load(cfg)
    res = cfg.resolve()
    dirs = Dirs(cfg.workdir)
    inst = res.instrument
    nsc = int(nsc if nsc is not None else cfg.reduction.get("fwdmod", {}).get("nsc", 21))
    rc = ReducedCube.load(dirs.proc / reducname, res.nlam)
    sdi = ReducedCube.load(dirs.proc / sdi_reducname, res.nlam) if sdi_reducname else None
    model = load_disk_model(diskmodel)
    if model.is_cube and model.data.shape[0] != res.nlam:
        raise ValueError(f"disk model cube has {model.data.shape[0]} channels, data {res.nlam}")
    dec = float(cfg.target.get("dec_deg", 0.0))
    shape = None
    for n in cfg.frames("science"):
        c = read_cube(dirs.reg / frame_filename(n, "reg_cal", prefix))
        shape = c.shape
        center = res.star_center or image_center(c.shape)
        rot = rotation_from_header(c.h1) - inst.north_offset_deg
        out = c.data.astype(float).copy()
        for il in range(res.nlam):
            m = model_in_frame(
                model.channel(il),
                rot,
                center,
                h0=c.h0,
                dec_deg=dec,
                lat_deg=inst.latitude_deg,
                nsc=nsc,
            )
            fov = np.isfinite(out[il]) & (out[il] != 0)  # #189: keep the data's own footprint
            out[il][fov] -= m[fov]
        c.h1["NEGDISK"] = (model.name, "negative disk model injected")
        write_cube(dirs.reg / frame_filename(n, "reg_cal_fc", prefix), out, c.h0, c.h1)
    suffix = None
    if rc.rsub > 0:
        imrsub(cfg, prad=True, fc=True, rmax=min(80.0, min(shape) / 2.0), prefix=prefix)
        suffix = "rsub"
    path = run_reduction(
        cfg,
        rc,
        reductype=reductype,
        sdi=sdi,
        refpath=refpath,
        fwdmod=False,
        outfile=outfile,
        suffix=suffix,
    )
    log.info("negative-disk reduction written to %s", path)
    return path
