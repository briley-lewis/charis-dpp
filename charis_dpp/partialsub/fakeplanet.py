"""Inject fake planets into the real data and re-reduce (``charis_{aloci,klip}_fakeplanet.pro``).

Unlike :mod:`fwdmod_planet` the planets ride on the real speckles, so the
recovered flux includes residual speckle noise.  Returns the recovered /
injected flux ratio per planet and channel and writes ``fakeplanet_res.txt``.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np

from charis_dpp.config.dataset import DatasetConfig
from charis_dpp.filt.imrsub import imrsub
from charis_dpp.io.fits import read_cube
from charis_dpp.partialsub.common import (
    ReducedCube,
    build_psf,
    inject_planets_into_sequence,
    read_targfile,
    run_reduction,
)
from charis_dpp.partialsub.fwdmod_planet import measure_throughput
from charis_dpp.partialsub.inject import Planet
from charis_dpp.prelims.setupdir import Dirs

log = logging.getLogger(__name__)

__all__ = ["fake_planet"]


def fake_planet(
    cfg: DatasetConfig | str | Path,
    *,
    reducname: str = "final.fits",
    planets: list[Planet] | None = None,
    targfile: str | Path | None = None,
    psf: str | Path = "empirical",
    planet_model: str | Path | None = None,
    reductype: int = 0,
    sdi_reducname: str | None = None,
    refpath: str | Path | None = None,
    outfile: str = "fake0.fits",
    prefix: str = "n",
) -> np.ndarray:
    cfg = cfg if isinstance(cfg, DatasetConfig) else DatasetConfig.load(cfg)
    res = cfg.resolve()
    dirs = Dirs(cfg.workdir)
    if planets is None:
        planets = read_targfile(targfile)
    rc = ReducedCube.load(dirs.proc / reducname, res.nlam)
    sdi = ReducedCube.load(dirs.proc / sdi_reducname, res.nlam) if sdi_reducname else None
    cal = read_cube(dirs.reg / f"{prefix}{cfg.frames('science')[0]:04d}reg_cal.fits")
    r_ap = (
        rc.r_ap
        if np.isfinite(rc.r_ap).all()
        else np.array([float(cal.h1[f"R_AP{l}"]) for l in range(res.nlam)])
    )
    psf_cube = build_psf(cfg, res, r_ap, kind=psf)
    input_spec = inject_planets_into_sequence(
        cfg,
        res,
        planets,
        psf_cube,
        suffix="reg_cal",
        prefix=prefix,
        into_data=True,
        model=planet_model,
    )
    suffix = None
    if rc.rsub > 0:
        imrsub(cfg, prad=True, fc=True, rmax=min(80.0, min(cal.shape) / 2.0), prefix=prefix)
        suffix = "rsub"
    res_path = run_reduction(
        cfg,
        rc,
        reductype=reductype,
        sdi=sdi,
        refpath=refpath,
        fwdmod=False,
        outfile=outfile,
        suffix=suffix,
    )
    thr = measure_throughput(read_cube(res_path).data, planets, input_spec, r_ap)
    np.savetxt(
        dirs.proc / "fakeplanet_res.txt",
        np.column_stack([res.wavelengths_um, thr.T]),
        fmt="%.6e",
        header="wavelength_um recovered/injected per planet",
    )
    return thr
