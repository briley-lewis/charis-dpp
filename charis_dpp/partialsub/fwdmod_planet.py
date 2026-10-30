"""Forward-model a planet through the PSF subtraction to measure its throughput
(``charis_{aloci,klip}_fwdmod_planet.pro``).

The model planet(s) are inserted into empty cubes, run through the same
reduction (same parameters, same reference frames) with ``fwdmod``, and the
per-channel aperture flux of the processed model divided by the injected
flux is the throughput, written to ``reduc/proc/synth_throughput<i>.txt``
(columns: wavelength [µm], throughput).

Fixes vs. IDL: all chain options are reachable (#162), intermediate outputs
no longer clobber ``final_fc.fits`` (#163), one position convention (#187).
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np

from charis_dpp.config.dataset import DatasetConfig
from charis_dpp.filt.imrsub import imrsub
from charis_dpp.geometry.apertures import aper_idl
from charis_dpp.io.fits import read_cube
from charis_dpp.partialsub.common import (
    ReducedCube,
    build_psf,
    inject_planets_into_sequence,
    read_targfile,
    run_reduction,
)
from charis_dpp.partialsub.inject import Planet
from charis_dpp.prelims.setupdir import Dirs

log = logging.getLogger(__name__)

__all__ = ["fwdmod_planet", "measure_throughput"]


def measure_throughput(
    res_cube: np.ndarray, planets: list[Planet], input_spec: np.ndarray, r_ap: np.ndarray
) -> np.ndarray:
    """``(nplanet, nλ)`` ratio of recovered aperture flux to injected flux."""
    out = np.full((len(planets), res_cube.shape[0]), np.nan)
    for j, pl in enumerate(planets):
        for l in range(res_cube.shape[0]):
            flux = aper_idl(res_cube[l], pl.x, pl.y, float(r_ap[l]), setskyval=0.0)[0]
            out[j, l] = flux / input_spec[j, l] if input_spec[j, l] != 0 else np.nan
    return out


def fwdmod_planet(
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
    empspectrum: np.ndarray | None = None,
    outfile: str = "res0.fits",
    prefix: str = "n",
    write_files: bool = True,
) -> np.ndarray:
    cfg = cfg if isinstance(cfg, DatasetConfig) else DatasetConfig.load(cfg)
    res = cfg.resolve()
    dirs = Dirs(cfg.workdir)
    if planets is None:
        if targfile is None:
            raise ValueError("give planets=[Planet(...)] or targfile=")
        planets = read_targfile(targfile)
    rc = ReducedCube.load(
        dirs.proc / reducname if not Path(reducname).is_absolute() else reducname, res.nlam
    )
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
        model=planet_model,
        empspectrum=empspectrum,
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
        fwdmod=True,
        outfile=outfile,
        suffix=suffix,
    )
    res_cube = read_cube(res_path).data
    thr = measure_throughput(res_cube, planets, input_spec, r_ap)
    if write_files:
        for j in range(len(planets)):
            np.savetxt(
                dirs.proc / f"synth_throughput{j}.txt",
                np.column_stack([res.wavelengths_um, thr[j]]),
                fmt="%.6f %.6e",
                header="wavelength_um throughput",
            )
    for j, pl in enumerate(planets):
        log.info(
            "planet %d at (%.1f, %.1f): median throughput %.3f", j, pl.x, pl.y, np.nanmedian(thr[j])
        )
    return thr
