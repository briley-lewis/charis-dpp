"""Negative-planet injection (``charis_{aloci,klip}_negspec_planet.pro``).

The extracted planet spectrum (``spectrum.dat`` columns wavelength, flux)
is injected with a negative sign at the planet position into the real
data and the reduction is re-run; a good spectrum/position nulls the
planet in ``reduc/proc/res0_np.fits``.  Fixes vs. IDL: the output that is
read back is the one that was written (#172); no ``fwdmod`` in the chains (#173).
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
    run_reduction,
)
from charis_dpp.partialsub.inject import Planet
from charis_dpp.prelims.setupdir import Dirs

log = logging.getLogger(__name__)

__all__ = ["negative_planet"]


def negative_planet(
    cfg: DatasetConfig | str | Path,
    *,
    reducname: str = "final.fits",
    planet: Planet | tuple[float, float],
    spectrum_file: str | Path = "spectrum.dat",
    psf: str | Path = "empirical",
    reductype: int = 0,
    sdi_reducname: str | None = None,
    refpath: str | Path | None = None,
    outfile: str = "res0_np.fits",
    prefix: str = "n",
) -> np.ndarray:
    cfg = cfg if isinstance(cfg, DatasetConfig) else DatasetConfig.load(cfg)
    res = cfg.resolve()
    dirs = Dirs(cfg.workdir)
    rc = ReducedCube.load(dirs.proc / reducname, res.nlam)
    sdi = ReducedCube.load(dirs.proc / sdi_reducname, res.nlam) if sdi_reducname else None
    sp = Path(spectrum_file) if Path(spectrum_file).exists() else dirs.proc / str(spectrum_file)
    arr = np.loadtxt(sp, ndmin=2)
    spec = np.interp(res.wavelengths_um, arr[:, 0], arr[:, 1])
    pl = planet if isinstance(planet, Planet) else Planet(float(planet[0]), float(planet[1]))
    pl.spectrum = spec
    cal = read_cube(dirs.reg / f"{prefix}{cfg.frames('science')[0]:04d}reg_cal.fits")
    r_ap = (
        rc.r_ap
        if np.isfinite(rc.r_ap).all()
        else np.array([float(cal.h1[f"R_AP{l}"]) for l in range(res.nlam)])
    )
    psf_cube = build_psf(cfg, res, r_ap, kind=psf)
    inject_planets_into_sequence(
        cfg, res, [pl], psf_cube, suffix="reg_cal", prefix=prefix, into_data=True, negative=True
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
    cube = read_cube(res_path).data
    residual = np.array(
        [aper_idl(cube[l], pl.x, pl.y, float(r_ap[l]), setskyval=0.0)[0] for l in range(res.nlam)]
    )
    np.savetxt(
        dirs.proc / "negplanet_residual.txt",
        np.column_stack([res.wavelengths_um, residual, spec]),
        fmt="%.6e",
        header="wavelength_um residual_flux injected_flux",
    )
    log.info("negative planet: median residual/injected = %.3f", np.nanmedian(residual / spec))
    return residual
