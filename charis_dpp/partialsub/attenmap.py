"""Attenuation (throughput) map from a grid of forward-modelled planets
(``charis_{aloci,klip}_attenmap_planet.pro``).

For ``ntc`` angular offsets, planets are placed on a spiral at radii
``lrmin..lrmax`` (step ``dlrstep``) and forward-modelled together; the
per-planet, per-channel throughput is azimuthally averaged into a radial
profile and expanded to 2-D.  Outputs ``reduc/proc/attenmap_lambda.fits``
``(nλ, ny, nx)`` and ``attenmap.fits`` (median over channels).

Fixes vs. IDL: all ``ntc`` offsets are kept (#169) and the per-planet results
are all written (#168).
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline

from charis_dpp.config.dataset import DatasetConfig
from charis_dpp.geometry.grids import image_center, radial_grid
from charis_dpp.io.fits import read_cube, write_cube, write_image
from charis_dpp.partialsub.fwdmod_planet import fwdmod_planet
from charis_dpp.partialsub.inject import Planet
from charis_dpp.prelims.setupdir import Dirs

log = logging.getLogger(__name__)

__all__ = ["attenuation_map", "planet_grid"]


def planet_grid(
    center, lrmin: float, lrmax: float, dlrstep: float, itc: int, ntc: int, log_contrast: float
) -> list[Planet]:
    radii = np.arange(lrmin, lrmax + 1e-6, dlrstep)
    n = radii.size
    planets = []
    for irc, rc in enumerate(radii):
        ang = np.deg2rad(itc * 360.0 / ntc / max(n, 1) + irc * 360.0 / n + 180.0 * (irc % 4))
        planets.append(
            Planet(center[0] + rc * np.cos(ang), center[1] + rc * np.sin(ang), log_contrast)
        )
    return planets


def attenuation_map(
    cfg: DatasetConfig | str | Path,
    *,
    reducname: str = "final.fits",
    contrast: float | None = None,
    ntc: int | None = None,
    lrmin: float | None = None,
    lrmax: float | None = None,
    dlrstep: float | None = None,
    psf: str | Path = "empirical",
    planet_model: str | Path | None = None,
    reductype: int = 0,
    sdi_reducname: str | None = None,
    refpath: str | Path | None = None,
) -> tuple[Path, Path]:
    cfg = cfg if isinstance(cfg, DatasetConfig) else DatasetConfig.load(cfg)
    res = cfg.resolve()
    dirs = Dirs(cfg.workdir)
    fm = cfg.reduction.get("fwdmod", {})
    contrast = float(contrast if contrast is not None else fm.get("contrast_log10", -6.0))
    ntc = int(ntc if ntc is not None else fm.get("ntc", 5))
    lrmin = float(lrmin if lrmin is not None else fm.get("lrmin_px", 10.0))
    lrmax = float(lrmax if lrmax is not None else fm.get("lrmax_px", 60.0))
    dlrstep = float(dlrstep if dlrstep is not None else fm.get("dlrstep_px", 2.5))
    first = read_cube(dirs.proc / reducname)
    shape = first.shape
    center = res.star_center or image_center(shape)
    lrmax = min(lrmax, min(shape) / 2.0 - 3 * float(res.fwhm_px.max()))
    nlam = res.nlam
    xs, ys, vals = [], [], []  # per planet: position, throughput (nλ)
    for itc in range(ntc):
        planets = planet_grid(center, lrmin, lrmax, dlrstep, itc, ntc, contrast)
        thr = fwdmod_planet(
            cfg,
            reducname=reducname,
            planets=planets,
            psf=psf,
            planet_model=planet_model,
            reductype=reductype,
            sdi_reducname=sdi_reducname,
            refpath=refpath,
            outfile=f"attenres{itc}.fits",
            write_files=False,
        )
        for pl, t in zip(planets, thr, strict=True):
            xs.append(pl.x)
            ys.append(pl.y)
            vals.append(t)
        np.savetxt(
            dirs.proc / f"atten{itc}.dat",
            np.column_stack([[p.x for p in planets], [p.y for p in planets], thr]),
            fmt="%.5e",
            header="x y throughput per channel",
        )
        log.info("attenuation map: offset %d/%d done", itc + 1, ntc)
    xs, ys, vals = np.array(xs), np.array(ys), np.array(vals)  # vals (nplanet, nλ)
    rp = np.hypot(xs - center[0], ys - center[1])
    r = radial_grid(shape, center)
    att_l = np.full((nlam,) + shape, np.nan)
    radii = np.unique(np.round(rp / dlrstep) * dlrstep)
    for il in range(nlam):
        prof = np.array([np.nanmedian(vals[np.abs(rp - rr) < dlrstep / 2, il]) for rr in radii])
        good = np.isfinite(prof)
        if good.sum() < 2:
            continue
        spl = CubicSpline(radii[good], prof[good], bc_type="natural")
        att_l[il] = np.clip(spl(np.clip(r, radii[good][0], radii[good][-1])), 0.0, 1.5)
    att = np.nanmedian(att_l, axis=0)
    p_l = dirs.proc / "attenmap_lambda.fits"
    p_c = dirs.proc / "attenmap.fits"
    write_cube(p_l, att_l, first.h0, first.h1)
    write_image(p_c, att)
    log.info("wrote %s and %s", p_l, p_c)
    return p_l, p_c
