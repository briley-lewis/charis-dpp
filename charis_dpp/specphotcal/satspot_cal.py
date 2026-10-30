"""Satellite-spot contrast from ladder cubes (``charis_satspot_cal.pro``).

In cubes where the star is unsaturated and the spots are on, the spot/star
flux ratio is measured per channel and fitted with ``P0 (λ/1.55)^-2``.  The
result is written to ``ladder/atten_fac_ladder`` (the IDL version created an
empty file — the fitted value was never written) and used by
``specphot_cal(calmethod=3)``.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from charis_dpp.config.dataset import DatasetConfig
from charis_dpp.geometry.apertures import aper_idl
from charis_dpp.image.filters import radial_profile_2d
from charis_dpp.io.filelists import frame_filename
from charis_dpp.io.fits import get_psf_center, get_sat_positions, read_cube
from charis_dpp.prelims.setupdir import Dirs
from charis_dpp.specphotcal.specphot_cal import measure_spot_fluxes

log = logging.getLogger(__name__)

__all__ = ["satspot_cal"]


def satspot_cal(
    cfg: DatasetConfig | str | Path,
    *,
    calmethod: int = 0,
    files: list[Path] | None = None,
    subskyannulus: bool = False,
    siglim: float | None = None,
    nopradcal: bool = False,
    meancomb: bool = False,
    ap_factor: float | None = None,
    prefix: str = "n",
    suffix: str = "reg",
) -> tuple[float, np.ndarray]:
    """Return ``(P0, ratio_per_channel)`` and write ``ladder/atten_fac_ladder``."""
    cfg = cfg if isinstance(cfg, DatasetConfig) else DatasetConfig.load(cfg)
    res = cfg.resolve()
    dirs = Dirs(cfg.workdir)
    sp = cfg.reduction.get("specphot", {})
    ap_factor = float(ap_factor if ap_factor is not None else sp.get("ap_factor", 1.0))
    siglim = float(siglim if siglim is not None else sp.get("siglim", 5.0))
    ref_um = float(res.instrument.sat_spots.reference_wavelength_um)
    if files is None:
        nums = cfg.frames("ladder")
        if not nums:
            raise ValueError("no ladder frames (observation.frames.ladder)")
        files = [dirs.reg / frame_filename(n, suffix, prefix) for n in nums]
    nlam, wav, fwhm = res.nlam, res.wavelengths_um, res.fwhm_px
    aperrad = 0.5 * fwhm * ap_factor
    cubes, poss, cens = [], [], []
    for p in files:
        c = read_cube(p)
        proc = np.empty_like(c.data)
        for il in range(nlam):
            im = c.data[il]
            if not nopradcal:
                im2 = np.where((im == 0) | ~np.isfinite(im), np.nan, im)
                pr = radial_profile_2d(im2, 1.0, 2 * fwhm[il], min(im.shape) // 2)
                im = np.where(np.isfinite(im) & (im != 0), im - pr, im)
            proc[il] = im
        cubes.append(proc)
        poss.append(get_sat_positions(c.h1, nlam))
        cens.append(get_psf_center(c.h1) or (c.shape[1] // 2, c.shape[0] // 2))
    if calmethod == 1 and len(cubes) > 1:
        cubes = [
            np.nanmean(np.stack(cubes), axis=0)
            if meancomb
            else np.nanmedian(np.stack(cubes), axis=0)
        ]
        poss = [np.nanmedian(np.stack(poss), axis=0)]
        cens = [tuple(np.nanmedian(np.array(cens), axis=0))]
    sat, esat, star = [], [], []
    for cube, pos, cen in zip(cubes, poss, cens, strict=True):
        ph = measure_spot_fluxes(
            cube,
            pos,
            aperrad,
            subskyannulus=subskyannulus,
            meancomb=True,
            expandederror=True,
            center=cen,
        )
        sat.append(np.nanmean(ph.per_spot, axis=1))
        esat.append(ph.err)
        st = np.array(
            [
                aper_idl(
                    cube[il],
                    cen[0],
                    cen[1],
                    aperrad[il],
                    skyrad=(2 * aperrad[il], 6 * aperrad[il]) if subskyannulus else None,
                    setskyval=None if subskyannulus else 0.0,
                )[0]
                for il in range(nlam)
            ]
        )
        star.append(st)
    satf = np.nanmedian(sat, axis=0)
    esatf = np.nanmedian(esat, axis=0)
    starf = np.nanmedian(star, axis=0)
    ratio = satf / starf
    eratio = esatf / starf
    good = (esatf / satf <= 1.0 / siglim) & np.isfinite(ratio) & (starf > 0)
    if good.sum() < 2:
        raise RuntimeError("too few channels with significant spot detections for the contrast fit")
    p0 = ratio[np.argmin(np.abs(wav - ref_um))]

    def resid(p):
        return (p[0] * (wav[good] / ref_um) ** -2 - ratio[good]) / np.maximum(eratio[good], 1e-12)

    fit = least_squares(resid, [p0], method="lm")
    P0 = float(fit.x[0])
    dirs.ladder.mkdir(exist_ok=True)
    out = dirs.ladder / "atten_fac_ladder"
    lines = [
        f"# P0 (spot/star at {ref_um} um), fit of P0*(lambda/{ref_um})^-2",
        f"{P0:.8e}",
        "# channel wavelength_um ratio err",
    ]
    lines += [f"{il:3d} {wav[il]:.5f} {ratio[il]:.6e} {eratio[il]:.6e}" for il in range(nlam)]
    out.write_text("\n".join(lines) + "\n")
    log.info(
        "spot/star contrast at %.2f um: %.4e (%.2f mag) -> %s", ref_um, P0, -2.5 * np.log10(P0), out
    )
    return P0, ratio
