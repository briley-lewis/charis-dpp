"""Re-evaluate the χ² of already forward-modelled disks with a different region
of interest or noise scaling (``charis_recalc_disk_fwdmod_chisq.pro``).

Fixes vs. IDL: the procedure actually runs (#191), the SNR routine exists (#192).
"""

from __future__ import annotations

import logging
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

from charis_dpp.config.dataset import DatasetConfig
from charis_dpp.geometry.apertures import sumaper_im
from charis_dpp.geometry.grids import image_center
from charis_dpp.io.fits import read_cube
from charis_dpp.partialsub.disk_fwdmod import reference_noise
from charis_dpp.partialsub.diskfit import (
    MODEL_KEYS,
    chisq_map,
    fit_outcomes_header,
    fit_outcomes_line,
    roi_mask,
)
from charis_dpp.prelims.setupdir import Dirs
from charis_dpp.tools.chisq import binned_chisq

log = logging.getLogger(__name__)

__all__ = ["recalc_disk_chisq"]


def recalc_disk_chisq(
    cfg: DatasetConfig | str | Path,
    *,
    reducname: str = "final.fits",
    models: list[str | Path] | None = None,
    roirange: tuple[float, float, float] | None = None,
    roioffset: tuple[float, float] = (0.0, 0.0),
    lrmin: float = 5.0,
    lrmax: float | None = None,
    snoise: float = 1.0,
) -> list[tuple[str, float, int]]:
    cfg = cfg if isinstance(cfg, DatasetConfig) else DatasetConfig.load(cfg)
    res = cfg.resolve()
    dirs = Dirs(cfg.workdir)
    refc = read_cube(dirs.proc / f"{Path(reducname).stem}_collapsed.fits")
    refcol = refc.data[0] if refc.data.ndim == 3 else refc.data
    shape = refcol.shape
    center = res.star_center or image_center(shape)
    medfwhm = float(np.median(res.fwhm_px))
    lrmax = float(lrmax if lrmax is not None else min(shape) / 2.0)
    refcol_con, sigdata = reference_noise(refcol, medfwhm, lrmax, center)
    sigdata = sigdata / snoise
    roi = roi_mask(
        shape, center=center, lrmin=lrmin, lrmax=lrmax, roirange=roirange, roioffset=roioffset
    )
    paths = (
        [Path(m) for m in models]
        if models
        else sorted((dirs.proc / "model_psfsub").glob("*_psfsubcol.fits"))
    )
    if not paths:
        raise FileNotFoundError("no forward-modelled disks found (run disk-fwdmod first)")
    stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%S")
    outcomes = dirs.proc / "model_stat" / f"NEW_fit_outcomes_{stamp}.dat"
    outcomes.parent.mkdir(parents=True, exist_ok=True)
    lines = [fit_outcomes_header()]
    out = []
    for p in paths:
        c = read_cube(p)
        imcol = c.data[0] if c.data.ndim == 3 else c.data
        params = {k: float(c.h1[k]) for k in MODEL_KEYS if k in c.h1}
        imcol_con = sumaper_im(imcol, medfwhm * 0.5, lrmax, center)
        chisq_red, nbins, _ = binned_chisq(chisq_map(refcol_con, imcol_con, sigdata, roi), medfwhm)
        name = p.name.replace("_psfsubcol.fits", "")
        lines.append(fit_outcomes_line(name, params, float(c.h1.get("SFAC", 1.0)), chisq_red))
        out.append((name, chisq_red, nbins))
        log.info("%s: reduced chi^2 %.3f over %d bins", name, chisq_red, nbins)
    outcomes.write_text("\n".join(lines) + "\n")
    return out
