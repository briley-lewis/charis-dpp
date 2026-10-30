"""ADI KLIP driver (``charis_adiklip.pro``)."""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np

from charis_dpp.config.dataset import DatasetConfig
from charis_dpp.psfsub.common import load_channel, prepare, ring_edges, write_results, zone_masks
from charis_dpp.psfsub.frame_selection import adi_references
from charis_dpp.psfsub.solvers.klip import klip_zone

log = logging.getLogger(__name__)

__all__ = ["adiklip"]


def adiklip(
    cfg: DatasetConfig | str | Path,
    *,
    npca: int | None = None,
    nzone: int | None = None,
    nfwhm: float | None = None,
    drsub=None,
    meansub: bool | None = None,
    zero: bool | None = None,
    meanadd: bool | None = None,
    rmin: float | None = None,
    rmax: float | None = None,
    rsubval: int | None = None,
    channel: int | None = None,
    fc: bool = False,
    fwdmod: bool = False,
    norot: bool = False,
    suffix: str | None = None,
    prefix: str = "n",
    outfile: str | None = None,
) -> tuple[Path, Path]:
    setup = prepare(
        cfg,
        algorithm="klip",
        rsubval=rsubval,
        fc=fc,
        fwdmod=fwdmod,
        prefix=prefix,
        suffix=suffix,
        outfile=outfile,
        out_suffix="_adisub",
        norot=norot,
        meanadd=meanadd,
        zero=zero,
    )
    red = setup.cfg.reduction
    kl = red.get("klip", {})
    npca = int(npca if npca is not None else kl.get("npca", 5))
    nzone = int(nzone if nzone is not None else kl.get("nzone", 1))
    meansub = bool(kl.get("meansub", True) if meansub is None else meansub)
    nfwhm = float(nfwhm if nfwhm is not None else red.get("nfwhm", 0.75)) or 0.75
    drsub = drsub if drsub is not None else red.get("dr_px", 5.0)
    rmin = float(rmin if rmin is not None else red.get("rmin_px", 5))
    rmax = float(rmax if rmax is not None else red.get("rmax_px", 65))
    params = {
        "klip_pca": npca,
        "klip_nfw": nfwhm,
        "klip_drs": float(np.atleast_1d(drsub)[0]),
        "klip_rmi": rmin,
        "klip_rma": rmax,
        "klip_nzo": nzone,
        "meansub": int(meansub),
        "zero": int(setup.zero),
        "meanadd": int(setup.meanadd),
    }
    rings = ring_edges(rmin, rmax, drsub)
    nlam = setup.nlam
    resid = np.full((setup.nfiles, nlam) + setup.shape, np.nan, dtype=np.float32)
    channels = [channel] if channel is not None else range(nlam)
    for il in channels:
        fwhm = float(setup.fwhm[il])
        data = load_channel(setup, il)
        model = load_channel(setup, il, model=True) if fwdmod else None
        flat = data.reshape(setup.nfiles, -1).T
        flat_m = model.reshape(setup.nfiles, -1).T if fwdmod else None
        out = np.full(flat.shape, np.nan)
        for ri, rf in rings:
            for it in range(nzone):
                sub_m, _ = zone_masks(
                    setup.shape, setup.center, ri, rf, ri, rf, nzone, it, pixmask=False
                )
                isub = np.flatnonzero(sub_m.ravel())
                if isub.size < 5:
                    continue
                out[isub] = klip_zone(
                    flat[isub],
                    ref_for=lambda nf, ri=ri, fwhm=fwhm: adi_references(
                        setup.pa, setup.dtpose, nf, ri, fwhm, nfwhm
                    ),
                    npca=npca,
                    meansub=meansub,
                    zero=setup.zero,
                    sub_model=flat_m[isub] if fwdmod else None,
                )
        resid[:, il] = out.T.reshape(setup.nfiles, *setup.shape)
        log.info("adiklip: channel %d/%d done", il + 1, nlam)
    return write_results(
        setup,
        resid,
        params=params,
        history=f"charis adiklip npca={npca} nzone={nzone} nfwhm={nfwhm} dr={drsub}",
    )
