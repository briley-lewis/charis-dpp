"""ADI A-LOCI driver (``charis_adialoci.pro``)."""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np

from charis_dpp.config.dataset import DatasetConfig
from charis_dpp.psfsub.common import (
    load_channel,
    loci_geometry,
    prepare,
    ring_edges,
    write_results,
    zone_masks,
)
from charis_dpp.psfsub.frame_selection import adi_references
from charis_dpp.psfsub.solvers.loci import loci_zone

log = logging.getLogger(__name__)

__all__ = ["adialoci"]


def adialoci(
    cfg: DatasetConfig | str | Path,
    *,
    nfwhm: float | None = None,
    drsub=None,
    na=None,
    geom: float | None = None,
    svd: float | None = None,
    pixmask: bool | None = None,
    nref: int | None = None,
    zero: bool | None = None,
    meanadd: bool | None = None,
    rmin: float | None = None,
    rmax: float | None = None,
    rsubval: int | None = None,
    channel: int | None = None,
    fc: bool = False,
    fwdmod: bool = False,
    postsdi: bool = False,
    norot: bool = False,
    suffix: str | None = None,
    prefix: str = "n",
    outfile: str | None = None,
) -> tuple[Path, Path]:
    setup = prepare(
        cfg,
        algorithm="aloci",
        rsubval=rsubval,
        fc=fc,
        fwdmod=fwdmod,
        postsdi=postsdi,
        prefix=prefix,
        suffix=suffix,
        outfile=outfile,
        out_suffix="_alocisub",
        norot=norot,
        meanadd=meanadd,
        zero=zero,
    )
    red = setup.cfg.reduction
    lo = red.get("loci", {})
    nfwhm = float(nfwhm if nfwhm is not None else red.get("nfwhm", 0.75)) or 0.75
    drsub = drsub if drsub is not None else red.get("dr_px", 5.0)
    rmin = float(rmin if rmin is not None else red.get("rmin_px", 5))
    rmax = float(rmax if rmax is not None else red.get("rmax_px", 65))
    na = na if na is not None else lo.get("na", 100.0)
    geom = float(geom if geom is not None else lo.get("geom", 1.0)) or 1.0
    svd_eps = float(svd if svd is not None else 10 ** float(lo.get("svd_log10", -6)))
    nref = int(nref if nref is not None else lo.get("nref", 0))
    pixmask = bool(lo.get("pixmask", False) if pixmask is None else pixmask)
    params = {
        "loci_nfw": nfwhm,
        "loci_na": float(np.atleast_1d(na)[0]),
        "loci_geo": geom,
        "loci_drs": float(np.atleast_1d(drsub)[0]),
        "loci_rmi": rmin,
        "loci_rma": rmax,
        "loci_nre": nref,
        "svd": svd_eps,
        "pixmaskv": int(pixmask),
        "zero": int(setup.zero),
        "meanadd": int(setup.meanadd),
    }
    rings = ring_edges(rmin, rmax, drsub)
    rim_max = min(
        1.1 * min(setup.shape) / 2.0, max(r for _, r in rings) + 2 * float(np.atleast_1d(drsub)[0])
    )
    nlam = setup.nlam
    resid = np.full((setup.nfiles, nlam) + setup.shape, np.nan, dtype=np.float32)
    channels = [channel] if channel is not None else range(nlam)
    for il in channels:
        fwhm = float(setup.fwhm[il])
        data = load_channel(setup, il)
        model = load_channel(setup, il, model=True) if fwdmod else None
        flat = data.reshape(setup.nfiles, -1).T  # (npix, nf)
        flat_m = model.reshape(setup.nfiles, -1).T if fwdmod else None
        out = np.full(flat.shape, np.nan)
        for ri, rf in rings:
            nt, r1opt, r2opt = loci_geometry(ri, rf - ri, fwhm, na, geom, rim_max)
            for it in range(nt):
                sub_m, opt_m = zone_masks(
                    setup.shape, setup.center, ri, rf, r1opt, r2opt, nt, it, pixmask=pixmask
                )
                isub = np.flatnonzero(sub_m.ravel())
                iopt = np.flatnonzero(opt_m.ravel())
                if isub.size < 2 or iopt.size < 5:
                    continue
                res = loci_zone(
                    flat[iopt],
                    flat[isub],
                    ref_for=lambda nf, ri=ri, fwhm=fwhm: adi_references(
                        setup.pa, setup.dtpose, nf, ri, fwhm, nfwhm
                    ),
                    nref=nref,
                    svd_eps=svd_eps,
                    zero=setup.zero,
                    opt_model=flat_m[iopt] if fwdmod else None,
                    sub_model=flat_m[isub] if fwdmod else None,
                )
                out[isub] = res
        resid[:, il] = out.T.reshape(setup.nfiles, *setup.shape)
        log.info("adialoci: channel %d/%d done (%d rings)", il + 1, nlam, len(rings))
    return write_results(
        setup,
        resid,
        params=params,
        history=f"charis adialoci nfwhm={nfwhm} dr={drsub} na={na} svd={svd_eps:g} nref={nref}",
    )
