"""SDI A-LOCI driver (``charis_sdialoci.pro``).

For each frame and each reference channel ``il`` the other channels are
radially rescaled so their speckles align with ``il`` (``align_speckles``),
and LOCI builds a reference for ``il`` from the channels whose speckles
move by more than ``nfwhm·fwhm`` at the annulus radius.  The optimisation
region is the whole annulus (or the annulus minus the subtraction sector
with ``pixmask``), as in the IDL code.

Fixes vs. IDL: the forward-model perturbation uses consistent pixel sets
(#138/#139); empty sectors are skipped (#140); no scratch files (#141).
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np

from charis_dpp.config.dataset import DatasetConfig
from charis_dpp.filt.imrsub import align_speckles
from charis_dpp.geometry.grids import angle_grid, radial_grid
from charis_dpp.psfsub.common import loci_geometry, prepare, ring_edges, write_results
from charis_dpp.psfsub.frame_selection import sdi_references
from charis_dpp.psfsub.solvers.loci import loci_zone

log = logging.getLogger(__name__)

__all__ = ["sdialoci"]


def sdialoci(
    cfg: DatasetConfig | str | Path,
    *,
    nfwhm: float | None = None,
    drsub: float | None = None,
    na: float | None = None,
    geom: float | None = None,
    svd: float | None = None,
    pixmask: bool | None = None,
    nref: int | None = None,
    zero: bool | None = None,
    meanadd: bool | None = None,
    rmin: float | None = None,
    rmax: float | None = None,
    rsubval: int | None = None,
    fc: bool = False,
    fwdmod: bool = False,
    postadi: bool = False,
    norot: bool = False,
    suffix: str | None = None,
    prefix: str = "n",
    outfile: str | None = None,
    channel: int | None = None,
) -> tuple[Path, Path]:
    setup = prepare(
        cfg,
        algorithm="sdialoci",
        rsubval=rsubval,
        fc=fc,
        fwdmod=fwdmod,
        postadi=postadi,
        prefix=prefix,
        suffix=suffix,
        outfile=outfile,
        out_suffix="_sdialocisub",
        norot=norot,
        meanadd=meanadd,
        zero=zero,
    )
    red = setup.cfg.reduction
    lo = red.get("loci", {})
    nfwhm = float(nfwhm if nfwhm is not None else red.get("nfwhm", 0.75)) or 0.75
    drsub = float(drsub if drsub is not None else red.get("dr_px", 5.0))
    rmin = float(rmin if rmin is not None else red.get("rmin_px", 5))
    rmax = float(rmax if rmax is not None else red.get("rmax_px", 65))
    na = float(na if na is not None else lo.get("na", 100.0))
    geom = float(geom if geom is not None else lo.get("geom", 1.0)) or 1.0
    svd_eps = float(svd if svd is not None else 10 ** float(lo.get("svd_log10", -6)))
    nref = int(nref if nref is not None else lo.get("nref", 0))
    pixmask = bool(lo.get("pixmask", False) if pixmask is None else pixmask)
    params = {
        "loci_nfw": nfwhm,
        "loci_na": na,
        "loci_geo": geom,
        "loci_drs": drsub,
        "loci_rmi": rmin,
        "loci_rma": rmax,
        "svd": svd_eps,
        "pixmaskv": int(pixmask),
        "zero": int(setup.zero),
        "meanadd": int(setup.meanadd),
        "sdiztype": 1,
    }
    rings = ring_edges(rmin, rmax, drsub)
    rim_max = min(1.1 * min(setup.shape) / 2.0, rmax + 2 * drsub)
    wav, fwhm = setup.wav, setup.fwhm
    nlam = setup.nlam
    r = radial_grid(setup.shape, setup.center)
    ang = angle_grid(setup.shape, setup.center) % (2 * np.pi)
    resid = np.full((setup.nfiles, nlam) + setup.shape, np.nan, dtype=np.float32)
    from charis_dpp.io.fits import read_cube

    channels = [channel] if channel is not None else range(nlam)
    for nf in range(setup.nfiles):
        cube = read_cube(setup.files[nf]).data
        model = read_cube(setup.files_fwd[nf]).data if fwdmod else None
        for il in channels:
            aligned = align_speckles(cube, wav, il, setup.center)
            aligned = np.where(
                np.isfinite(aligned) | (np.arange(nlam)[:, None, None] == il), aligned, 0.0
            )
            flat = aligned.reshape(nlam, -1).T  # (npix, nλ): "frames" are channels
            flat_m = (
                align_speckles(model, wav, il, setup.center).reshape(nlam, -1).T if fwdmod else None
            )
            out = np.full(flat.shape[0], np.nan)
            for ri, rf in rings:
                nt, r1opt, r2opt = loci_geometry(ri, rf - ri, float(fwhm.max()), na, geom, rim_max)
                ring = (r >= ri) & (r < rf)
                dt = 2 * np.pi / nt
                for it in range(nt):
                    sector = ring & (ang >= it * dt) & (ang < (it + 1) * dt)
                    isub = np.flatnonzero(sector.ravel())
                    iopt = np.flatnonzero((ring & ~sector).ravel() if pixmask else ring.ravel())
                    if isub.size < 5 or iopt.size < 5:
                        continue
                    res = loci_zone(
                        flat[iopt],
                        flat[isub],
                        ref_for=lambda ch, ri=ri, il=il: sdi_references(
                            wav, il, ri, float(fwhm[il]), nfwhm
                        ),
                        nref=nref,
                        svd_eps=svd_eps,
                        zero=setup.zero,
                        opt_model=flat_m[iopt] if fwdmod else None,
                        sub_model=flat_m[isub] if fwdmod else None,
                        min_pix=5,
                    )
                    out[isub] = res[:, il]
            resid[nf, il] = out.reshape(setup.shape)
        log.info("sdialoci: frame %d/%d done", nf + 1, setup.nfiles)
    return write_results(
        setup,
        resid,
        params=params,
        history=f"charis sdialoci nfwhm={nfwhm} dr={drsub} na={na} svd={svd_eps:g}",
    )
