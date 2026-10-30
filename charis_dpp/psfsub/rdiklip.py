"""RDI KLIP driver (``charis_rdiklip_2.pro``): KL basis from a reference-star library.

The reference cubes come from another dataset directory (``refpath``, the
same ``reduc/{reg,rsub}`` layout; all cubes with the chosen suffix).  Zones
are full annuli; science and reference pixels are mean-subtracted
separately; the forward model keeps the over-subtraction term only.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
from astropy.io import fits

from charis_dpp.config.dataset import DatasetConfig
from charis_dpp.geometry.masks import annulus
from charis_dpp.psfsub.common import load_channel, outlier_mask, prepare, ring_edges, write_results
from charis_dpp.psfsub.solvers.klip import klip_basis

log = logging.getLogger(__name__)

__all__ = ["rdiklip"]


def _reference_files(refpath: Path, rsubval: int, suffix: str) -> list[Path]:
    sub = "rsub" if rsubval > 0 else "reg"
    files = sorted((Path(refpath) / "reduc" / sub).glob(f"*{suffix}.fits"))
    if not files:
        raise FileNotFoundError(f"no reference cubes *{suffix}.fits in {refpath}/reduc/{sub}")
    return files


def _load_ref_channel(files: list[Path], il: int, shape) -> np.ndarray:
    out = np.empty((len(files),) + shape)
    for i, p in enumerate(files):
        with fits.open(p, memmap=True) as hdul:
            hdu = hdul[1] if len(hdul) > 1 and hdul[1].data is not None else hdul[0]
            out[i] = np.asarray(hdu.data[il], dtype=float)
    return out


def rdiklip(
    cfg: DatasetConfig | str | Path,
    *,
    refpath: str | Path,
    npca: int | None = None,
    psfrefsuffix: str | None = None,
    drsub: float | None = None,
    meansub: bool | None = None,
    zero: bool | None = None,
    meanadd: bool | None = None,
    rmin: float | None = None,
    rmax: float | None = None,
    rsubval: int | None = None,
    nfwhm: float | None = None,  # unused (kept for CLI symmetry)
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
        algorithm="rdiklip",
        rsubval=rsubval,
        fc=fc,
        fwdmod=fwdmod,
        prefix=prefix,
        suffix=suffix,
        outfile=outfile,
        out_suffix="_rdisub",
        norot=norot,
        meanadd=meanadd,
        zero=zero,
    )
    red = setup.cfg.reduction
    kl = red.get("klip", {})
    npca = int(npca if npca is not None else kl.get("npca", 5))
    meansub = bool(kl.get("meansub", True) if meansub is None else meansub)
    drsub = float(drsub if drsub is not None else red.get("dr_px", 5.0))
    rmin = float(rmin if rmin is not None else red.get("rmin_px", 5))
    rmax = float(
        rmax if rmax is not None else min(red.get("rmax_px", 65), 1.1 * min(setup.shape) / 2)
    )
    sci_suffix = setup.files[0].name[5:-5]  # n####<suffix>.fits
    ref_files = _reference_files(refpath, setup.rsubval, psfrefsuffix or sci_suffix)
    log.info("rdiklip: %d reference cubes from %s", len(ref_files), refpath)
    params = {
        "klip_pca": npca,
        "klip_drs": drsub,
        "klip_rmi": rmin,
        "klip_rma": rmax,
        "klip_nzo": 1,
        "meansub": int(meansub),
        "zero": int(setup.zero),
        "meanadd": int(setup.meanadd),
        "rdi": 1,
        "nrefcub": len(ref_files),
    }
    rings = ring_edges(rmin, rmax, drsub)
    nlam = setup.nlam
    resid = np.full((setup.nfiles, nlam) + setup.shape, np.nan, dtype=np.float32)
    channels = [channel] if channel is not None else range(nlam)
    for il in channels:
        sci = load_channel(setup, il).reshape(setup.nfiles, -1).T
        ref = _load_ref_channel(ref_files, il, setup.shape).reshape(len(ref_files), -1).T
        model = load_channel(setup, il, model=True).reshape(setup.nfiles, -1).T if fwdmod else None
        out = np.full(sci.shape, np.nan)
        for ri, rf in rings:
            isub = np.flatnonzero(annulus(setup.shape, ri, rf, setup.center).ravel())
            if isub.size < 5:
                continue
            S = sci[isub].copy()
            R = ref[isub].copy()
            M = model[isub].copy() if fwdmod else None
            if meansub:
                S -= np.nanmean(S, axis=0, keepdims=True)
                R -= np.nanmean(R, axis=0, keepdims=True)
                if fwdmod:
                    M -= np.nanmean(M, axis=0, keepdims=True)
            good = outlier_mask(S) & np.isfinite(R).all(axis=1)
            if good.sum() < 5 or R.shape[1] < 2:
                continue
            try:
                Z, lam, _ = klip_basis(np.nan_to_num(R[good]), npca)
            except np.linalg.LinAlgError:
                continue
            Z = Z[np.isfinite(lam)]
            zone = np.full(S.shape, np.nan)
            for nf in range(setup.nfiles):
                tgt = np.nan_to_num(S[good, nf]) if not fwdmod else np.nan_to_num(M[good, nf])
                proj = (Z @ tgt) @ Z
                d = np.full(isub.size, np.nan)
                d[good] = (S[good, nf] if not fwdmod else M[good, nf]) - proj
                if setup.zero:
                    d -= np.nanmedian(d)
                zone[:, nf] = d
            out[isub] = zone
        resid[:, il] = out.T.reshape(setup.nfiles, *setup.shape)
        log.info("rdiklip: channel %d/%d done", il + 1, nlam)
    return write_results(
        setup,
        resid,
        params=params,
        history=f"charis rdiklip npca={npca} refs={len(ref_files)} dr={drsub}",
    )
