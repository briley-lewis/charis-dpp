"""Shared pieces of the forward-modelling drivers (the ~250 lines each IDL file repeated)."""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from astropy.io.fits import Header

from charis_dpp.config.dataset import DatasetConfig, Resolved
from charis_dpp.config.paths import models_dir
from charis_dpp.geometry.grids import image_center
from charis_dpp.image.wcs import rotation_from_header
from charis_dpp.io.filelists import frame_filename
from charis_dpp.io.fits import Cube, read_cube, write_cube
from charis_dpp.partialsub.inject import Planet, frame_position, insert_planet, normalise_psf
from charis_dpp.planetspec.planet_photcal import planet_spectrum_mjy, star_spectrum_from_header
from charis_dpp.prelims.setupdir import Dirs
from charis_dpp.specphotcal.synth_psf import gaussian_psf_cube

log = logging.getLogger(__name__)

__all__ = [
    "ReducedCube",
    "build_psf",
    "inject_planets_into_sequence",
    "run_reduction",
    "read_targfile",
    "planet_model_path",
]


@dataclass
class ReducedCube:
    path: Path
    cube: Cube
    algorithm: str  # 'aloci' | 'klip' | 'sdialoci' | 'rdiklip'
    params: dict = field(default_factory=dict)
    rotang: float = 0.0
    r_ap: np.ndarray | None = None
    fstar: np.ndarray | None = None
    rsub: int = 0

    @classmethod
    def load(cls, path: str | Path, nlam: int | None = None) -> ReducedCube:
        c = read_cube(path)
        h = c.h1
        nlam = nlam or c.nlam
        algo = str(h.get("PSFALGO", "")).strip()
        if not algo:
            algo = "aloci" if "LOCI_NFW" in h else ("klip" if "KLIP_PCA" in h else "unknown")
        p: dict = {}
        if algo in ("aloci", "sdialoci"):
            p = {
                "nfwhm": float(h.get("LOCI_NFW", 0.75)),
                "na": float(h.get("LOCI_NA", 100)),
                "geom": float(h.get("LOCI_GEO", 1.0)),
                "drsub": float(h.get("LOCI_DRS", 5)),
                "rmin": float(h.get("LOCI_RMI", 5)),
                "rmax": float(h.get("LOCI_RMA", 65)),
                "nref": int(h.get("LOCI_NRE", 0)),
                "svd": float(h.get("SVD", 1e-99)),
                "pixmask": bool(int(h.get("PIXMASKV", 0))),
                "zero": bool(int(h.get("ZERO", 0))),
                "meanadd": bool(int(h.get("MEANADD", 0))),
            }
        elif algo in ("klip", "rdiklip"):
            p = {
                "npca": int(h.get("KLIP_PCA", 5)),
                "nfwhm": float(h.get("KLIP_NFW", 0.75)),
                "drsub": float(h.get("KLIP_DRS", 5)),
                "rmin": float(h.get("KLIP_RMI", 5)),
                "rmax": float(h.get("KLIP_RMA", 65)),
                "nzone": int(h.get("KLIP_NZO", 1)),
                "meansub": bool(int(h.get("MEANSUB", 1))),
                "zero": bool(int(h.get("ZERO", 0))),
                "meanadd": bool(int(h.get("MEANADD", 0))),
            }
        r_ap = np.array([float(h.get(f"R_AP{l}", np.nan)) for l in range(nlam)])
        return cls(
            Path(path),
            c,
            algo,
            p,
            float(h.get("ROTANG", 0.0)),
            r_ap,
            star_spectrum_from_header(h, nlam),
            int(h.get("RSUB", 0)),
        )


def build_psf(
    cfg: DatasetConfig,
    res: Resolved,
    r_ap: np.ndarray,
    *,
    kind: str | Path = "empirical",
    dpsf: int = 21,
) -> np.ndarray:
    """Intrinsic PSF cube ``(nλ, dpsf, dpsf)`` normalised to unit flux in ``r_ap`` per channel."""
    dirs = Dirs(cfg.workdir)
    if isinstance(kind, str | Path) and str(kind) not in ("empirical", "gaussian"):
        psf = read_cube(kind).data
    elif str(kind) == "gaussian":
        psf = gaussian_psf_cube(res.fwhm_px, dpsf)
    else:
        p = dirs.psfmodel / "psfcube_med.fits"
        if not p.exists():
            raise FileNotFoundError(f"{p} not found: run make-emp-psf (or pass psf='gaussian')")
        psf = read_cube(p).data
    if psf.shape[-1] != dpsf:
        from charis_dpp.image.interp import subarr

        psf = np.stack([subarr(s, dpsf, fill=0.0) for s in psf])
    if psf.shape[0] != res.nlam:
        raise ValueError(f"PSF cube has {psf.shape[0]} channels, data {res.nlam}")
    return normalise_psf(psf, np.where(np.isfinite(r_ap), r_ap, 0.5 * res.fwhm_px))


def planet_model_path(cfg: DatasetConfig, name: str | Path | None) -> Path:
    name = Path(name or cfg.reduction.get("planet_model", "L0_g.fits"))
    if name.exists():
        return name
    cand = models_dir() / "planetmodels" / name
    if cand.exists():
        return cand
    raise FileNotFoundError(
        f"planet model {name} not found (looked in {models_dir() / 'planetmodels'})"
    )


def read_targfile(path: str | Path) -> list[Planet]:
    """``x y log_contrast`` per line (the order ``charis_check_for_detection`` writes)."""
    arr = np.atleast_2d(np.loadtxt(path, ndmin=2))
    return [Planet(float(r[0]), float(r[1]), float(r[2]) if r.size > 2 else -4.0) for r in arr]


def inject_planets_into_sequence(
    cfg: DatasetConfig,
    res: Resolved,
    planets: list[Planet],
    psf: np.ndarray,
    *,
    suffix: str,
    prefix: str = "n",
    into_data: bool = False,
    negative: bool = False,
    model: str | Path | None = None,
    empspectrum: np.ndarray | None = None,
    nsc: int | None = None,
) -> np.ndarray:
    """Write ``n####<suffix>_fc.fits`` for every science frame; returns input spectra ``(nplanet, nλ)`` (median over frames).

    With ``into_data`` the planets are added to the real cubes (fake planets);
    otherwise to an empty cube (pure model for forward modelling).
    """
    dirs = Dirs(cfg.workdir)
    inst = res.instrument
    nums = cfg.frames("science")
    nsc = int(nsc or cfg.reduction.get("fwdmod", {}).get("nsc", 21))
    dec = float(cfg.target.get("dec_deg", 0.0))
    specs = np.full((len(planets), len(nums), res.nlam), np.nan)
    for i, n in enumerate(nums):
        c = read_cube(dirs.reg / frame_filename(n, suffix, prefix))
        center = res.star_center or image_center(c.shape)
        rot = rotation_from_header(c.h1) - inst.north_offset_deg
        out = c.data.copy() if into_data else np.full(c.data.shape, 1e-15)
        star = star_spectrum_from_header(c.h1, c.nlam)
        for j, pl in enumerate(planets):
            if pl.spectrum is not None:
                spec = np.asarray(pl.spectrum, dtype=float)
            elif empspectrum is not None:
                spec = np.asarray(empspectrum, dtype=float)
            else:
                spec = planet_spectrum_mjy(
                    planet_model_path(cfg, model),
                    star,
                    res.wavelengths_um,
                    res.band,
                    pl.log_contrast,
                )
            specs[j, i] = spec
            px, py = frame_position(pl.x, pl.y, rot, center)
            out = insert_planet(
                out,
                px,
                py,
                spec,
                psf,
                center=center,
                negative=negative,
                h0=c.h0,
                dec_deg=dec,
                lat_deg=inst.latitude_deg,
                nsc=nsc,
            )
        c.h1["FAKEPLAN"] = (len(planets), "number of injected synthetic planets")
        write_cube(dirs.reg / frame_filename(n, suffix + "_fc", prefix), out, c.h0, c.h1)
    return np.nanmedian(specs, axis=1)


def run_reduction(
    cfg: DatasetConfig,
    rc: ReducedCube,
    *,
    reductype: int = 0,
    sdi: ReducedCube | None = None,
    refpath=None,
    fwdmod: bool = True,
    outfile: str = "res0.fits",
    suffix: str | None = None,
) -> Path:
    """Re-run the PSF subtraction recorded in ``rc`` on the ``_fc`` cubes (forward model or fake planets)."""
    from charis_dpp.psfsub.adialoci import adialoci
    from charis_dpp.psfsub.adiklip import adiklip
    from charis_dpp.psfsub.rdiklip import rdiklip
    from charis_dpp.psfsub.sdialoci import sdialoci

    p = dict(rc.params)
    common = dict(rsubval=rc.rsub, fc=True, fwdmod=fwdmod, suffix=suffix)
    stem = Path(outfile).stem
    if rc.algorithm == "aloci":
        if reductype == 0:
            return adialoci(cfg, **p, **common, outfile=outfile)[0]
        if reductype == 1:
            return sdialoci(cfg, **p, **common, outfile=outfile)[0]
        if sdi is None:
            raise ValueError("reductype 2/3 need the SDI reduced cube (sdi_reducname)")
        if reductype == 2:
            adialoci(cfg, **p, **common, norot=True, outfile=f"{stem}_adi_tmp.fits")
            return sdialoci(
                cfg,
                **sdi.params,
                fc=True,
                fwdmod=fwdmod,
                postadi=True,
                suffix="_alocisub",
                outfile=outfile,
            )[0]
        if reductype == 3:
            sdialoci(
                cfg,
                **sdi.params,
                fc=True,
                fwdmod=fwdmod,
                norot=True,
                rsubval=rc.rsub,
                outfile=f"{stem}_sdi_tmp.fits",
            )
            return adialoci(
                cfg,
                **p,
                fc=True,
                fwdmod=fwdmod,
                postsdi=True,
                suffix="_sdialocisub",
                outfile=outfile,
            )[0]
    if rc.algorithm == "klip":
        return adiklip(cfg, **p, **common, outfile=outfile)[0]
    if rc.algorithm == "rdiklip":
        if refpath is None:
            raise ValueError("rdiklip forward model needs refpath")
        p.pop("nzone", None)
        return rdiklip(cfg, refpath=refpath, **p, **common, outfile=outfile)[0]
    raise ValueError(f"cannot re-run reduction of algorithm {rc.algorithm!r}")


def collapsed_header(h: Header) -> Header:
    return h.copy()
