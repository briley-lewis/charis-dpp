"""Spectrophotometric calibration (``charis_specphot_cal.pro`` + ``_unsat.pro``).

Each registered cube is multiplied, channel by channel, by
``FSCALE = F_star(λ) · atten(λ) / F_spot(λ)`` where ``F_spot`` is the measured
satellite-spot flux (or the unsaturated stellar flux with ``unsat=True``,
``atten = 1``), ``atten`` the spot/star contrast and ``F_star`` the model
stellar flux density from :mod:`photcal_calc`.

``calmethod``: 0 each cube calibrated with its own spots; 1 the median cube;
2 user-supplied calibration cube(s); 3 as 0 with the per-channel attenuation
from ``ladder/atten_fac_ladder``; 4 the ladder cubes; 5 copy ``FSCALE``/
``FSTAR``/``CERR`` from another star's calibrated cube.

Header cards written (consumed downstream): ``FLUXUNIT``, ``FSTAR_<l>``,
``R_AP<l>``, ``SKY_IN``/``SKY_OUT``, ``FSCALE<l>``, ``CERR<l>``.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from charis_dpp.astro.units import canonical_units
from charis_dpp.cals.sat_contrast import modulation_from_header, sat_contrast
from charis_dpp.config.dataset import DatasetConfig, Resolved
from charis_dpp.geometry.apertures import aper_idl
from charis_dpp.geometry.grids import radial_grid
from charis_dpp.image.filters import idl_filter_image, radial_profile_2d
from charis_dpp.image.snr import snr_map
from charis_dpp.image.stats import resistant_mean, robust_sigma
from charis_dpp.io.filelists import frame_filename
from charis_dpp.io.fits import Cube, get_psf_center, get_sat_positions, read_cube, write_cube
from charis_dpp.prelims.setupdir import Dirs
from charis_dpp.specphotcal.photcal_calc import stellar_flux_density

log = logging.getLogger(__name__)

__all__ = ["specphot_cal", "measure_spot_fluxes", "SpotPhotometry", "read_fscale"]

SPOT_SELECT = {
    None: None,
    "all": None,
    "left": (0, 3),
    "right": (1, 2),
    "upperleft": (0,),
    "lowerleft": (3,),
    "upperright": (1,),
    "lowerright": (2,),
}


@dataclass
class SpotPhotometry:
    flux: np.ndarray  # (nλ,) combined spot flux
    err: np.ndarray  # (nλ,)
    per_spot: np.ndarray  # (nλ, 4)
    noise: np.ndarray  # (nλ, 4)
    aperture_px: np.ndarray  # (nλ,)


def _prep_cal_slice(
    im: np.ndarray, fwhm: float, *, nopradcal: bool, filtcal: bool, rmax: float, rmin: float = 1.0
) -> np.ndarray:
    out = im.copy()
    if not nopradcal:
        im2 = np.where((im == 0) | ~np.isfinite(im), np.nan, im)
        pr = radial_profile_2d(im2, 1.0, rmin, rmax)
        good = np.isfinite(im) & (im != 0)
        out[good] = im[good] - pr[good]
    if filtcal:
        out = out - idl_filter_image(out, median=int(round(10 * fwhm)))
    return out


def _plane_fit_background(
    im: np.ndarray, x: float, y: float, r_src: float, r_bkg: float, center
) -> np.ndarray:
    """Subtract a plane fitted to the annular background around a spot (``fitbckgd``)."""
    r_spot = radial_grid(im.shape, (x, y))
    r_star = radial_grid(im.shape, center)
    src_dist = np.hypot(x - center[0], y - center[1])
    region = (r_spot <= r_bkg) & (np.abs(r_star - src_dist) < r_src)
    bkg = region & (r_spot > r_src) & np.isfinite(im)
    if bkg.sum() < 10:
        return im
    yy, xx = np.nonzero(bkg)
    A = np.column_stack([np.ones(xx.size), xx, yy])
    coeff, *_ = np.linalg.lstsq(A, im[bkg], rcond=None)
    out = im.copy()
    yy2, xx2 = np.nonzero(region)
    out[region] = im[region] - (coeff[0] + coeff[1] * xx2 + coeff[2] * yy2)
    return out


def measure_spot_fluxes(
    cube: np.ndarray,
    positions: np.ndarray,
    aperture_px: np.ndarray,
    *,
    subskyannulus: bool = False,
    sky_annulus: tuple[float, float] = (2.0, 6.0),
    fitbckgd: bool = False,
    meancomb: bool = False,
    expandederror: bool = False,
    spot_select: str | None = None,
    center: tuple[float, float] | None = None,
    sys_frac: tuple[float, float] = (0.01, 0.02),
) -> SpotPhotometry:
    """Aperture photometry of the spots (``positions`` ``(nλ, nspot, 2)``) in every channel.

    The noise at each spot is the local robust azimuthal scatter of the
    aperture-summed image with the spots masked (small-sample corrected), as
    in ``charis_snratio_sub``; the spot-2 index bug of the IDL code (#67) is
    fixed.  The error adds a 1 % (median) / 2 % (mean) systematic term, or the
    robust scatter of the four spots with ``expandederror``.
    """
    nlam, ny, nx = cube.shape
    nspot = positions.shape[1]
    per = np.full((nlam, nspot), np.nan)
    noise = np.full((nlam, nspot), np.nan)
    flux = np.full(nlam, np.nan)
    err = np.full(nlam, np.nan)
    sel = SPOT_SELECT.get(spot_select)
    if spot_select is not None and spot_select not in SPOT_SELECT:
        raise ValueError(f"spot_select must be one of {list(SPOT_SELECT)}")
    xc, yc = center if center is not None else (nx // 2, ny // 2)
    for il in range(nlam):
        im = cube[il]
        ap = float(aperture_px[il])
        pos = positions[il]
        if fitbckgd:
            for q in range(nspot):
                if np.isfinite(pos[q]).all():
                    im = _plane_fit_background(im, pos[q, 0], pos[q, 1], 5 * ap, 15 * ap, (xc, yc))
        for q in range(nspot):
            if not np.isfinite(pos[q]).all():
                continue
            if subskyannulus:
                per[il, q] = aper_idl(
                    im,
                    pos[q, 0],
                    pos[q, 1],
                    ap,
                    skyrad=(sky_annulus[0] * ap, sky_annulus[1] * ap),
                    setskyval=None,
                )[0]
            else:
                per[il, q] = aper_idl(im, pos[q, 0], pos[q, 1], ap, setskyval=0.0)[0]
        masked = im.copy()
        for q in range(nspot):
            if np.isfinite(pos[q]).all():
                masked[radial_grid(im.shape, (pos[q, 0], pos[q, 1])) <= ap] = np.nan
        nm = snr_map(masked, 2 * ap, max(nx, ny) // 2, center=(xc, yc), finite=True).noisemap
        for q in range(nspot):
            if np.isfinite(pos[q]).all():
                iy, ix = int(round(pos[q, 1])), int(round(pos[q, 0]))
                if 0 <= iy < ny and 0 <= ix < nx:
                    noise[il, q] = nm[iy, ix]
        vals = per[il]
        nz = noise[il]
        if sel is None:
            good = np.isfinite(vals)
            if meancomb:
                flux[il] = resistant_mean(vals[good], 3.0) if good.sum() > 1 else np.nanmean(vals)
                stat_in = resistant_mean(nz[good], 3.0) if good.sum() > 1 else np.nanmean(nz)
            else:
                flux[il] = np.nanmedian(vals)
                stat_in = np.nanmedian(nz)
        else:
            flux[il] = np.nanmean(vals[list(sel)])
            stat_in = (
                np.nanstd(vals[list(sel)], ddof=1) if len(sel) > 1 else np.nanmean(nz[list(sel)])
            )
        if expandederror:
            sys = robust_sigma(vals[np.isfinite(vals)])
        else:
            sys = (sys_frac[1] if meancomb else sys_frac[0]) * flux[il]
        err[il] = np.sqrt(stat_in**2 + sys**2)
    return SpotPhotometry(flux, err, per, noise, np.asarray(aperture_px, dtype=float))


def read_fscale(h1) -> tuple[np.ndarray, np.ndarray, np.ndarray, str]:
    """``(fscale, cerr, fstar, units)`` from a calibrated cube's extension header."""
    nlam = int(h1.get("NAXIS3", 0))
    if nlam == 0:
        nlam = len([k for k in h1 if k.startswith("FSCALE")])
    fscale = np.array([float(h1[f"FSCALE{l}"]) for l in range(nlam)])
    cerr = np.array([float(h1.get(f"CERR{l}", np.nan)) for l in range(nlam)])
    fstar = np.array([float(h1.get(f"FSTAR_{l}", np.nan)) for l in range(nlam)])
    return fscale, cerr, fstar, str(h1.get("FLUXUNIT", "mJy")).strip()


def _spot_positions_for(c: Cube, nlam: int, unsat: bool) -> np.ndarray:
    if unsat:
        xc, yc = get_psf_center(c.h1) or (c.shape[1] // 2, c.shape[0] // 2)
        return np.tile(np.array([[xc, yc]], dtype=float), (nlam, 1, 1))
    sats = get_sat_positions(c.h1, nlam)
    if not np.isfinite(sats).any():
        raise ValueError(
            f"{c.path}: no SATS<l>_<j> cards — run registration with method='sats' first"
        )
    return sats


def specphot_cal(
    cfg: DatasetConfig | str | Path,
    *,
    calmethod: int = 0,
    unsat: bool = False,
    calcube: list[Path] | None = None,
    datacube: list[Path] | None = None,
    modamp: float | None = None,
    atten_fac: float | np.ndarray | None = None,
    subskyannulus: bool = False,
    spot_select: str | None = None,
    fitbckgd: bool = False,
    expandederror: bool = False,
    filtcal_slice: bool = False,
    nopradcal: bool = False,
    meancomb: bool = False,
    fixradius: float | None = None,
    ap_factor: float | None = None,
    filtername: str | None = None,
    starlib: str | None = None,
    empspectrum: str | Path | None = None,
    av: float | None = None,
    fluxunits: str = "mJy",
    prefix: str = "n",
    suffix: str = "reg",
    out_suffix: str | None = None,
) -> list[Path]:
    cfg = cfg if isinstance(cfg, DatasetConfig) else DatasetConfig.load(cfg)
    res: Resolved = cfg.resolve()
    inst = res.instrument
    dirs = Dirs(cfg.workdir)
    sp = cfg.reduction.get("specphot", {})
    ap_factor = float(ap_factor if ap_factor is not None else sp.get("ap_factor", 1.0))
    starlib = starlib or sp.get("starlib", "pickles")
    av = float(av if av is not None else sp.get("av", 0.0))
    units = canonical_units(fluxunits)
    sky_ann = tuple(
        sp.get("sky_annulus_unsat" if unsat else "sky_annulus", (4, 10) if unsat else (2, 6))
    )
    nlam = res.nlam
    wav = res.wavelengths_um
    aperrad = (
        0.5 * res.fwhm_px * ap_factor if fixradius is None else np.full(nlam, float(fixradius))
    )

    # --- which cubes calibrate, which cubes get calibrated ----------------------------------
    sci = [dirs.reg / frame_filename(n, suffix, prefix) for n in cfg.frames("science")]
    if datacube:
        targets = [Path(p) for p in datacube]
    else:
        targets = sci
    if calmethod in (0, 1, 3):
        cal_paths = sci
    elif calmethod == 4:
        cal_paths = [dirs.reg / frame_filename(n, suffix, prefix) for n in cfg.frames("ladder")]
    elif calmethod in (2, 5):
        if not calcube:
            raise ValueError(f"calmethod={calmethod} needs calcube=[...]")
        cal_paths = [Path(p) for p in calcube]
    else:
        raise ValueError("calmethod must be 0..5")
    if not cal_paths or not targets:
        raise ValueError("no calibration or target cubes (check observation.frames)")
    missing = [p for p in cal_paths + targets if not p.exists()]
    if missing:
        raise FileNotFoundError(f"missing cubes: {missing[:3]}{'...' if len(missing) > 3 else ''}")
    first = read_cube(cal_paths[0])
    mjd = float(first.get("MJD", inst.mjd))

    # --- calmethod 5: copy the calibration of another star ----------------------------------
    if calmethod == 5:
        rows = [read_fscale(read_cube(p).h1) for p in cal_paths]
        fscale = np.nanmedian([r[0] for r in rows], axis=0)
        cerr = np.nanmedian([r[1] for r in rows], axis=0)
        fstar = np.nanmedian([r[2] for r in rows], axis=0)
        units = rows[0][3]
        return _write_calibrated(
            targets,
            np.tile(fscale, (len(targets), 1)),
            cerr,
            fstar,
            aperrad,
            sky_ann,
            units,
            dirs,
            out_suffix or f"{suffix}_cal",
            datacube is not None,
        )

    # --- attenuation ---------------------------------------------------------------------------
    if unsat:
        atten = np.ones(nlam)
    elif calmethod == 3 or atten_fac is not None:
        if atten_fac is None:
            f = dirs.ladder / "atten_fac_ladder"
            if not f.exists():
                raise FileNotFoundError(f"{f} not found: run satspot_cal first")
            atten_fac = _read_atten_file(f, nlam)
        atten, _ = sat_contrast(
            25.0, wav, inst=inst, manual=-2.5 * np.log10(np.asarray(atten_fac, dtype=float))
        )
    else:
        mod = modamp if modamp is not None else cfg.observation.get("modulation_nm")
        if mod is None:
            mod = modulation_from_header(first.h0, first.h1)
        atten, _ = sat_contrast(float(mod), wav, mjd, inst=inst)

    # --- load and pre-process calibration cubes ----------------------------------------------
    rmax = min(first.shape) // 2
    cal_cubes = []
    cal_pos = []
    for p in cal_paths:
        c = read_cube(p)
        proc = np.stack(
            [
                _prep_cal_slice(
                    c.data[il],
                    res.fwhm_px[il],
                    nopradcal=nopradcal,
                    filtcal=filtcal_slice,
                    rmax=rmax,
                )
                for il in range(nlam)
            ]
        )
        cal_cubes.append(proc)
        cal_pos.append(_spot_positions_for(c, nlam, unsat))
    if calmethod in (1, 2, 4) and len(cal_cubes) > 1:
        stack = np.stack(cal_cubes)
        cal_cubes = [np.nanmean(stack, axis=0) if meancomb else np.nanmedian(stack, axis=0)]
        cal_pos = [np.nanmedian(np.stack(cal_pos), axis=0)]
    elif calmethod in (1, 2, 4):
        cal_cubes, cal_pos = cal_cubes[:1], cal_pos[:1]

    # --- spot photometry -------------------------------------------------------------------
    phots = [
        measure_spot_fluxes(
            cube,
            pos,
            aperrad,
            subskyannulus=subskyannulus,
            sky_annulus=sky_ann,
            fitbckgd=fitbckgd,
            meancomb=meancomb,
            expandederror=expandederror,
            spot_select=spot_select,
            center=res.star_center,
        )
        for cube, pos in zip(cal_cubes, cal_pos, strict=True)
    ]
    if len(phots) > 1 and calmethod not in (0, 3):
        satflux = np.nanmedian([p.flux for p in phots], axis=0)
        esat = np.nanmedian([p.err for p in phots], axis=0)
        per_target = np.tile(satflux, (len(targets), 1))
    elif calmethod in (0, 3):
        if len(phots) != len(targets):
            raise ValueError(
                "calmethod 0/3 calibrates each cube with its own spots: calibration and target lists must match"
            )
        per_target = np.stack([p.flux for p in phots])
        satflux = np.nanmedian(per_target, axis=0)
        esat = np.nanmedian([p.err for p in phots], axis=0)
    else:
        satflux, esat = phots[0].flux, phots[0].err
        per_target = np.tile(satflux, (len(targets), 1))
    snr_sats = satflux / esat

    # --- stellar flux --------------------------------------------------------------------------
    band = res.band
    filt_band = filtername or band
    mags = cfg.target.get("mags", {})
    magkey = {"J": "J", "H": "H", "K": "K", "lowres": "H"}.get(filt_band, "H")
    if magkey not in mags:
        raise ValueError(f"target.mags lacks the {magkey} magnitude needed for calibration")
    spt = str(cfg.target.get("sptype", "")).strip()
    if len(spt) == 2:
        spt += "V"
    dloglam = first.get("DLOGLAM")
    fstar = stellar_flux_density(
        spt,
        float(mags[magkey]),
        band,
        wav,
        library=starlib,
        av=av,
        empspectrum=empspectrum,
        units=units,
        dloglam=float(dloglam) if dloglam else None,
    )
    fscale = fstar[None, :] * atten[None, :] / per_target
    cerr = 1.0 / snr_sats
    log.info(
        "specphot_cal: calmethod=%d unsat=%s units=%s  median FSCALE %.4g",
        calmethod,
        unsat,
        units,
        np.nanmedian(fscale),
    )
    return _write_calibrated(
        targets,
        fscale,
        cerr,
        fstar,
        aperrad,
        sky_ann,
        units,
        dirs,
        out_suffix or f"{suffix}_cal",
        datacube is not None,
        extra={"ATTEN": atten, "SATFLUX": satflux},
    )


def _read_atten_file(path: Path, nlam: int) -> np.ndarray:
    vals = [float(x) for line in path.read_text().splitlines() for x in line.split() if _isnum(x)]
    if not vals:
        raise ValueError(f"{path} is empty")
    arr = np.array(vals)
    return arr if arr.size in (1, nlam) else arr[:1]


def _isnum(s: str) -> bool:
    try:
        float(s)
        return True
    except ValueError:
        return False


def _write_calibrated(
    targets,
    fscale,
    cerr,
    fstar,
    aperrad,
    sky_ann,
    units,
    dirs,
    out_suffix,
    explicit_paths,
    extra=None,
) -> list[Path]:
    outs = []
    for i, p in enumerate(targets):
        c = read_cube(p)
        data = c.data * fscale[i][:, None, None]
        h1 = c.h1
        h1["FLUXUNIT"] = (units, "Flux Density Units in Cube")
        for l in range(data.shape[0]):
            h1[f"FSTAR_{l}"] = (float(fstar[l]), " Star Flux")
            h1[f"R_AP{l}"] = (float(aperrad[l]), "Aperture Radius (Pixels)")
            h1[f"FSCALE{l}"] = (float(fscale[i][l]), f"scale to convert counts to {units}")
            h1[f"CERR{l}"] = (float(cerr[l]), f"Cal Fractional error for slice {l}")
            if extra:
                h1[f"ATTEN{l}"] = (float(extra["ATTEN"][l]), "spot/star attenuation")
        h1["SKY_IN"] = (
            float(sky_ann[0]),
            "Inner Radius for Sky Annulus (in Units of Aperture Radius)",
        )
        h1["SKY_OUT"] = (
            float(sky_ann[1]),
            "Outer Radius for Sky Annulus (in Units of Aperture Radius)",
        )
        if explicit_paths:
            out = p.with_name(p.stem + "_cal.fits")
        else:
            from charis_dpp.io.filelists import frame_number

            out = dirs.reg / frame_filename(frame_number(p), out_suffix)
        write_cube(out, data, c.h0, h1)
        outs.append(out)
    return outs
