"""5σ contrast curves (``charis_calc_final_contrast.pro``, ``charis_calc_raw_contrast.pro``).

Fixes vs. IDL: the ``maskcoord`` keyword works (#230), ``finite``/``zero``/``filt``
are honoured (#231), sub-band channel sets come from the wavelength table
(#232), the typos that crashed ``meanadd`` and ``databb`` are gone (#228/#229);
raw contrast masks each channel with its own aperture (#233) and uses the
correct spot index (#234).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from charis_dpp.cals.sat_contrast import modulation_from_header, sat_contrast
from charis_dpp.config.dataset import DatasetConfig
from charis_dpp.geometry.apertures import aper_idl
from charis_dpp.geometry.grids import image_center, radial_grid
from charis_dpp.image.filters import radial_profile, radial_profile_2d
from charis_dpp.image.snr import snr_map
from charis_dpp.image.stats import combine_frames
from charis_dpp.io.filelists import frame_filename
from charis_dpp.io.fits import get_sat_positions, read_cube, read_image, write_cube, write_image
from charis_dpp.planetspec.planet_photcal import star_spectrum_from_header
from charis_dpp.prelims.setupdir import Dirs

log = logging.getLogger(__name__)

__all__ = ["ContrastCurves", "final_contrast", "raw_contrast"]


@dataclass
class ContrastCurves:
    sep_arcsec: np.ndarray
    broadband: np.ndarray
    bands: dict[str, np.ndarray]
    per_channel: np.ndarray  # (nλ, nsep)
    outdir: Path


def _noise(im: np.ndarray, fwhm: float, center, *, zero=True, finite=True, filt=True) -> np.ndarray:
    return snr_map(
        im, fwhm, int(min(im.shape) / 2), center=center, zero=zero, finite=finite, highpass=filt
    ).noisemap


def final_contrast(
    cfg: DatasetConfig | str | Path,
    datacube: str | Path = "final.fits",
    *,
    attenmap: str | Path | None = "attenmap.fits",
    attenmap_lambda: str | Path | None = "attenmap_lambda.fits",
    maskcoord: tuple[float, float] | None = None,
    meanadd: bool = False,
    finite: bool = True,
    zero: bool = True,
    filt: bool = True,
    rmax_px: float | None = None,
    plot: bool = True,
) -> ContrastCurves:
    cfg = cfg if isinstance(cfg, DatasetConfig) else DatasetConfig.load(cfg)
    res = cfg.resolve()
    dirs = Dirs(cfg.workdir)
    path = (
        Path(datacube)
        if Path(datacube).is_absolute() or Path(datacube).exists()
        else dirs.proc / datacube
    )
    c = read_cube(path)
    cube = c.data
    nlam = c.nlam
    fwhm = res.fwhm_px
    pixscale = res.instrument.pixscale_arcsec
    center = res.star_center or image_center(c.shape)
    outdir = dirs.proc / "final_contrast"
    outdir.mkdir(parents=True, exist_ok=True)
    method = "mean" if meanadd else "median"

    def _atten(p):
        if p is None:
            return None
        pp = Path(p) if Path(p).is_absolute() or Path(p).exists() else dirs.proc / str(p)
        return read_image(pp) if pp.exists() else None

    def _positive(a):
        # throughput ≤ 0 (extrapolated inside the innermost measured radius, or a
        # failed recovery) is not a measurement: no contrast there rather than a negative one
        return None if a is None else np.where(np.isfinite(a) & (a > 0), a, np.nan)

    att_l = _positive(_atten(attenmap_lambda))
    att_c = _positive(_atten(attenmap))
    if att_l is None and attenmap_lambda is not None:
        log.warning(
            "no attenuation map found (%s); contrast curves are not throughput-corrected",
            attenmap_lambda,
        )
    starflux = star_spectrum_from_header(c.h1, nlam)
    if not np.isfinite(starflux).any():
        raise ValueError("cube has no FSTAR_<l> cards; run specphot-cal first")
    mask = None
    if maskcoord is not None:
        mask = radial_grid(c.shape, maskcoord) <= 1.5 * float(fwhm.max())
    rmax = float(rmax_px or min(c.shape) // 2 - 1)

    def masked(im):
        return np.where(mask, np.nan, im) if mask is not None else im

    contrast_slice = np.empty_like(cube)
    for il in range(nlam):
        noise = _noise(masked(cube[il]), fwhm[il], center, zero=zero, finite=finite, filt=filt)
        if att_l is not None:
            with np.errstate(divide="ignore", invalid="ignore"):
                noise = noise / att_l[il]
        contrast_slice[il] = 5 * noise / starflux[il]
    bands = {}
    for name, idx in res.sub_band_channels.items():
        if idx.size == 0:
            continue
        im = combine_frames(cube[idx], method, axis=0)
        noise = _noise(
            masked(im), float(np.median(fwhm[idx])), center, zero=zero, finite=finite, filt=filt
        )
        if att_l is not None:
            noise = noise / combine_frames(att_l[idx], method, axis=0)
        bands[name] = 5 * noise / float(np.median(starflux[idx]))
    col = combine_frames(cube, method, axis=0)
    noise = _noise(masked(col), float(np.median(fwhm)), center, zero=zero, finite=finite, filt=filt)
    if att_c is not None:
        noise = noise / att_c
    contrast_bb = 5 * noise / float(np.median(starflux))
    write_cube(outdir / "contrast_slice.fits", contrast_slice, c.h0, c.h1)
    write_image(outdir / "contrast_bb.fits", contrast_bb, c.h1)
    sep, pr_bb = radial_profile(contrast_bb, 1.0, 1.0, rmax, center)
    np.savetxt(
        outdir / "broadband_contrast.txt", np.column_stack([sep * pixscale, pr_bb]), fmt="%.5f %.6e"
    )
    curves = {}
    for name, cmap in bands.items():
        _, pr = radial_profile(cmap, 1.0, 1.0, rmax, center)
        curves[name] = pr
        np.savetxt(
            outdir / f"{name.lower()}_contrast.txt",
            np.column_stack([sep * pixscale, pr]),
            fmt="%.5f %.6e",
        )
    per = np.empty((nlam, sep.size))
    for il in range(nlam):
        _, per[il] = radial_profile(contrast_slice[il], 1.0, 1.0, rmax, center)
        np.savetxt(
            outdir / f"slice_{il}_contrast.txt",
            np.column_stack([sep * pixscale, per[il]]),
            fmt="%.5f %.6e",
        )
    if plot:
        from charis_dpp.plotting import plot_contrast

        plot_contrast(
            sep * pixscale,
            {
                **{f"ch {il}": per[il] for il in range(0, nlam, max(nlam // 6, 1))},
                **curves,
                "broadband": pr_bb,
            },
            title=path.name,
            savefig=outdir / "contrast_slice.png",
        )
    log.info("final contrast written to %s", outdir)
    return ContrastCurves(sep * pixscale, pr_bb, curves, per, outdir)


def raw_contrast(
    cfg: DatasetConfig | str | Path,
    *,
    files: list[Path] | None = None,
    modamp: float | None = None,
    sub: bool = True,
    ap_factor: float | None = None,
    prefix: str = "n",
    suffix: str = "reg_cal",
) -> ContrastCurves:
    """Raw (un-PSF-subtracted) 5σ contrast normalised by the satellite-spot flux / attenuation."""
    cfg = cfg if isinstance(cfg, DatasetConfig) else DatasetConfig.load(cfg)
    res = cfg.resolve()
    dirs = Dirs(cfg.workdir)
    inst = res.instrument
    nums = cfg.frames("science")
    files = files or [dirs.reg / frame_filename(n, suffix, prefix) for n in nums]
    ap_factor = float(
        ap_factor
        if ap_factor is not None
        else cfg.reduction.get("specphot", {}).get("ap_factor", 1.0)
    )
    first = read_cube(files[0])
    nlam = first.nlam
    wav, fwhm = res.wavelengths_um, res.fwhm_px
    aperrad = 0.5 * fwhm * ap_factor
    center = res.star_center or image_center(first.shape)
    pixscale = inst.pixscale_arcsec
    mod = (
        modamp
        if modamp is not None
        else cfg.observation.get("modulation_nm") or modulation_from_header(first.h0, first.h1)
    )
    atten, _ = sat_contrast(float(mod), wav, float(first.get("MJD", inst.mjd)), inst=inst)
    outdir = dirs.proc / ("raw_contrasts" if sub else "raw_contrasts_rms")
    outdir.mkdir(parents=True, exist_ok=True)
    cubes = [read_cube(p) for p in files]
    stack = np.stack([c.data for c in cubes])
    med = np.nanmedian(stack, axis=0)
    sats = np.nanmedian(np.stack([get_sat_positions(c.h1, nlam) for c in cubes]), axis=0)
    rmax = min(first.shape) // 2 - 1
    contrast = np.empty_like(med)
    for il in range(nlam):
        im = med[il].copy()
        im = (
            im
            - radial_profile_2d(
                np.where(np.isfinite(im) & (im != 0), im, np.nan), 1.0, 0.0, rmax, center
            )
            if sub
            else np.abs(im)
        )
        satflux = np.nanmedian(
            [
                aper_idl(med[il], sx, sy, aperrad[il], setskyval=0.0)[0]
                for sx, sy in sats[il]
                if np.isfinite(sx)
            ]
        )
        masked = im.copy()
        for sx, sy in sats[il]:
            if np.isfinite(sx):
                masked[radial_grid(im.shape, (sx, sy)) <= aperrad[il]] = np.nan
        noise = snr_map(masked, 2 * aperrad[il], rmax, center=center, finite=True).noisemap
        contrast[il] = 5 * noise / (satflux / atten[il])
    write_cube(outdir / "raw_contrast_slice.fits", contrast, first.h0, first.h1)
    sep, _ = radial_profile(contrast[0], 1.0, 1.0, rmax, center)
    per = np.empty((nlam, sep.size))
    curves = {}
    for il in range(nlam):
        _, per[il] = radial_profile(contrast[il], 1.0, 1.0, rmax, center)
        np.savetxt(
            outdir / f"slice_{il}_contrast.txt",
            np.column_stack([sep * pixscale, per[il]]),
            fmt="%.5f %.6e",
        )
    for name, idx in res.sub_band_channels.items():
        if idx.size:
            curves[name] = np.nanmedian(per[idx], axis=0)
    bb = np.nanmedian(per, axis=0)
    np.savetxt(
        outdir / "broadband_contrast.txt", np.column_stack([sep * pixscale, bb]), fmt="%.5f %.6e"
    )
    return ContrastCurves(sep * pixscale, bb, curves, per, outdir)
