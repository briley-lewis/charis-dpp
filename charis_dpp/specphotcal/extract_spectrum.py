"""Extract a 1-D spectrum of a point source (or a disk surface brightness) from a
reduced cube (``charis_extract_1d_spectrum[_disk].pro``).

Per channel: optional radial-profile subtraction, optional local plane
background fit, optional high-pass; aperture photometry with ``R_AP<l>``
(sky annulus ``[2, bklim]·R_AP`` or zero sky); SNR from the local robust
noise; error ``sqrt((F/SNR)² + (F·CERR)²)``; optional throughput correction.

Fixes vs. IDL: the dead calibration call and GUI prompts are gone (#100/#103),
the throughput file is an explicit argument, ``fitrange`` works (#101),
``breakout`` catches failed centroids (#102), frame geometry comes from the
data (#104), disk mode uses aperture photometry per arcsec² consistently
(#107/#108), the plot label uses the real units and the east sign (#105).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from charis_dpp.config.dataset import DatasetConfig
from charis_dpp.geometry.apertures import aper_idl, aperture_weights
from charis_dpp.geometry.grids import image_center, radial_grid
from charis_dpp.image.centroid import cntrd, gcntrd
from charis_dpp.image.filters import idl_filter_image, radial_profile_2d
from charis_dpp.image.snr import snr_map
from charis_dpp.io.fits import read_cube, write_image
from charis_dpp.prelims.setupdir import Dirs

log = logging.getLogger(__name__)

__all__ = ["ExtractedSpectrum", "extract_spectrum", "plane_background"]


@dataclass
class ExtractedSpectrum:
    wavelength_um: np.ndarray
    flux: np.ndarray
    err: np.ndarray
    snr: np.ndarray
    position: tuple[float, float]
    snr_collapsed: float
    units: str = "mJy"
    throughput: np.ndarray | None = None
    path: Path | None = None
    extra: dict = field(default_factory=dict)

    def offsets_arcsec(self, center: tuple[float, float], pixscale: float) -> tuple[float, float]:
        """(ΔRA east-positive, ΔDec north-positive) in arcsec for a north-up, east-left image."""
        return (
            -(self.position[0] - center[0]) * pixscale,
            (self.position[1] - center[1]) * pixscale,
        )


def plane_background(
    im: np.ndarray, x: float, y: float, r_src: float, r_out: float, bklir: float, center
) -> np.ndarray:
    """Subtract, inside the source region, a plane fitted to the surrounding annular background."""
    r_s = radial_grid(im.shape, (x, y))
    r_c = radial_grid(im.shape, center)
    sep = np.hypot(x - center[0], y - center[1])
    ring = (np.abs(r_c - sep) < bklir * r_src / 2.0) if bklir > 0 else np.ones_like(im, bool)
    bkg = (r_s <= r_out) & (r_s > r_src) & ring & np.isfinite(im)
    src = r_s <= r_src
    if bkg.sum() < 10:
        return im
    yy, xx = np.nonzero(bkg)
    A = np.column_stack([np.ones(xx.size), xx, yy])
    coeff, *_ = np.linalg.lstsq(A, im[bkg], rcond=None)
    out = im.copy()
    ys, xs = np.nonzero(src)
    out[src] = im[src] - (coeff[0] + coeff[1] * xs + coeff[2] * ys)
    return out


def _centroid(
    im: np.ndarray, x: float, y: float, fwhm: float, order=("cntrd", "gcntrd")
) -> tuple[float, float]:
    for name in order:
        fn = cntrd if name == "cntrd" else gcntrd
        cx, cy = fn(im, x, y, fwhm)
        if np.isfinite(cx) and np.isfinite(cy):
            return cx, cy
    return np.nan, np.nan


def extract_spectrum(
    cfg: DatasetConfig | str | Path,
    datacube: str | Path = "final.fits",
    coords: tuple[float, float] = (0.0, 0.0),
    *,
    fitcollapse: bool = True,
    fitslice: bool = False,
    fitbckgd: bool = False,
    extended: bool = False,
    throughput: str | Path | np.ndarray | None = None,
    fitrange: list[int] | None = None,
    nozerosky: bool = False,
    noprad: bool = False,
    filt_slice: bool = False,
    filtsnr: bool = False,
    nocalerror: bool = False,
    bklim: float = 4.0,
    bklir: float = 2.0,
    mode: str = "point",
    outputspec: str | None = None,
    plot: bool = True,
    r_ap: float | None = None,
) -> ExtractedSpectrum:
    cfg = cfg if isinstance(cfg, DatasetConfig) else DatasetConfig.load(cfg)
    res = cfg.resolve()
    dirs = Dirs(cfg.workdir)
    path = (
        Path(datacube)
        if Path(datacube).is_absolute() or Path(datacube).exists()
        else dirs.proc / datacube
    )
    c = read_cube(path)
    data = c.data.copy()
    nlam = c.nlam
    wav = res.wavelengths_um
    fwhm = res.fwhm_px
    center = res.star_center or image_center(c.shape)
    units = str(c.h1.get("FLUXUNIT", "mJy")).strip()
    aperrad = (
        np.array([float(c.h1.get(f"R_AP{l}", 0.5 * fwhm[l])) for l in range(nlam)])
        if r_ap is None
        else np.full(nlam, float(r_ap))
    )
    cerr = np.array([float(c.h1.get(f"CERR{l}", 0.0)) for l in range(nlam)])
    fitrange = list(range(nlam)) if fitrange is None else list(fitrange)
    data_col = np.nanmedian(data[fitrange], axis=0)
    fwhm_col = float(np.median(fwhm))
    apmed = float(np.median(aperrad[fitrange]))
    xposf, yposf = float(coords[0]), float(coords[1])
    if fitcollapse:
        image = data_col
        if fitbckgd:
            image = plane_background(
                data_col, xposf, yposf, 2 * apmed, bklim * apmed, bklir, center
            )
            xposf, yposf = _centroid(image, xposf, yposf, fwhm_col, ("gcntrd", "cntrd"))
        else:
            xposf, yposf = _centroid(image, xposf, yposf, fwhm_col, ("cntrd", "gcntrd"))
        if not np.isfinite(xposf):
            raise RuntimeError(f"centroid failed near {coords}")
        if extended:
            m = radial_grid(image.shape, (xposf, yposf)) <= 2 * apmed
            yy, xx = np.nonzero(m)
            w = image[m]
            xposf, yposf = float(np.sum(xx * w) / np.sum(w)), float(np.sum(yy * w) / np.sum(w))
    flux = np.full(nlam, np.nan)
    err = np.full(nlam, np.nan)
    snr = np.full(nlam, np.nan)
    pixscale = res.instrument.pixscale_arcsec
    for il in range(nlam):
        im = data[il]
        x, y = xposf, yposf
        if fitslice:
            x, y = _centroid(im, xposf, yposf, fwhm[il], ("gcntrd", "cntrd"))
            if not np.isfinite(x):
                x, y = xposf, yposf
        if not noprad:
            im = im - radial_profile_2d(im, 1.0, 0.0, None, center)
        if fitbckgd:
            im = plane_background(im, x, y, 2 * aperrad[il], bklim * aperrad[il], bklir, center)
        if filt_slice:
            im = im - idl_filter_image(im, median=int(round(10 * fwhm[il])))
        ap = aperrad[il]
        if nozerosky:
            f, *_ = aper_idl(im, x, y, ap, skyrad=(2 * ap, bklim * ap), setskyval=None)
        else:
            f, *_ = aper_idl(im, x, y, ap, setskyval=0.0)
        if mode == "disk":
            area = float(np.sum(aperture_weights(im.shape, x, y, ap))) * pixscale**2
            f = f / area  # per arcsec^2
        r = snr_map(
            im,
            fwhm[il],
            max(c.shape) // 2,
            coord=(x, y),
            fixpos=True,
            zero=True,
            finite=True,
            highpass=filtsnr,
            highpass_factor=3.0,
        )
        s = r.snr if np.isfinite(r.snr) and r.snr != 0 else np.nan
        flux[il] = f
        snr[il] = s
        calterm = 0.0 if nocalerror else f * cerr[il]
        err[il] = np.sqrt((f / s) ** 2 + calterm**2) if np.isfinite(s) else np.nan
    rcol = snr_map(
        data_col,
        fwhm_col,
        max(c.shape) // 2,
        coord=(xposf, yposf),
        fixpos=True,
        zero=True,
        finite=True,
        highpass=filtsnr,
        highpass_factor=3.0,
    )
    write_image(dirs.proc / "snrmapcol.fits", rcol.snrmap)
    thr = None
    if throughput is not None:
        if isinstance(throughput, np.ndarray):
            thr = np.asarray(throughput, dtype=float)
        else:
            p = Path(throughput) if Path(throughput).exists() else dirs.proc / str(throughput)
            thr = np.loadtxt(p, ndmin=2)[:, 1]
        if thr.size != nlam:
            raise ValueError(f"throughput has {thr.size} entries, cube has {nlam} channels")
        flux = flux / thr
        err = err / thr
    name = outputspec or ("spectrum_disk.dat" if mode == "disk" else "spectrum.dat")
    out = dirs.proc / name if not Path(name).is_absolute() else Path(name)
    np.savetxt(
        out,
        np.column_stack([wav, flux, err, snr]),
        fmt="%.6f %.6e %.6e %.3f",
        header=f"wavelength_um flux[{units}{'/arcsec2' if mode == 'disk' else ''}] err snr",
    )
    result = ExtractedSpectrum(
        wav, flux, err, snr, (xposf, yposf), float(rcol.snr), units, thr, out
    )
    if plot:
        from charis_dpp.plotting import plot_spectrum

        good = np.isfinite(flux)
        plot_spectrum(
            wav[good],
            flux[good],
            err[good],
            units=units if mode == "point" else f"{units}/arcsec^2",
            title=f"{path.name} at ({xposf:.1f}, {yposf:.1f})  SNR_col={rcol.snr:.1f}",
            savefig=out.with_suffix(".png"),
        )
    dx, dy = result.offsets_arcsec(center, pixscale)
    log.info(
        "spectrum at (%.2f, %.2f) px = (%+.3f E, %+.3f N) arcsec; collapsed SNR %.1f -> %s",
        xposf,
        yposf,
        dx,
        dy,
        rcol.snr,
        out,
    )
    return result
