"""Cube registration (``reg/charis_register_cube.pro``).

Every slice of every science cube is shifted so the star lands on the frame
centre.  Three methods:

``sats``   (default) centroid the four satellite spots in the reference
           channel, predict their positions in every other channel from the
           wavelength scaling, re-centroid, and take the spot-centre median;
``psf``    centroid the (unsaturated) stellar core in every channel;
``calsat`` find the centre of a ladder cube with spots as in ``sats``, then
           register each science cube to it by halo cross-correlation.

In all cases the per-channel centres are fitted with a 2nd-order polynomial
in channel index (robust, trimming ``fittrim`` channels at the ends) before
the shifts are applied.  Outputs ``reduc/reg/n####reg.fits`` with ``SATS<l>_<j>``,
``PSFC_<l>``, ``SHIFT_<l>``, ``PSFCENTX/Y`` cards and ``reduc.log``.

Fixes vs. IDL (BUGS.md): #35 method names, #36 sky-subtracted file detection,
#37/#46 reference slice and frame size derived from the data, #38 ``smask``,
#39 frames whose centroiding fails are skipped instead of shifted by 9999 px,
#40 NaN-channel fallback uses the right wavelength ratio, #41 y guess index,
#42 no ``stop`` in method 2, #43 ``keepprev`` uses the previous frame, #51
``fitgauss`` implemented, #59 ``reduc.log`` keyed by frame number.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from charis_dpp.astro.time_coords import ten
from charis_dpp.config.dataset import DatasetConfig, Resolved
from charis_dpp.geometry.apertures import sumaper_im
from charis_dpp.geometry.grids import image_center, radial_grid
from charis_dpp.image.centroid import (
    center_of_symmetry,
    cntrd,
    fft_xcorr_shift,
    fit_gaussian2d,
    gcntrd,
    xcorr_shift,
)
from charis_dpp.image.filters import idl_filter_image, idl_smooth, radial_profile_2d
from charis_dpp.image.interp import shift_sub
from charis_dpp.image.stats import robust_poly_fit
from charis_dpp.io.filelists import frame_filename
from charis_dpp.io.fits import Cube, read_cube, set_sat_positions, write_cube
from charis_dpp.io.reduclog import ReducLog
from charis_dpp.prelims.setupdir import Dirs

log = logging.getLogger(__name__)

__all__ = ["register_cube", "RegistrationOptions", "fit_centers_vs_channel", "centroid_spot"]


@dataclass
class RegistrationOptions:
    method: str = "sats"
    guessoffsets: tuple[float, float] = (0.0, 0.0)
    refslice: int | None = None
    medbox: int = 15
    rsub: bool = False
    fitgauss: bool = False
    fittrim: int = 4
    smoothval: int = 3
    nosmooth: bool = False
    fwhmlim: float = 7.0
    xcorr: bool = False
    smallsteps: bool = False
    revise: bool = False
    keepprev: bool = False
    splitpsf: bool = False
    refcube: int = 0
    rin_corr_arcsec: float = 0.2
    rout_corr_arcsec: float = 0.5
    astrogrid_lambda_d: float | None = None
    ladder: bool = False
    smask: list[tuple[float, float]] | None = None
    prefix: str = "n"
    suffix: str | None = None
    out_suffix: str = "reg"
    reduclog_path: Path | None = None  # default reduc.log (reduc_ladder.log for ladder frames)
    # method='xcorr': how the reference frame's absolute star centre is found:
    # 'symmetry' (180-deg symmetry of the halo), 'centroid' (unsaturated star core),
    # 'given' (use `center` as is)
    xcorr_center: str = "symmetry"
    center: tuple[float, float] | None = None  # target centre; default frame centre
    extra: dict = field(default_factory=dict)


# ----------------------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------------------
def centroid_spot(
    im: np.ndarray,
    x: float,
    y: float,
    fwhmlim: float,
    *,
    fitgauss: bool = False,
    fwhm: float | None = None,
    revise: bool = False,
) -> tuple[float, float]:
    """gcntrd (or 2-D Gaussian fit) with the IDL fallbacks: cntrd/keepcenter, then cntrd/extended box."""
    if fitgauss:
        f = fit_gaussian2d(im, x, y, fwhm or fwhmlim)
        cx, cy = f["x"], f["y"]
    else:
        cx, cy = gcntrd(im, x, y, fwhmlim)
    if not (np.isfinite(cx) and np.isfinite(cy)):
        cx, cy = cntrd(im, x, y, fwhmlim, keepcenter=True)
    if not (np.isfinite(cx) and np.isfinite(cy)) and revise:
        cx, cy = cntrd(im, x, y, fwhmlim * 1.5)
    return cx, cy


def fit_centers_vs_channel(centers: np.ndarray, refslice: int, fittrim: int) -> np.ndarray | None:
    """Robust quadratic fit of (x, y) centre vs channel; returns fitted centres ``(nλ, 2)`` or ``None``.

    Channels ``ceil(fittrim/2)`` at the blue end and ``floor(fittrim/2)`` at the
    red end are excluded from the fit, as are channels whose centroid failed.
    """
    nlam = centers.shape[0]
    lo = int(np.ceil(fittrim / 2.0))
    hi = nlam - int(np.floor(fittrim / 2.0))
    sl = np.arange(nlam, dtype=float)
    use = np.zeros(nlam, dtype=bool)
    use[lo:hi] = True
    use &= np.isfinite(centers).all(axis=1)
    if use.sum() < 3 or not np.isfinite(centers[refslice]).all():
        # fall back: fit through whatever is finite, relative to its median
        if use.sum() < 3:
            return None
        ref = np.nanmedian(centers[use], axis=0)
    else:
        ref = centers[refslice]
    out = np.zeros((nlam, 2))
    for k in range(2):
        c = robust_poly_fit(sl[use] - refslice, centers[use, k] - ref[k], 2, numit=15)
        out[:, k] = ref[k] + np.polyval(c[::-1], sl - refslice)
    return out


def _spots_polar(cens0: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    c0 = np.nanmean(cens0, axis=0)
    d = cens0 - c0
    return c0, np.hypot(d[:, 0], d[:, 1]), np.arctan2(d[:, 1], d[:, 0])


def _predict_spots(cens0: np.ndarray, scale: float) -> np.ndarray:
    """Spot positions in a channel whose wavelength is ``scale`` × the reference one."""
    c0, r, th = _spots_polar(cens0)
    return np.stack([c0[0] + r * scale * np.cos(th), c0[1] + r * scale * np.sin(th)], axis=1)


def _prep_slice(im: np.ndarray, opts: RegistrationOptions, rmax: float, smooth: bool) -> np.ndarray:
    a = im - idl_filter_image(im, median=opts.medbox)
    if opts.rsub:
        a = a - radial_profile_2d(a, 2.0, 1.0, rmax)
    if smooth and not opts.nosmooth and opts.smoothval > 1:
        a = idl_smooth(a, opts.smoothval)
    if opts.smask:
        r_all = None
        for sx, sy in opts.smask:
            cx, cy = gcntrd(a, sx, sy, opts.fwhmlim)
            if not np.isfinite(cx):
                cx, cy = sx, sy
            m = radial_grid(a.shape, (cx, cy)) <= 2
            r_all = m if r_all is None else (r_all | m)
        if r_all is not None:
            a = a.copy()
            a[r_all] = np.nan
    return a


def _wrap_pa(pa: float, dec_deg: float | None, lat_deg: float) -> float:
    """IDL heuristic: for targets south of the zenith wrap PA into (-180, 180]."""
    if dec_deg is not None and dec_deg < lat_deg and pa > 180:
        return pa - 360.0
    return pa


# ----------------------------------------------------------------------------------------
# per-frame measurements
# ----------------------------------------------------------------------------------------
def measure_spot_centers(
    cube: np.ndarray,
    guesses: np.ndarray,
    res: Resolved,
    opts: RegistrationOptions,
    refslice: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Spot centroids ``(nλ, 4, 2)`` and star centres ``(nλ, 2)`` (NaN where a spot failed)."""
    nlam, ny, nx = cube.shape
    wav = res.wavelengths_um
    rmax = min(nx, ny) // 2
    ref_im = _prep_slice(cube[refslice], opts, rmax, smooth=False)
    cens0 = np.array(
        [
            centroid_spot(
                ref_im,
                gx,
                gy,
                opts.fwhmlim,
                fitgauss=opts.fitgauss,
                fwhm=res.fwhm_px[refslice],
                revise=True,
            )
            for gx, gy in guesses
        ]
    )
    if not np.isfinite(cens0).all():
        bad = np.flatnonzero(~np.isfinite(cens0).all(axis=1))
        raise RuntimeError(
            f"satellite spot(s) {bad.tolist()} not found in reference channel {refslice} near {guesses[bad].round(1).tolist()}"
        )
    cens = np.full((nlam, 4, 2), np.nan)
    psfcens = np.full((nlam, 2), np.nan)
    init = np.median(cens0, axis=0)
    for il in range(nlam):
        pred = _predict_spots(cens0, wav[il] / wav[refslice])
        a = _prep_slice(cube[il], opts, rmax, smooth=True)
        if not opts.xcorr:
            c = np.array(
                [
                    centroid_spot(
                        a,
                        px,
                        py,
                        opts.fwhmlim,
                        fitgauss=opts.fitgauss,
                        fwhm=res.fwhm_px[il],
                        revise=opts.revise,
                    )
                    for px, py in pred
                ]
            )
            if np.isfinite(c).all():
                cens[il] = c
                psfcens[il] = np.median(c, axis=0)
        else:
            roi = np.zeros((ny, nx), dtype=bool)
            for px, py in pred:
                roi |= radial_grid((ny, nx), (px, py)) <= 2.0
            xc, yc = center_of_symmetry(a, init[0], init[1], roi, smallsteps=opts.smallsteps)
            psfcens[il] = (xc, yc)
            cens[il] = pred - np.median(pred, axis=0) + np.array([xc, yc])
    return cens, psfcens


def measure_psf_centers(
    cube: np.ndarray,
    guess: tuple[float, float],
    res: Resolved,
    opts: RegistrationOptions,
    prev: np.ndarray | None,
) -> np.ndarray:
    nlam = cube.shape[0]
    out = np.full((nlam, 2), np.nan)
    gx, gy = guess
    for il in range(nlam):
        im = cube[il]
        cx, cy = centroid_spot(
            im,
            gx,
            gy,
            opts.fwhmlim,
            fitgauss=opts.fitgauss,
            fwhm=res.fwhm_px[il],
            revise=opts.revise,
        )
        if opts.revise and opts.splitpsf:
            az = sumaper_im(idl_smooth(im, 3), opts.fwhmlim / 2.0, 80)
            cx, cy = cntrd(az, gx, gy, 4 * opts.fwhmlim)
        if opts.revise and opts.keepprev and prev is not None and np.isfinite(prev[il]).all():
            cx, cy = prev[il]
        if np.isfinite(cx) and np.isfinite(cy):
            out[il] = (cx, cy)
    return out


# ----------------------------------------------------------------------------------------
# apply
# ----------------------------------------------------------------------------------------
def _apply_registration(
    c: Cube,
    fitted: np.ndarray,
    target: tuple[float, float],
    *,
    cens: np.ndarray | None,
    res: Resolved,
    refslice: int,
) -> tuple[np.ndarray, float, float]:
    """Shift every slice so ``fitted[il]`` lands on ``target``; update headers; return (data, xc, yc)."""
    nlam = c.data.shape[0]
    shifts = np.array(target)[None, :] - fitted  # (dx, dy) applied per slice
    out = np.empty_like(c.data)
    for il in range(nlam):
        out[il] = shift_sub(c.data[il], shifts[il, 0], shifts[il, 1])
    psfc = fitted + shifts  # == target for every slice
    if cens is not None:
        cens = cens + shifts[:, None, :]
        # channels whose spots failed: rescale the geometry of the first finite channel
        good = np.isfinite(cens).all(axis=(1, 2))
        if good.any() and not good.all():
            g0 = int(np.flatnonzero(good)[0])
            base = cens[g0]
            for il in np.flatnonzero(~good):
                cens[il] = _predict_spots(base, res.wavelengths_um[il] / res.wavelengths_um[g0])
        set_sat_positions(c.h1, cens)
        for il in range(nlam):
            c.h1[f"PSFC_{il}"] = (
                f"{np.median(cens[il, :, 0]):.3f} {np.median(cens[il, :, 1]):.3f}",
                f"PSF Center of slice {il}",
            )
    else:
        for il in range(nlam):
            c.h1[f"PSFC_{il}"] = (
                f"{psfc[il, 0]:.3f} {psfc[il, 1]:.3f}",
                f"PSF Center of slice {il}",
            )
    for il in range(nlam):
        c.h1[f"SHIFT_{il}"] = (
            f"{shifts[il, 0]:.3f} {shifts[il, 1]:.3f}",
            f"Image shift for slice {il}",
        )
    xc, yc = float(np.median(psfc[:, 0])), float(np.median(psfc[:, 1]))
    c.h1["PSFCENTX"] = (xc, "Mean PSF center X")
    c.h1["PSFCENTY"] = (yc, "Mean PSF center Y")
    c.h1["REFSLICE"] = (refslice, "reference channel for registration")
    return out, xc, yc


# ----------------------------------------------------------------------------------------
# driver
# ----------------------------------------------------------------------------------------
def register_cube(cfg: DatasetConfig | str | Path, **kwargs) -> ReducLog:
    """Register the science (or ladder) cubes; see module docstring.  Keyword arguments
    override :class:`RegistrationOptions` fields and the ``[reduction.register]`` config."""
    cfg = cfg if isinstance(cfg, DatasetConfig) else DatasetConfig.load(cfg)
    rc = cfg.reduction.get("register", {})
    opts = RegistrationOptions(
        method=str(kwargs.pop("method", rc.get("method", "sats"))),
        medbox=int(kwargs.pop("medbox", rc.get("medbox", 15))),
        fittrim=max(int(kwargs.pop("fittrim", rc.get("fittrim", 4))), 0),
        smoothval=int(kwargs.pop("smoothval", rc.get("smoothval", 3))),
        fwhmlim=float(kwargs.pop("fwhmlim", rc.get("fwhmlim", 7.0))),
        rin_corr_arcsec=float(kwargs.pop("rin_corr", rc.get("rin_corr_arcsec", 0.2))),
        rout_corr_arcsec=float(kwargs.pop("rout_corr", rc.get("rout_corr_arcsec", 0.5))),
        astrogrid_lambda_d=kwargs.pop("astrogrid", cfg.observation.get("astrogrid_lambda_d")),
        xcorr_center=str(kwargs.pop("xcorr_center", rc.get("xcorr_center", "symmetry"))),
    )
    for k, v in kwargs.items():
        if hasattr(opts, k):
            setattr(opts, k, v)
        else:
            raise TypeError(f"register_cube: unknown option {k!r}")
    opts.method = {"0": "sats", "1": "psf", "2": "calsat", "calsats": "calsat", "3": "xcorr"}.get(
        str(opts.method), opts.method
    )
    if opts.method not in ("sats", "psf", "calsat", "xcorr"):
        raise ValueError(f"method must be sats/psf/calsat/xcorr, got {opts.method!r}")

    res = cfg.resolve()
    inst = res.instrument
    dirs = Dirs(cfg.workdir)
    dirs.reg.mkdir(parents=True, exist_ok=True)
    nums = cfg.frames("ladder" if opts.ladder else "science")
    if not nums:
        raise ValueError("no frames to register (observation.frames)")
    if opts.suffix is None:
        opts.suffix = (
            "e_skysub"
            if (dirs.prep / frame_filename(nums[0], "e_skysub", opts.prefix)).exists()
            else "e"
        )
        if (
            opts.suffix == "e"
            and (dirs.prep / frame_filename(nums[0], "e_flat_skysub", opts.prefix)).exists()
        ):
            opts.suffix = "e_flat_skysub"
    refslice = res.reference_channel if opts.refslice is None else int(opts.refslice)
    first = read_cube(dirs.prep / frame_filename(nums[0], opts.suffix, opts.prefix))
    ny, nx = first.shape
    xc0, yc0 = image_center((ny, nx)) if opts.center is None else opts.center
    target = (float(xc0), float(yc0))
    dec = cfg.target.get("dec_deg")
    rsat = float(cfg.reduction.get("rsat_px", -1))
    log.info(
        "register_cube: method=%s suffix=%s refslice=%d frames=%d",
        opts.method,
        opts.suffix,
        refslice,
        len(nums),
    )

    # satellite-spot guesses (absolute, for this frame size)
    guesses = (
        inst.spot_guesses_px(res.band, opts.astrogrid_lambda_d)
        + np.array([xc0, yc0])
        + np.array(opts.guessoffsets)
    )
    log_path = (
        Path(opts.reduclog_path)
        if opts.reduclog_path is not None
        else (dirs.reduclog if not opts.ladder else dirs.root / "reduc_ladder.log")
    )
    reduclog = ReducLog(path=log_path)
    if not opts.ladder and log_path.exists():
        try:
            reduclog = ReducLog.load(log_path)
        except Exception:
            reduclog = ReducLog(path=log_path)

    ref_state: dict = {}
    if opts.method == "calsat":
        ref_state = _calsat_reference(cfg, res, dirs, opts, guesses, refslice, target)
    elif opts.method == "xcorr":
        ref_state = _xcorr_reference(cfg, res, dirs, opts, nums, refslice, target)

    prev: np.ndarray | None = None
    for n in nums:
        src = dirs.prep / frame_filename(n, opts.suffix, opts.prefix)
        c = read_cube(src)
        pa0 = _wrap_pa(float(c.get("PA", c.get("PARANG", 0.0))), dec, inst.latitude_deg)
        ha0 = ten(c.get("HA", 0.0))
        try:
            if opts.method == "sats":
                cens, psfcens = measure_spot_centers(c.data, guesses, res, opts, refslice)
            elif opts.method == "psf":
                psfcens = measure_psf_centers(
                    c.data,
                    (xc0 + opts.guessoffsets[0], yc0 + opts.guessoffsets[1]),
                    res,
                    opts,
                    prev,
                )
                cens = None
            elif opts.method == "xcorr":
                psfcens = _xcorr_measure(c.data, ref_state)
                cens = None
            else:
                psfcens = _calsat_measure(c.data, ref_state, res, opts, refslice)
                cens = None
            fitted = fit_centers_vs_channel(psfcens, refslice, opts.fittrim)
            if fitted is None:
                raise RuntimeError("too few channels with a valid centroid")
        except RuntimeError as exc:
            log.error("frame %d: registration failed (%s); skipped", n, exc)
            continue
        data, xc, yc = _apply_registration(c, fitted, target, cens=cens, res=res, refslice=refslice)
        prev = psfcens
        c.h1["REGMETH"] = (opts.method, "registration method")
        for k, v, cm in inst.to_header_cards():
            c.h1[k] = (v, cm)
        c.h0.add_history(f"charis register: method={opts.method} refslice={refslice}")
        out = dirs.reg / frame_filename(n, opts.out_suffix, opts.prefix)
        write_cube(out, data, c.h0, c.h1)
        meas = np.nanmedian(psfcens, axis=0)
        reduclog.add(n, float(meas[0]), float(meas[1]), rsat, ha0, pa0)
        log.info(
            "frame %4d: centre (%.3f, %.3f) -> (%.1f, %.1f)  PA %.2f",
            n,
            meas[0],
            meas[1],
            xc,
            yc,
            pa0,
        )
    reduclog.save()
    cfg.set("geometry.star_center", [target[0], target[1]])
    cfg.set("geometry.reference_channel", refslice)
    if cfg.path is not None:
        cfg.save()
    return reduclog


# ----------------------------------------------------------------------------------------
# calsat (method 2)
# ----------------------------------------------------------------------------------------
def _calsat_reference(cfg, res, dirs, opts, guesses, refslice, target) -> dict:
    lad = cfg.frames("ladder")
    if not lad:
        raise ValueError("method='calsat' needs ladder frames (observation.frames.ladder)")
    ref = read_cube(dirs.prep / frame_filename(lad[0], opts.suffix, opts.prefix))
    o = RegistrationOptions(**{**opts.__dict__, "smask": None})
    o.nosmooth = True
    cens, psfcens = measure_spot_centers(ref.data, guesses, res, o, refslice)
    fitted = fit_centers_vs_channel(psfcens, refslice, opts.fittrim)
    if fitted is None:
        raise RuntimeError("calsat: could not centroid the ladder cube")
    data, xc, yc = _apply_registration(ref, fitted, target, cens=cens, res=res, refslice=refslice)
    write_cube(dirs.reg / "refimage.fits", data, ref.h0, ref.h1)
    pixscale = res.instrument.pixscale_arcsec
    r = radial_grid(data.shape[-2:], (xc, yc))
    roi = (r > opts.rin_corr_arcsec / pixscale) & (r <= opts.rout_corr_arcsec / pixscale)
    return {"ref": data, "roi": roi, "center": (xc, yc)}


def _calsat_measure(cube, ref_state, res, opts, refslice) -> np.ndarray:
    ref, roi = ref_state["ref"], ref_state["roi"]
    nlam = cube.shape[0]
    out = np.full((nlam, 2), np.nan)
    for il in range(nlam):
        tgt = np.nan_to_num(cube[il])
        out[il] = xcorr_shift(np.nan_to_num(ref[il]), tgt, roi & np.isfinite(cube[il]), 0.0, 0.0)
    return out


# ----------------------------------------------------------------------------------------
# xcorr (method 3): no satellite spots needed
# ----------------------------------------------------------------------------------------
def _xcorr_reference(cfg, res, dirs, opts, nums, refslice, target) -> dict:
    """Reference cube for cross-correlation registration: the first science frame (or
    ``opts.refcube``-th), centred per channel on the star.

    The star centre of the reference is found from the 180° symmetry of the
    PSF halo inside ``[rin_corr, rout_corr]`` arcsec (``charis_getcenter``),
    starting from ``opts.center`` (or the frame centre).  Every other frame is
    then aligned to this reference by FFT phase cross-correlation, so the
    method works for data without satellite spots (``[reduction.register]``
    ``method = "xcorr"``); the absolute centre is only as good as the reference's.
    """
    ref_num = nums[int(opts.refcube)] if 0 <= int(opts.refcube) < len(nums) else nums[0]
    ref = read_cube(dirs.prep / frame_filename(ref_num, opts.suffix, opts.prefix))
    ny, nx = ref.shape
    pixscale = res.instrument.pixscale_arcsec
    x0, y0 = target if opts.center is not None else image_center((ny, nx))
    cx0, cy0 = opts.center if opts.center is not None else (x0, y0)
    centres = np.full((ref.nlam, 2), np.nan)
    data = np.empty_like(ref.data)
    if opts.xcorr_center not in ("symmetry", "centroid", "given"):
        raise ValueError("xcorr_center must be 'symmetry', 'centroid' or 'given'")
    for il in range(ref.nlam):
        im = np.where(np.isfinite(ref.data[il]), ref.data[il], 0.0)
        if opts.xcorr_center == "given":
            cx, cy = cx0, cy0
        elif opts.xcorr_center == "centroid":
            cx, cy = centroid_spot(im, cx0, cy0, max(opts.fwhmlim, 2 * res.fwhm_px[il]))
        else:
            r = radial_grid((ny, nx), (cx0, cy0))
            roi = (r > opts.rin_corr_arcsec / pixscale) & (r <= opts.rout_corr_arcsec / pixscale)
            cx, cy = center_of_symmetry(im, cx0, cy0, roi, smallsteps=opts.smallsteps)
        if not (np.isfinite(cx) and np.isfinite(cy)):
            cx, cy = cx0, cy0
        if np.hypot(cx - cx0, cy - cy0) > 5 * res.fwhm_px[il]:
            log.warning(
                "xcorr reference channel %d: centre (%.1f, %.1f) is far from the guess (%.1f, %.1f); "
                "check reduction.register.xcorr_center / center",
                il,
                cx,
                cy,
                cx0,
                cy0,
            )
        centres[il] = (cx, cy)
    fitted = fit_centers_vs_channel(centres, refslice, opts.fittrim)
    if fitted is None:
        raise RuntimeError("xcorr: could not centre the reference cube")
    for il in range(ref.nlam):
        data[il] = shift_sub(ref.data[il], target[0] - fitted[il, 0], target[1] - fitted[il, 1])
    write_cube(dirs.reg / "refimage.fits", data, ref.h0, ref.h1)
    # correlate everything inside 3*rout_corr: masking the core would plant a fixed
    # zeroed disk in both images and bias the cross-correlation towards zero shift
    r = radial_grid((ny, nx), target)
    roi = r <= 3.0 * opts.rout_corr_arcsec / pixscale
    return {"ref": data, "roi": roi, "center": target}


def _xcorr_measure(cube: np.ndarray, ref_state: dict) -> np.ndarray:
    """Per-channel star position in ``cube`` from the shift that aligns it to the reference."""
    ref, roi, (xc, yc) = ref_state["ref"], ref_state["roi"], ref_state["center"]
    nlam = cube.shape[0]
    out = np.full((nlam, 2), np.nan)
    for il in range(nlam):
        dx, dy = fft_xcorr_shift(ref[il], cube[il], roi & np.isfinite(cube[il]))
        out[il] = (xc - dx, yc - dy)
    return out
