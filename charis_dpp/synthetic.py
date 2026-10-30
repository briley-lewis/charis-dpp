"""Synthetic CHARIS datasets for tests and demonstrations.

:func:`make_dataset` writes raw-style ``CRSA<seq>_cube.fits`` files, a
dataset configuration and a ``truth.json`` describing everything that was
injected: star centre per frame, parallactic angles, satellite-spot and planet
positions, fluxes and the north-up rotation of each frame.  The optical model
is deliberately simple (Gaussian PSF core + power-law halo, pupil-fixed
speckles that scale with wavelength, Gaussian spots/planet) but every header
keyword the pipeline reads is populated, so the full chain can run on it.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.time import Time

from charis_dpp.cals.sat_contrast import sat_contrast
from charis_dpp.config.dataset import DatasetConfig
from charis_dpp.config.instrument import instrument_for, to_mjd
from charis_dpp.geometry.grids import image_center, radial_grid
from charis_dpp.io.catalogs import StarInfo
from charis_dpp.io.wavelengths import wavelengths_um
from charis_dpp.pdi.stokes import polarisation_angle_grid

__all__ = ["SyntheticConfig", "make_dataset", "psf_image", "DEFAULT_STAR"]

DEFAULT_STAR = StarInfo(
    name="SYNTH 1",
    ra_deg=346.86965,
    dec_deg=21.13425,
    sptype="A5V",
    plx_mas=25.04,
    plx_err_mas=0.85,
    pmra_mas_yr=107.88,
    pmra_err_mas_yr=0.75,
    pmdec_mas_yr=-50.08,
    pmdec_err_mas_yr=0.64,
    mags={"J": 5.383, "H": 5.280, "K": 5.240},
    mag_errors={"J": 0.027, "H": 0.018, "K": 0.018},
)


@dataclass
class Planet:
    sep_px: float = 22.0
    pa_deg: float = 60.0  # position angle in the north-up frame, from +x towards +y
    log_contrast: float = -4.0  # planet/star flux ratio (log10), flat spectrum unless `spectrum`
    spectrum: list[float] | None = None  # per-channel contrast multipliers


@dataclass
class SyntheticConfig:
    nframes: int = 12
    shape: tuple[int, int] = (101, 101)
    band: str = "lowres"
    date: str = "2017-08-14"
    utc_start: str = "08:00:00"
    frame_spacing_s: float = 40.0
    pa_start_deg: float = -30.0
    pa_step_deg: float = 5.0
    reads: int = 20  # LASTRD - FIRSTRD + 1
    star_flux: float = 2.0e5  # peak counts/s at the reference wavelength
    fwhm_factor: float = 1.03  # PSF FWHM in units of lambda/D
    halo_amp: float = 0.02
    halo_index: float = -2.0
    n_speckles: int = 30
    speckle_contrast: float = 5e-4
    speckle_r_range: tuple[float, float] = (6.0, 38.0)
    sat_spots: bool = True
    astrogrid_lambda_d: float = 15.9
    modulation_nm: float = 25.0
    planet: Planet | None = field(default_factory=Planet)
    dither_px: float = 0.8
    sky_level: float = 5.0
    noise: float = 2.0
    nsky: int = 3
    seed: int = 0
    object_name: str = "SYNTH 1"
    # polarimetry: two Wollaston beams at ±polcent offset, an HWP cycling 0/22.5/45/67.5 deg
    # per frame, a tangentially polarised ring and a VAMPIRES-style HWP log
    pdi: bool = False
    pol_ring_r_px: float = 18.0
    pol_ring_width_px: float = 2.5
    pol_ring_contrast: float = 2e-3  # ring peak / star peak
    pol_fraction: float = 0.5


def psf_image(
    shape,
    x: float,
    y: float,
    fwhm: float,
    peak: float,
    halo_amp: float = 0.0,
    halo_index: float = -2.0,
) -> np.ndarray:
    """Gaussian core of the given FWHM plus a power-law halo ``halo_amp·peak·(r/fwhm)^index`` beyond 1 FWHM."""
    ny, nx = shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    r2 = (xx - x) ** 2 + (yy - y) ** 2
    sig = fwhm / 2.35482
    im = peak * np.exp(-0.5 * r2 / sig**2)
    if halo_amp > 0:
        r = np.sqrt(r2)
        halo = halo_amp * peak * np.where(r > fwhm, (np.maximum(r, 1e-3) / fwhm) ** halo_index, 1.0)
        im = im + halo
    return im


def _header_common(
    cfg: SyntheticConfig,
    inst,
    seq: int,
    mjd: float,
    parang: float,
    ha_h: float,
    npa: float,
    star: StarInfo,
):
    from astropy import units as u
    from astropy.coordinates import Angle

    t = Time(mjd, format="mjd")
    h = fits.Header()
    h["OBJECT"] = cfg.object_name
    h["ORIGNAME"] = (f"{seq:08d}", "raw sequence number")
    h["FILTNAME"] = "Broadband" if cfg.band == "lowres" else cfg.band
    h["CAL_BAND"] = cfg.band
    exptime = cfg.reads * inst.read_time_s
    t_end = Time(mjd + exptime / 86400.0, format="mjd")
    h["UTC-DATE"] = t.iso[:10]
    h["UTC-TIME"] = (t_end.iso[11:23], "end of the exposure (as in CHARIS cubes)")
    h["MJD"] = mjd
    h["EXPTIME"] = exptime
    h["FIRSTRD"] = 1
    h["LASTRD"] = cfg.reads
    h["RA"] = Angle(star.ra_deg * u.deg).to_string(unit=u.hour, sep=":", precision=2, pad=True)
    h["DEC"] = Angle(star.dec_deg * u.deg).to_string(
        unit=u.deg, sep=":", precision=1, pad=True, alwayssign=True
    )
    h["PARANG"] = (parang, "parallactic angle [deg]")
    h["HA"] = Angle(ha_h * u.hour).to_string(
        unit=u.hour, sep=":", precision=1, pad=True, alwayssign=True
    )
    h["ALTITUDE"] = 60.0
    h["AZIMUTH"] = 180.0 + parang
    h["X_GRDAMP"] = (cfg.modulation_nm * 1e-3 if cfg.sat_spots else 0.0, "sat spot modulation [um]")
    h["Y_GRDAMP"] = h["X_GRDAMP"]
    h["X_GRDSEP"] = (
        cfg.astrogrid_lambda_d if cfg.sat_spots else 0.0,
        "sat spot separation [lambda/D]",
    )
    s = inst.pixscale_arcsec / 3600.0
    a = np.deg2rad(npa)
    h["XPIXSCAL"] = s
    h["YPIXSCAL"] = s
    h["CD1_1"] = np.cos(a) * s
    h["CD1_2"] = np.sin(a) * s
    h["CD2_1"] = -np.sin(a) * s
    h["CD2_2"] = np.cos(a) * s
    h["CRPIX1"] = cfg.shape[1] // 2 + 1
    h["CRPIX2"] = cfg.shape[0] // 2 + 1
    h["CRVAL1"] = star.ra_deg
    h["CRVAL2"] = star.dec_deg
    h["CTYPE1"] = "RA---TAN"
    h["CTYPE2"] = "DEC--TAN"
    h["TELESCOP"] = "Subaru"
    h["INSTRUME"] = "CHARIS"
    h["DLOGLAM"] = 0.0324
    return h


def make_dataset(
    outdir: str | Path, cfg: SyntheticConfig | None = None, star: StarInfo = DEFAULT_STAR
) -> Path:
    """Write the synthetic dataset to ``outdir``; returns the dataset config path."""
    cfg = cfg or SyntheticConfig()
    outdir = Path(outdir)
    raw = outdir / "data" / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(cfg.seed)
    inst = instrument_for(cfg.date)
    wav = wavelengths_um(cfg.band)
    nlam = wav.size
    ny, nx = cfg.shape
    xc0, yc0 = image_center(cfg.shape)
    lod = np.asarray(inst.lambda_over_d_px(wav))
    fwhm = cfg.fwhm_factor * lod
    ref = int(np.argmin(np.abs(wav - np.sqrt(wav[0] * wav[-1]))))
    mjd0 = to_mjd(f"{cfg.date}T{cfg.utc_start}")
    north_offset = inst.north_offset_deg
    atten, _ = sat_contrast(cfg.modulation_nm, wav, inst=inst)
    star_spec = cfg.star_flux * np.ones(nlam)  # flat peak spectrum

    # pupil-fixed speckles, positions at the reference wavelength
    sp_r = rng.uniform(*cfg.speckle_r_range, cfg.n_speckles)
    sp_t = rng.uniform(0, 2 * np.pi, cfg.n_speckles)
    sp_a = rng.uniform(0.3, 1.0, cfg.n_speckles) * cfg.speckle_contrast

    # satellite spots: instrument-fixed, point-symmetric about the star (real spots are
    # diffraction copies); the first angle follows the default guess so the guesses work
    g = inst.spot_guesses_px(cfg.band)
    spot_ang = np.arctan2(g[0, 1], g[0, 0]) + np.deg2rad([0.0, -90.0, 180.0, 90.0])

    truth: dict = {
        "config": asdict(cfg),
        "band": cfg.band,
        "nlam": nlam,
        "wavelengths_um": wav.tolist(),
        "fwhm_px": fwhm.tolist(),
        "reference_channel": ref,
        "north_offset_deg": north_offset,
        "epoch": inst.epoch_label,
        "star_peak": star_spec.tolist(),
        "sat_atten": atten.tolist(),
        "frames": [],
        "sky_frames": [],
    }
    planet = cfg.planet
    if planet is not None:
        pspec = np.array(planet.spectrum) if planet.spectrum is not None else np.ones(nlam)
        planet_flux = star_spec * 10**planet.log_contrast * pspec
        px_sky = xc0 + planet.sep_px * np.cos(np.deg2rad(planet.pa_deg))
        py_sky = yc0 + planet.sep_px * np.sin(np.deg2rad(planet.pa_deg))
        truth["planet"] = {
            "x_northup": float(px_sky),
            "y_northup": float(py_sky),
            "sep_px": planet.sep_px,
            "pa_deg": planet.pa_deg,
            "log_contrast": planet.log_contrast,
            "contrast": (10**planet.log_contrast * pspec).tolist(),
            "peak": planet_flux.tolist(),
        }

    hwp_angles = [float(a) for a in inst.pdi.hwp_angles_deg]
    polcent = np.array(inst.pdi.polcent_offset_px, dtype=float) * (
        nx / float(inst.detector.cube_size_px)
    )
    hwp_log_rows: list[tuple[float, float, float]] = []
    if cfg.pdi:
        truth["pdi"] = {
            "polcent_offset_px": polcent.tolist(),
            "ring_r_px": cfg.pol_ring_r_px,
            "ring_width_px": cfg.pol_ring_width_px,
            "ring_contrast": cfg.pol_ring_contrast,
            "pol_fraction": cfg.pol_fraction,
            "hwp_angles_deg": hwp_angles,
        }
    seq = 1
    science_nums, sky_nums = [], []
    for n in range(cfg.nframes):
        parang = cfg.pa_start_deg + n * cfg.pa_step_deg
        rot_n = parang + north_offset  # rotation the pipeline must apply for north-up
        npa = (
            rot_n + north_offset
        )  # header CD such that rotation_from_header - north_offset == rot_n
        mjd = mjd0 + n * cfg.frame_spacing_s / 86400.0
        ha_h = -1.0 + n * cfg.frame_spacing_s / 3600.0
        dx, dy = rng.uniform(-cfg.dither_px, cfg.dither_px, 2)
        xc, yc = xc0 + dx, yc0 + dy
        cube = np.zeros((nlam, ny, nx))
        spots_xy = np.zeros((nlam, 4, 2))
        spots_left = np.zeros((nlam, 4, 2))
        planet_xy = None
        hwp_angle = None
        if cfg.pdi:
            hwp_angle = float(hwp_angles[n % len(hwp_angles)])
            hwp_log_rows.append((hwp_angle, mjd - 5.0 / 86400.0, mjd - 4.0 / 86400.0))
        if planet is not None:
            t = np.deg2rad(rot_n)
            ct, st = np.cos(t), np.sin(t)
            # inverse of rotat's forward map: p_det = c + R^T (p_sky - c)
            ux, uy = px_sky - xc0, py_sky - yc0
            pdx = ct * ux - st * uy
            pdy = st * ux + ct * uy
            planet_xy = (xc + pdx, yc + pdy)
        for l in range(nlam):

            def unpolarised(xb, yb, l=l, spots=None, planet_xy=planet_xy, xc=xc, yc=yc):
                im = psf_image(
                    cfg.shape, xb, yb, fwhm[l], star_spec[l], cfg.halo_amp, cfg.halo_index
                )
                scale = wav[l] / wav[ref]
                for k in range(cfg.n_speckles):
                    sx = xb + sp_r[k] * scale * np.cos(sp_t[k])
                    sy = yb + sp_r[k] * scale * np.sin(sp_t[k])
                    im += psf_image(cfg.shape, sx, sy, fwhm[l], sp_a[k] * star_spec[l])
                if cfg.sat_spots:
                    rs = cfg.astrogrid_lambda_d * lod[l]
                    for j in range(4):
                        sxj = xb + rs * np.cos(spot_ang[j])
                        syj = yb + rs * np.sin(spot_ang[j])
                        if spots is not None:
                            spots[l, j] = (sxj, syj)
                        im += psf_image(cfg.shape, sxj, syj, fwhm[l], atten[l] * star_spec[l])
                if planet_xy is not None:
                    im += psf_image(
                        cfg.shape,
                        xb + planet_xy[0] - xc,
                        yb + planet_xy[1] - yc,
                        fwhm[l],
                        planet_flux[l],
                    )
                return im

            if not cfg.pdi:
                im = unpolarised(xc, yc, spots=spots_xy)
            else:
                # tangentially polarised ring about each beam centre (rotation invariant)
                def ring(xb, yb, l=l):
                    r = radial_grid(cfg.shape, (xb, yb))
                    phi = polarisation_angle_grid(cfg.shape, (xb, yb))
                    i_ring = (
                        cfg.pol_ring_contrast
                        * star_spec[l]
                        * np.exp(-0.5 * ((r - cfg.pol_ring_r_px) / cfg.pol_ring_width_px) ** 2)
                    )
                    q = -cfg.pol_fraction * i_ring * np.cos(2 * phi)
                    u = -cfg.pol_fraction * i_ring * np.sin(2 * phi)
                    return i_ring, q, u

                th = np.deg2rad(hwp_angle)
                im = np.zeros(cfg.shape)
                for sign, spots in ((-1.0, spots_left), (1.0, spots_xy)):
                    xb, yb = xc + sign * polcent[0], yc + sign * polcent[1]
                    i_ring, q, u = ring(xb, yb)
                    q_meas = q * np.cos(4 * th) + u * np.sin(4 * th)  # ideal HWP + Wollaston
                    im += 0.5 * (unpolarised(xb, yb, spots=spots) + i_ring + sign * q_meas)
            im += cfg.sky_level + rng.normal(0, cfg.noise, im.shape)
            cube[l] = im
        h0 = _header_common(cfg, inst, seq, mjd, parang, ha_h, npa, star)
        _write_raw(raw / f"CRSA{seq:08d}_cube.fits", cube, h0, cfg.noise, hwp_angle)
        truth["frames"].append(
            {
                "seq": seq,
                "filenum": seq,
                "mjd": mjd,
                "parang": parang,
                "ha_hours": ha_h,
                "rot_northup_deg": rot_n,
                "center": [float(xc), float(yc)],
                "spots": spots_xy.tolist(),
                "spots_left": spots_left.tolist() if cfg.pdi else None,
                "planet": None if planet_xy is None else [float(planet_xy[0]), float(planet_xy[1])],
                "hwp_angle": hwp_angle,
            }
        )
        science_nums.append(seq)
        seq += 1
    for n in range(cfg.nsky):
        mjd = mjd0 + (cfg.nframes + n) * cfg.frame_spacing_s / 86400.0
        parang = cfg.pa_start_deg + (cfg.nframes + n) * cfg.pa_step_deg
        cube = cfg.sky_level + rng.normal(0, cfg.noise, (nlam, ny, nx))
        h0 = _header_common(cfg, inst, seq, mjd, parang, 0.0, parang + 2 * north_offset, star)
        _write_raw(raw / f"CRSA{seq:08d}_cube.fits", cube, h0, cfg.noise)
        truth["sky_frames"].append({"seq": seq, "filenum": seq, "mjd": mjd})
        sky_nums.append(seq)
        seq += 1

    if cfg.pdi:
        lines = ["angle,ut_start,ut_end,target"]
        for a, t0, t1 in hwp_log_rows:
            lines.append(
                f"{a},{Time(t0, format='mjd').iso[11:19]},{Time(t1, format='mjd').iso[11:19]},"
                f"{cfg.object_name}"
            )
        (outdir / f"vampHWPLog_{cfg.date.replace('-', '')}.csv").write_text("\n".join(lines) + "\n")
    ds = DatasetConfig.new(
        target=star.to_config(),
        observation={
            "date": cfg.date,
            "mjd": mjd0,
            "band": cfg.band,
            "astrogrid_lambda_d": cfg.astrogrid_lambda_d,
            "modulation_nm": cfg.modulation_nm,
            "frames": {"science": _rng(science_nums), "sky": _rng(sky_nums), "ladder": ""},
        },
        reduction={"rmin_px": 5, "rmax_px": int(min(nx, ny) // 2 - 6), "fwhm_px": float(fwhm[ref])},
    )
    ds.set("geometry.frame_shape", [ny, nx])
    name = f"{cfg.object_name.replace(' ', '')}_{cfg.band}"
    cfg_path = ds.save(outdir / f"{name}.toml")
    (outdir / "truth.json").write_text(json.dumps(truth, indent=1))
    return cfg_path


def _rng(nums: list[int]) -> str:
    from charis_dpp.io.filelists import format_range

    return format_range(nums)


def _write_raw(
    path: Path, cube: np.ndarray, h0: fits.Header, noise: float, hwp_angle: float | None = None
) -> None:
    primary = fits.PrimaryHDU(header=h0)
    h1 = h0.copy()
    data = fits.ImageHDU(data=cube.astype(np.float32), header=h1, name="DATA")
    ivar = fits.ImageHDU(data=np.full(cube.shape, 1.0 / noise**2, dtype=np.float32), name="IVAR")
    h3 = fits.Header()
    h3["D_IMRANG"] = (0.0, "derotator angle [deg]")
    if hwp_angle is not None:
        h3["RET-ANG1"] = (hwp_angle, "HWP angle [deg]")
    h3["MJD-STR"] = h0["MJD"]
    h3["MJD-END"] = h0["MJD"] + h0["EXPTIME"] / 86400.0
    h3["PARANG"] = h0["PARANG"]
    raw3 = fits.ImageHDU(header=h3, name="RAWHDR")
    fits.HDUList([primary, data, ivar, raw3]).writeto(path, overwrite=True)
