"""Automatic end-to-end reduction (``scripts/charis_autoreduce.pro``).

Steps: newobs → imprep → subtract-sky → register → specphot-cal →
make-emp-psf → imrsub → adialoci → detect → fwdmod-planet → extract-spectrum
→ attenmap → final-contrast.  ``start``/``stop`` select a sub-range; the
target name, date and band are read from the first raw cube when not given.
"""

from __future__ import annotations

import logging
from pathlib import Path

from astropy.io import fits

from charis_dpp.cli.registry import register
from charis_dpp.config.dataset import DatasetConfig, canonical_band
from charis_dpp.io.catalogs import StarInfo

log = logging.getLogger(__name__)

STEPS = [
    "newobs",
    "imprep",
    "skysub",
    "register",
    "specphotcal",
    "emppsf",
    "imrsub",
    "adialoci",
    "detect",
    "fwdmod",
    "extract",
    "attenmap",
    "contrast",
]


def _star_from_catalog_or_header(target: str, h) -> StarInfo:
    """SIMBAD/VizieR lookup; without network (or an unknown name) fall back to the header coordinates.

    The dataset config then carries no magnitudes or parallax: fill them in by
    hand before ``specphot-cal`` (the calibration needs the star's magnitude).
    """
    from charis_dpp.io.catalogs import query_star

    try:
        return query_star(target)
    except (LookupError, OSError, RuntimeError) as exc:
        log.warning("catalog lookup for %r failed (%s); using the header coordinates", target, exc)
        from astropy import units as u
        from astropy.coordinates import Angle

        ra = h.get("RA_DEG", h.get("RA", 0.0))
        dec = h.get("DEC_DEG", h.get("DEC", 0.0))
        if isinstance(ra, str):
            ra = Angle(ra, unit=u.hourangle).deg
        if isinstance(dec, str):
            dec = Angle(dec, unit=u.deg).deg
        return StarInfo(target, float(ra), float(dec))


def autoreduce(
    workdir: str | Path = ".",
    *,
    target: str | None = None,
    date: str | None = None,
    band: str | None = None,
    star: StarInfo | None = None,
    config: str | Path | None = None,
    start: str = "newobs",
    stop: str = "contrast",
    modamp: float | None = None,
    astrogrid: float | None = None,
    starlib: str = "pickles",
    rsub: bool = True,
    ntc: int = 5,
    psf: str = "empirical",
) -> Path:
    workdir = Path(workdir)
    steps = STEPS[STEPS.index(start) : STEPS.index(stop) + 1]
    cfg_path = Path(config) if config else None
    if "newobs" in steps:
        raw = sorted(workdir.glob("CRSA*cube.fits")) or sorted(
            (workdir / "data" / "raw").glob("CRSA*cube.fits")
        )
        if not raw:
            raise FileNotFoundError("no CRSA*cube.fits in the working directory or data/raw")
        h = fits.getheader(raw[len(raw) // 2])
        target = target or str(h.get("OBJECT", "target")).strip()
        band = (
            band or canonical_band(str(h.get("CAL_BAND", h.get("FILTNAME", "lowres")))) or "lowres"
        )
        date = date or str(h.get("UTC-DATE", "")).replace("-", "")
        from charis_dpp.prelims.newobs import newobs

        if star is None:
            star = _star_from_catalog_or_header(target, h)
        if astrogrid is None and float(h.get("X_GRDSEP", 0.0) or 0.0) > 0:
            astrogrid = float(h["X_GRDSEP"])
        cfg_path = newobs(
            date,
            target,
            band,
            workdir=workdir,
            star=star,
            modulation_nm=modamp,
            astrogrid_lambda_d=astrogrid,
        )
    if cfg_path is None:
        cands = sorted(workdir.glob("*.toml"))
        if len(cands) != 1:
            raise ValueError("give config= (several or no dataset .toml files found)")
        cfg_path = cands[0]
    cfg = DatasetConfig.load(cfg_path)
    if "imprep" in steps:
        from charis_dpp.prelims.imprep import imprep

        imprep(cfg)
    if "skysub" in steps and cfg.frames("sky"):
        from charis_dpp.skysub import subtract_sky

        subtract_sky(cfg)
    if "register" in steps:
        from charis_dpp.reg import register_cube

        register_cube(cfg)
    if "specphotcal" in steps:
        from charis_dpp.specphotcal import specphot_cal

        specphot_cal(cfg, calmethod=0, modamp=modamp, starlib=starlib)
    if "emppsf" in steps:
        from charis_dpp.specphotcal import make_empirical_psf

        make_empirical_psf(cfg)
    if "imrsub" in steps and rsub:
        from charis_dpp.filt import imrsub

        imrsub(cfg, prad=True)
    if "adialoci" in steps:
        from charis_dpp.psfsub import adialoci

        adialoci(cfg, rsubval=int(rsub), outfile="aloci.fits")
    dets = []
    if "detect" in steps:
        from charis_dpp.tools.detection import find_detections

        dets = find_detections(
            cfg, "aloci.fits", outfile=workdir / "reduc" / "proc" / "aloci_det.targ"
        )
    if "fwdmod" in steps and dets:
        from charis_dpp.partialsub.fwdmod_planet import fwdmod_planet
        from charis_dpp.partialsub.inject import Planet

        fwdmod_planet(
            cfg,
            reducname="aloci.fits",
            planets=[Planet(d.x, d.y, d.log_contrast) for d in dets],
            psf=psf,
        )
    if "extract" in steps and dets:
        from charis_dpp.specphotcal.extract_spectrum import extract_spectrum

        for i, d in enumerate(dets):
            extract_spectrum(
                cfg,
                "aloci.fits",
                (d.x, d.y),
                throughput=f"synth_throughput{i}.txt",
                fitcollapse=True,
                filtsnr=True,
                outputspec=f"aloci_auto_spec_{i}.dat",
            )
    if "attenmap" in steps:
        from charis_dpp.partialsub.attenmap import attenuation_map

        attenuation_map(cfg, reducname="aloci.fits", ntc=ntc, psf=psf)
    if "contrast" in steps:
        from charis_dpp.tools.contrast import final_contrast

        final_contrast(cfg, "aloci.fits")
    return cfg_path


@register("autoreduce", help="Run the standard reduction chain end to end")
def _add(parser):
    parser.add_argument("--target")
    parser.add_argument("--date")
    parser.add_argument("--band")
    parser.add_argument("--config", type=Path, help="existing dataset .toml (skips newobs)")
    parser.add_argument("--start", choices=STEPS, default="newobs")
    parser.add_argument("--stop", choices=STEPS, default="contrast")
    parser.add_argument("--modamp", type=float, help="sat-spot modulation [nm] (default: header)")
    parser.add_argument(
        "--astrogrid", type=float, help="sat-spot separation [lambda/D] (default: header X_GRDSEP)"
    )
    parser.add_argument("--starlib", default="pickles")
    parser.add_argument("--no-rsub", action="store_true")
    parser.add_argument("--ntc", type=int, default=5)
    parser.add_argument("--psf", default="empirical")

    def run(args) -> int:
        p = autoreduce(
            args.workdir,
            target=args.target,
            date=args.date,
            band=args.band,
            config=args.config,
            start=args.start,
            stop=args.stop,
            modamp=args.modamp,
            astrogrid=args.astrogrid,
            starlib=args.starlib,
            rsub=not args.no_rsub,
            ntc=args.ntc,
            psf=args.psf,
        )
        print(f"autoreduce finished; config {p}")
        return 0

    return run
