"""``charis detect`` / ``final-contrast`` / ``raw-contrast`` (+ atmosfit, empcomp, prop-mot when available)."""

from __future__ import annotations

from pathlib import Path

from charis_dpp.cli.registry import register


@register("detect", help="Find point sources above an SNR threshold in a reduced cube")
def _det(parser):
    parser.add_argument("config", type=Path)
    parser.add_argument("--datacube", default="final.fits")
    parser.add_argument("--threshold", type=float)
    parser.add_argument("--excl", type=float)
    parser.add_argument("--rangeconstant", type=float)
    parser.add_argument(
        "--outfile", type=Path, help="target file (x y log contrast) for fwdmod-planet"
    )
    parser.add_argument("--snrmap", type=Path)

    def run(args) -> int:
        from charis_dpp.tools.detection import find_detections

        dets = find_detections(
            args.config,
            args.datacube,
            threshold=args.threshold,
            excl=args.excl,
            rangeconstant=args.rangeconstant,
            outfile=args.outfile,
            snrmap_out=args.snrmap,
        )
        for d in dets:
            print(f"({d.x:7.2f}, {d.y:7.2f})  SNR {d.snr:6.1f}  log contrast {d.log_contrast:6.2f}")
        print(f"{len(dets)} detection(s)")
        return 0

    return run


@register("final-contrast", help="5-sigma contrast curves of a PSF-subtracted cube")
def _fc(parser):
    parser.add_argument("config", type=Path)
    parser.add_argument("--datacube", default="final.fits")
    parser.add_argument("--attenmap", default="attenmap.fits")
    parser.add_argument("--attenmap-lambda", default="attenmap_lambda.fits")
    parser.add_argument("--noattenmap", action="store_true")
    parser.add_argument("--maskcoord", nargs=2, type=float, metavar=("X", "Y"))
    parser.add_argument("--meanadd", action="store_true")
    parser.add_argument("--nofinite", action="store_true")
    parser.add_argument("--nozero", action="store_true")
    parser.add_argument("--nofilt", action="store_true")

    def run(args) -> int:
        from charis_dpp.tools.contrast import final_contrast

        c = final_contrast(
            args.config,
            args.datacube,
            attenmap=None if args.noattenmap else args.attenmap,
            attenmap_lambda=None if args.noattenmap else args.attenmap_lambda,
            maskcoord=tuple(args.maskcoord) if args.maskcoord else None,
            meanadd=args.meanadd,
            finite=not args.nofinite,
            zero=not args.nozero,
            filt=not args.nofilt,
        )
        print(f"wrote contrast curves to {c.outdir}")
        return 0

    return run


@register(
    "raw-contrast", help="Raw (no PSF subtraction) 5-sigma contrast from the registered cubes"
)
def _rc(parser):
    parser.add_argument("config", type=Path)
    parser.add_argument("--modamp", type=float)
    parser.add_argument("--nosub", action="store_true", help="no radial-profile subtraction")
    parser.add_argument("--ap-factor", type=float)
    parser.add_argument("--suffix", default="reg_cal")

    def run(args) -> int:
        from charis_dpp.tools.contrast import raw_contrast

        c = raw_contrast(
            args.config,
            modamp=args.modamp,
            sub=not args.nosub,
            ap_factor=args.ap_factor,
            suffix=args.suffix,
        )
        print(f"wrote raw contrast curves to {c.outdir}")
        return 0

    return run


@register(
    "spec-covar",
    help="Spectral covariance of the residual noise at a separation (Greco & Brandt 2016)",
)
def _cov(parser):
    parser.add_argument("config", type=Path)
    parser.add_argument("--datacube", default="final.fits")
    parser.add_argument("--rho", type=float, help="separation in lambda/D")
    parser.add_argument("--width", type=float)
    parser.add_argument("--step", type=float)
    parser.add_argument("--prad", action="store_true")
    parser.add_argument("--filt", action="store_true")
    parser.add_argument("--normnoise", action="store_true")
    parser.add_argument("--maskcoord", nargs=2, type=float, metavar=("X", "Y"))

    def run(args) -> int:
        from charis_dpp.specphotcal.spec_covar import spectral_covariance

        r = spectral_covariance(
            args.config,
            args.datacube,
            rho=args.rho,
            width=args.width,
            step=args.step,
            prad=args.prad,
            filt=args.filt,
            normnoise=args.normnoise,
            maskcoord=tuple(args.maskcoord) if args.maskcoord else None,
        )
        print(
            f"A_rho {r.a_rho:.3f}  A_lam {r.a_lam:.3f}  A_del {r.a_del:.3f}  "
            f"sig_rho {r.sig_rho:.3f}  sig_lam {r.sig_lam:.3f}  -> {r.params_file}"
        )
        return 0

    return run


@register("astrom-grid", help="Astrometry from extractions over a grid of starting positions")
def _ag(parser):
    parser.add_argument("config", type=Path)
    parser.add_argument("--datacube", default="final.fits")
    parser.add_argument("--xrange", nargs=2, type=float, required=True)
    parser.add_argument("--yrange", nargs=2, type=float, required=True)
    parser.add_argument("--deltapos", type=float, default=1.0)
    parser.add_argument("--throughput")

    def run(args) -> int:
        from charis_dpp.specphotcal.astrom_grid import astrometry_grid

        r = astrometry_grid(
            args.config,
            args.datacube,
            tuple(args.xrange),
            tuple(args.yrange),
            deltapos=args.deltapos,
            throughput=args.throughput,
        )
        print(
            f"x = {r['x_mean']:.3f} +/- {r['x_std']:.3f}, y = {r['y_mean']:.3f} +/- {r['y_std']:.3f} ({r['n']} fits)"
        )
        return 0

    return run


@register("atmosfit", help="Fit a BT-Settl model grid to an extracted spectrum")
def _atm(parser):
    parser.add_argument(
        "spectrum", type=Path, help="wavelength flux err file (extract-spectrum output)"
    )
    parser.add_argument("--grid-dir", type=Path, help="default: $CHARIS_DPP_BTSETTL")
    parser.add_argument("--covar", type=Path, help="covar_params file from spec-covar")
    parser.add_argument("--phot", type=Path, help="photometry file: filter wavelength flux err")
    parser.add_argument("--band", choices=["j", "h", "k", "hk", "jhk"])
    parser.add_argument("--distance", type=float, default=10.0, help="distance [pc]")
    parser.add_argument("--outdir", default="btsettl_fit")
    parser.add_argument("--plot", action="store_true")

    def run(args) -> int:
        from charis_dpp.tools.atmosfit import fit_model_grid

        res = fit_model_grid(
            args.spectrum,
            grid_dir=args.grid_dir,
            covar=args.covar,
            phot=args.phot,
            band=args.band,
            distance_pc=args.distance,
            outdir=args.outdir,
            plot=args.plot,
        )
        for m in res[:10]:
            print(
                f"{m.name:40s} Teff {m.teff:5.0f} logg {m.logg:.1f} chi2_red {m.reduced:7.3f} R {m.radius_rjup:.2f} RJ"
            )
        return 0

    return run


@register("empcomp", help="Compare an extracted spectrum with the Montreal/Bonnefoy libraries")
def _emp(parser):
    parser.add_argument("spectrum", type=Path)
    parser.add_argument("--library", choices=["montreal", "bonnefoy"], default="montreal")
    parser.add_argument("--library-dir", type=Path)
    parser.add_argument("--covar", type=Path)
    parser.add_argument("--band", choices=["j", "h", "k", "hk", "jhk"])
    parser.add_argument("--outdir")
    parser.add_argument("--plot", action="store_true")

    def run(args) -> int:
        from charis_dpp.tools.empcomp import fit_empirical_library

        res = fit_empirical_library(
            args.spectrum,
            library=args.library,
            library_dir=args.library_dir,
            covar=args.covar,
            band=args.band,
            outdir=args.outdir,
            plot=args.plot,
        )
        for m in res[:10]:
            print(f"{m.name:30s} {m.spt:12s} grav {m.gravity} chi2_red {m.reduced:7.3f}")
        return 0

    return run


@register("prop-mot", help="Background-object track from the star's proper motion and parallax")
def _pm(parser):
    parser.add_argument(
        "config", type=Path, help="dataset .toml (target coordinates, parallax, proper motion)"
    )
    parser.add_argument("--start", required=True, help="start date, e.g. 2018-10-18")
    parser.add_argument("--nyears", type=float, default=1.0)
    parser.add_argument("--position", nargs=2, type=float, default=(0.0, 0.0), metavar=("E", "N"))
    parser.add_argument("--posfile", type=Path, help="measured positions: e n sig_e sig_n date")
    parser.add_argument("--plot", default="prop_mot.png")

    def run(args) -> int:
        from charis_dpp.config.dataset import DatasetConfig
        from charis_dpp.io.catalogs import StarInfo
        from charis_dpp.tools.astrometry import background_track

        cfg = DatasetConfig.load(args.config)
        t = cfg.target
        star = StarInfo(
            t.get("name", ""),
            float(t["ra_deg"]),
            float(t["dec_deg"]),
            plx_mas=float(t.get("plx_mas", float("nan"))),
            pmra_mas_yr=float(t.get("pmra_mas_yr", 0.0)),
            pmdec_mas_yr=float(t.get("pmdec_mas_yr", 0.0)),
        )
        tr = background_track(
            star,
            args.start,
            nyears=args.nyears,
            position=tuple(args.position),
            posfile=args.posfile,
            plot=args.plot,
        )
        print(
            f"{tr.times.size} days; final offset ({tr.d_east_arcsec[-1]:+.3f}, {tr.d_north_arcsec[-1]:+.3f}) arcsec"
        )
        if tr.chisq is not None:
            print(f"background hypothesis chi^2 = {tr.chisq:.2f}")
        return 0

    return run


@register("encircled-energy", help="Encircled-energy curve of a PSF cube")
def _ee(parser):
    parser.add_argument("cube", type=Path, nargs="?", default=Path("psfmodel/psfcube_med.fits"))
    parser.add_argument("--position", nargs=2, type=float)
    parser.add_argument("--dpsf", type=int, default=41)
    parser.add_argument("--channel", type=int)
    parser.add_argument("--prad", action="store_true")
    parser.add_argument("--background", type=float, default=0.0)
    parser.add_argument("--plot", default="encircled_energy.png")

    def run(args) -> int:
        from charis_dpp.tools.imagequality import encircled_energy

        r = encircled_energy(
            args.cube,
            position=tuple(args.position) if args.position else None,
            dpsf=args.dpsf,
            channel=args.channel,
            prad=args.prad,
            background=args.background,
            plot=args.plot,
        )
        il = r.ee.shape[0] // 2 if args.channel is None else args.channel
        for rad, e in zip(r.radii_px[::4], r.ee[il][::4], strict=True):
            print(f"r = {rad:5.1f} px  EE = {e:.3f}")
        return 0

    return run


@register("resample-disk-model", help="Resample a disk model image to the CHARIS pixel scale")
def _rdm(parser):
    parser.add_argument("model", type=Path)
    parser.add_argument("out", type=Path)
    parser.add_argument(
        "--input-pixscale", type=float, required=True, help="arcsec/px of the model"
    )
    parser.add_argument(
        "--output-pixscale", type=float, help="default: instrument pixel scale (today)"
    )
    parser.add_argument("--size", type=int, help="output size in pixels (default: input)")

    def run(args) -> int:
        from charis_dpp.config.instrument import instrument_for
        from charis_dpp.io.fits import write_image
        from charis_dpp.tools.diskmodel_resample import resample_disk_model

        out_scale = args.output_pixscale
        if out_scale is None:
            from astropy.time import Time

            out_scale = instrument_for(Time.now().mjd).pixscale_arcsec
        im = resample_disk_model(args.model, args.input_pixscale, out_scale, out_size=args.size)
        write_image(args.out, im)
        print(f"wrote {args.out} ({im.shape[-1]} px at {out_scale} arcsec/px)")
        return 0

    return run


@register(
    "compare-fits", help="Compare two FITS products (relative RMS, correlation, finite agreement)"
)
def _cmp(parser):
    parser.add_argument("a", type=Path)
    parser.add_argument("b", type=Path)
    parser.add_argument("--border", type=int, default=2)
    parser.add_argument("--ext", type=int)

    def run(args) -> int:
        from charis_dpp.testing.compare import compare_fits

        r = compare_fits(args.a, args.b, border=args.border, ext=args.ext)
        print(r)
        return 0 if r.ok() else 1

    return run
