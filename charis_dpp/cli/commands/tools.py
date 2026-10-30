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
