"""``charis register`` / ``psf-quality``."""

from __future__ import annotations

from pathlib import Path

from charis_dpp.cli.registry import register


@register(
    "register", help="Register every slice of every cube to the frame centre; writes reduc.log"
)
def _add(parser):
    parser.add_argument("config", type=Path)
    parser.add_argument(
        "--method",
        default=None,
        help="sats (default: satellite spots) | psf (unsaturated star) | calsat (ladder reference) | xcorr (cross-correlation, no spots needed)",
    )
    parser.add_argument(
        "--xcorr-center",
        choices=["symmetry", "centroid", "given"],
        help="method=xcorr: how the reference frame's star centre is found (default symmetry)",
    )
    parser.add_argument(
        "--guessoffsets",
        nargs=2,
        type=float,
        metavar=("DX", "DY"),
        help="offset of the star from the frame centre [px]",
    )
    parser.add_argument(
        "--astrogrid", type=float, help="satellite-spot spacing in lambda/D (default: nominal 15.9)"
    )
    parser.add_argument("--refslice", type=int)
    parser.add_argument("--medbox", type=int)
    parser.add_argument("--fittrim", type=int)
    parser.add_argument("--smoothval", type=int)
    parser.add_argument("--fwhmlim", type=float)
    parser.add_argument(
        "--rsub", action="store_true", help="subtract a radial profile before centroiding"
    )
    parser.add_argument(
        "--fitgauss", action="store_true", help="2-D Gaussian fit instead of gcntrd"
    )
    parser.add_argument(
        "--xcorr", action="store_true", help="symmetry-based centre instead of spot centroids"
    )
    parser.add_argument("--smallsteps", action="store_true")
    parser.add_argument("--revise", action="store_true")
    parser.add_argument("--keepprev", action="store_true")
    parser.add_argument("--splitpsf", action="store_true")
    parser.add_argument("--nosmooth", action="store_true")
    parser.add_argument("--refcube", type=int, default=0)
    parser.add_argument(
        "--rin-corr",
        type=float,
        help="calsat: inner radius of the cross-correlation annulus [arcsec]",
    )
    parser.add_argument("--rout-corr", type=float, help="calsat: outer radius [arcsec]")
    parser.add_argument(
        "--ladder", action="store_true", help="register the ladder (spot calibration) frames"
    )
    parser.add_argument("--suffix", help="input suffix (default: e_skysub if present, else e)")
    parser.add_argument("--smask", type=Path, help="text file of x y positions of sources to mask")
    parser.add_argument("--checkquality", action="store_true")

    def run(args) -> int:
        import numpy as np

        from charis_dpp.reg import check_psf_quality, register_cube

        kw = {}
        for k in (
            "method",
            "refslice",
            "medbox",
            "fittrim",
            "smoothval",
            "fwhmlim",
            "refcube",
            "suffix",
        ):
            v = getattr(args, k)
            if v is not None:
                kw[k] = v
        for k in (
            "rsub",
            "fitgauss",
            "xcorr",
            "smallsteps",
            "revise",
            "keepprev",
            "splitpsf",
            "nosmooth",
            "ladder",
        ):
            if getattr(args, k):
                kw[k] = True
        if args.xcorr_center:
            kw["xcorr_center"] = args.xcorr_center
        if args.guessoffsets:
            kw["guessoffsets"] = tuple(args.guessoffsets)
        if args.astrogrid is not None:
            kw["astrogrid"] = args.astrogrid
        if args.rin_corr is not None:
            kw["rin_corr"] = args.rin_corr
        if args.rout_corr is not None:
            kw["rout_corr"] = args.rout_corr
        if args.smask:
            kw["smask"] = [tuple(r) for r in np.loadtxt(args.smask, ndmin=2)[:, :2]]
        log = register_cube(args.config, **kw)
        print(f"registered {len(log)} cubes -> {log.path}")
        if args.checkquality:
            check_psf_quality(args.config)
        return 0

    return run


@register("psf-quality", help="Satellite-spot / halo ratio per registered cube")
def _pq(parser):
    parser.add_argument("config", type=Path)
    parser.add_argument("--all-slices", action="store_true")
    parser.add_argument("--suffix", default="reg")

    def run(args) -> int:
        from charis_dpp.reg import check_psf_quality

        r = check_psf_quality(args.config, suffix=args.suffix, all_slices=args.all_slices)
        print(f"median sat/halo = {float(__import__('numpy').nanmedian(r)):.3f}")
        return 0

    return run
