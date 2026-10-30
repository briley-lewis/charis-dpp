"""``charis adialoci`` / ``adiklip`` / ``sdialoci`` / ``rdiklip`` / ``ssdi``."""

from __future__ import annotations

from pathlib import Path

from charis_dpp.cli.registry import register


def _common(parser):
    parser.add_argument("config", type=Path)
    parser.add_argument("--nfwhm", type=float, help="rotation gap in PSF FWHM")
    parser.add_argument(
        "--drsub", type=float, nargs="+", help="annulus width [px] (or w_in w_out r_break r_scale)"
    )
    parser.add_argument("--rmin", type=float)
    parser.add_argument("--rmax", type=float)
    parser.add_argument("--rsubval", type=int, help="1: use the spatially filtered (rsub) cubes")
    parser.add_argument("--zero", action="store_true")
    parser.add_argument("--meanadd", action="store_true")
    parser.add_argument("--channel", type=int, help="process a single channel")
    parser.add_argument("--fc", action="store_true", help="synthetic-companion (_fc) cubes")
    parser.add_argument(
        "--fwdmod",
        action="store_true",
        help="forward-model the _fc cubes through the data reduction",
    )
    parser.add_argument(
        "--norot", action="store_true", help="do not rotate north-up (for a later SDI step)"
    )
    parser.add_argument("--suffix")
    parser.add_argument("--outfile", help="output cube name in reduc/proc (default final.fits)")


def _kw(args, names):
    kw = {}
    for n in names:
        v = getattr(args, n, None)
        if v is not None and v is not False:
            kw[n] = v
    if kw.get("drsub") is not None and len(kw["drsub"]) == 1:
        kw["drsub"] = kw["drsub"][0]
    return kw


_COMMON = [
    "nfwhm",
    "drsub",
    "rmin",
    "rmax",
    "rsubval",
    "zero",
    "meanadd",
    "channel",
    "fc",
    "fwdmod",
    "norot",
    "suffix",
    "outfile",
]


@register("adialoci", help="ADI PSF subtraction with A-LOCI")
def _adialoci(parser):
    _common(parser)
    parser.add_argument("--na", type=float, help="optimisation area in PSF cores")
    parser.add_argument("--geom", type=float)
    parser.add_argument("--svd", type=float, help="SVD cutoff (absolute)")
    parser.add_argument("--nref", type=int)
    parser.add_argument("--pixmask", action="store_true")
    parser.add_argument("--postsdi", action="store_true")

    def run(args) -> int:
        from charis_dpp.psfsub import adialoci

        out, col = adialoci(
            args.config, **_kw(args, _COMMON + ["na", "geom", "svd", "nref", "pixmask", "postsdi"])
        )
        print(f"wrote {out} and {col}")
        return 0

    return run


@register("adiklip", help="ADI PSF subtraction with KLIP")
def _adiklip(parser):
    _common(parser)
    parser.add_argument("--npca", type=int)
    parser.add_argument("--nzone", type=int)
    parser.add_argument("--nomeansub", action="store_true")

    def run(args) -> int:
        from charis_dpp.psfsub import adiklip

        kw = _kw(args, _COMMON + ["npca", "nzone"])
        if args.nomeansub:
            kw["meansub"] = False
        out, col = adiklip(args.config, **kw)
        print(f"wrote {out} and {col}")
        return 0

    return run


@register("sdialoci", help="SDI PSF subtraction with A-LOCI")
def _sdialoci(parser):
    _common(parser)
    parser.add_argument("--na", type=float)
    parser.add_argument("--geom", type=float)
    parser.add_argument("--svd", type=float)
    parser.add_argument("--nref", type=int)
    parser.add_argument("--pixmask", action="store_true")
    parser.add_argument("--postadi", action="store_true")

    def run(args) -> int:
        from charis_dpp.psfsub.sdialoci import sdialoci

        out, col = sdialoci(
            args.config, **_kw(args, _COMMON + ["na", "geom", "svd", "nref", "pixmask", "postadi"])
        )
        print(f"wrote {out} and {col}")
        return 0

    return run


@register("rdiklip", help="RDI PSF subtraction with KLIP using a reference-star library")
def _rdiklip(parser):
    _common(parser)
    parser.add_argument(
        "--refpath", type=Path, required=True, help="reference dataset directory (with reduc/)"
    )
    parser.add_argument("--npca", type=int)
    parser.add_argument(
        "--psfrefsuffix", help="suffix of the reference cubes (default: same as science)"
    )
    parser.add_argument("--nomeansub", action="store_true")

    def run(args) -> int:
        from charis_dpp.psfsub.rdiklip import rdiklip

        kw = _kw(args, _COMMON + ["npca"])
        if args.nomeansub:
            kw["meansub"] = False
        out, col = rdiklip(args.config, refpath=args.refpath, psfrefsuffix=args.psfrefsuffix, **kw)
        print(f"wrote {out} and {col}")
        return 0

    return run


@register("ssdi", help="Simple classical SDI of PSF-subtracted (or filtered) cubes")
def _ssdi(parser):
    parser.add_argument("config", type=Path)
    parser.add_argument("--source", choices=["aloci", "klip", "rsub"], default="aloci")
    parser.add_argument("--weight", action="store_true")
    parser.add_argument("--filt", action="store_true")
    parser.add_argument("--northup", action="store_true")

    def run(args) -> int:
        from charis_dpp.psfsub.ssdi import ssdi

        out = ssdi(
            args.config,
            source=args.source,
            weight=args.weight,
            filt=args.filt,
            northup_final=args.northup,
        )
        print(f"wrote {out}")
        return 0

    return run
