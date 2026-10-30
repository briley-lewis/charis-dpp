"""``charis imrsub``."""

from __future__ import annotations

from pathlib import Path

from charis_dpp.cli.registry import register


@register(
    "imrsub",
    help="Spatial filtering (radial profile / median box / spider mask) + classical ADI products",
)
def _add(parser):
    parser.add_argument("config", type=Path)
    parser.add_argument(
        "--prad", action="store_true", help="subtract a radial profile from each slice"
    )
    parser.add_argument(
        "--pang", action="store_true", help="subtract an azimuthal profile [untested in IDL]"
    )
    parser.add_argument(
        "--psfsub", action="store_true", help="subtract the sequence-median cube before derotating"
    )
    parser.add_argument("--sdi", action="store_true", help="simple SDI before derotating")
    parser.add_argument("--medbox", type=int, help="moving-box median high-pass width [px]")
    parser.add_argument("--mask", action="store_true", help="mask the spider vanes")
    parser.add_argument("--rmax", type=float)
    parser.add_argument(
        "--fc", action="store_true", help="operate on the *_fc (synthetic-companion) cubes"
    )
    parser.add_argument("--nonorthup", action="store_true")
    parser.add_argument("--suffix")

    def run(args) -> int:
        from charis_dpp.filt import imrsub

        out = imrsub(
            args.config,
            prad=args.prad,
            pang=args.pang,
            psfsub=args.psfsub,
            sdi=args.sdi,
            medbox=args.medbox,
            mask=args.mask,
            rmax=args.rmax,
            fc=args.fc,
            northup_final=not args.nonorthup,
            suffix=args.suffix,
        )
        print(f"wrote {len(out)} filtered cubes")
        return 0

    return run
