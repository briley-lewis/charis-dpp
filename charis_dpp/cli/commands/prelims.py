"""``charis newobs`` / ``imprep`` / ``set-files`` / ``add-files`` / ``del-files``."""

from __future__ import annotations

import json
from pathlib import Path

from charis_dpp.cli.registry import register


def _cfg_arg(parser):
    parser.add_argument("config", type=Path, help="dataset .toml (or legacy .info)")


@register("newobs", help="Create the directory layout and dataset config for a new observation")
def _newobs(parser):
    parser.add_argument("--date", required=True, help="YYYYMMDD or YYYY-MM-DD")
    parser.add_argument("--object", required=True, dest="object_name")
    parser.add_argument("--filter", required=True, dest="band", help="low/lowres/J/H/K")
    parser.add_argument("--astrogrid", type=float, help="satellite-spot spacing in lambda/D")
    parser.add_argument("--modamp", type=float, help="satellite-spot modulation amplitude [nm]")
    parser.add_argument(
        "--star-json", type=Path, help="offline StarInfo JSON instead of a SIMBAD/VizieR lookup"
    )

    def run(args) -> int:
        from charis_dpp.io.catalogs import StarInfo
        from charis_dpp.prelims.newobs import newobs

        star = StarInfo(**json.loads(args.star_json.read_text())) if args.star_json else None
        p = newobs(
            args.date,
            args.object_name,
            args.band,
            workdir=args.workdir,
            star=star,
            astrogrid_lambda_d=args.astrogrid,
            modulation_nm=args.modamp,
        )
        print(f"wrote {p}")
        return 0

    return run


@register(
    "imprep", help="Normalise raw cube headers, classify frames, write reduc/prep and obslog.txt"
)
def _imprep(parser):
    _cfg_arg(parser)
    parser.add_argument("--fixpix", action="store_true", help="repair bad pixels")
    parser.add_argument("--pdi", action="store_true", help="PDI data (two beams)")
    parser.add_argument("--exptime", type=float, help="override single-read time [s]")
    parser.add_argument("--coadds", type=int, help="override number of reads")

    def run(args) -> int:
        from charis_dpp.prelims.imprep import imprep

        log = imprep(
            args.config,
            fixpix=args.fixpix,
            pdi=args.pdi,
            manual_exptime=args.exptime,
            manual_coadds=args.coadds,
        )
        print(
            f"prepared {len(log.entries)} cubes; science {log.frames('science')}; sky {log.frames('sky')}"
        )
        return 0

    return run


def _frames_cmd(name, help_text, mode):
    @register(name, help=help_text)
    def _add(parser):
        _cfg_arg(parser)
        parser.add_argument("kind", choices=["science", "sky", "ladder"])
        parser.add_argument("frames", help="range spec, e.g. 1-5,8")

        def run(args) -> int:
            from charis_dpp.config.dataset import DatasetConfig
            from charis_dpp.io.filelists import parse_range

            cfg = DatasetConfig.load(args.config)
            cur = set(cfg.frames(args.kind))
            new = set(parse_range(args.frames))
            if mode == "set":
                cur = new
            elif mode == "add":
                cur |= new
            else:
                cur -= new
            cfg.set_frames(args.kind, sorted(cur))
            cfg.save()
            print(f"{args.kind}: {cfg.observation['frames'][args.kind]}")
            return 0

        return run


_frames_cmd("set-files", "Set the science/sky/ladder frame list", "set")
_frames_cmd("add-files", "Add frames to the science/sky/ladder list", "add")
_frames_cmd("del-files", "Remove frames from the science/sky/ladder list", "del")
