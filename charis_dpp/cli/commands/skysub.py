"""``charis subtract-sky``."""

from __future__ import annotations

from pathlib import Path

from charis_dpp.cli.registry import register


@register("subtract-sky", help="Subtract a master sky cube from the science cubes")
def _add(parser):
    parser.add_argument("config", type=Path)
    parser.add_argument(
        "--scalesky",
        action="store_true",
        help="scale by the outer-region median of the reddest channel",
    )
    parser.add_argument(
        "--scaleindchan", action="store_true", help="scale each channel independently"
    )
    parser.add_argument(
        "--skyrad", type=float, help="radius beyond which the sky level is measured [px]"
    )
    parser.add_argument(
        "--skylim", type=float, help="minimum level (cts/s) for per-channel scaling"
    )
    parser.add_argument("--ladder", action="store_true")
    parser.add_argument(
        "--sky-files",
        nargs="+",
        type=Path,
        help="explicit sky cubes (instead of observation.frames.sky)",
    )
    parser.add_argument("--suffix", help="input suffix (default: e, or e_flat if PDI flats exist)")
    parser.add_argument("--save-master", action="store_true")

    def run(args) -> int:
        from charis_dpp.skysub import subtract_sky

        out = subtract_sky(
            args.config,
            scalesky=args.scalesky,
            scaleindchan=args.scaleindchan,
            skyrad=args.skyrad,
            skylim=args.skylim,
            ladder=args.ladder,
            sky_files=args.sky_files,
            suffix=args.suffix,
            save_master=args.save_master,
        )
        print(f"wrote {len(out)} sky-subtracted cubes")
        return 0

    return run
