"""``charis make-synthetic OUTDIR`` — write a synthetic CHARIS dataset."""

from __future__ import annotations

from pathlib import Path

from charis_dpp.cli.registry import register


@register(
    "make-synthetic",
    help="Write a synthetic CHARIS dataset (raw cubes + config + truth.json)",
    group="misc",
)
def _add(parser):
    parser.add_argument("outdir", type=Path)
    parser.add_argument("--nframes", type=int, default=12)
    parser.add_argument("--size", type=int, default=101, help="frame size in pixels (square)")
    parser.add_argument("--band", default="lowres")
    parser.add_argument("--date", default="2017-08-14")
    parser.add_argument("--pa-step", type=float, default=5.0)
    parser.add_argument("--no-planet", action="store_true")
    parser.add_argument("--planet-sep", type=float, default=22.0)
    parser.add_argument("--planet-pa", type=float, default=60.0)
    parser.add_argument("--planet-contrast", type=float, default=-4.0, help="log10 contrast")
    parser.add_argument("--seed", type=int, default=0)

    def run(args) -> int:
        from charis_dpp.synthetic import Planet, SyntheticConfig, make_dataset

        planet = (
            None
            if args.no_planet
            else Planet(args.planet_sep, args.planet_pa, args.planet_contrast)
        )
        cfg = SyntheticConfig(
            nframes=args.nframes,
            shape=(args.size, args.size),
            band=args.band,
            date=args.date,
            pa_step_deg=args.pa_step,
            planet=planet,
            seed=args.seed,
        )
        p = make_dataset(args.outdir, cfg)
        print(f"wrote synthetic dataset; config: {p}")
        return 0

    return run
