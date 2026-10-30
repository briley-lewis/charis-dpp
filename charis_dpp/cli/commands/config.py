"""``charis migrate-info``, ``charis show-config``, ``charis constants``."""

from __future__ import annotations

from pathlib import Path

from charis_dpp.cli.registry import register


@register(
    "migrate-info", help="Convert a legacy .info parameter file to a dataset TOML", group="config"
)
def _migrate(parser):
    parser.add_argument("info", type=Path, help="legacy <target>_<band>.info file")
    parser.add_argument("--out", type=Path, help="output .toml (default: same stem)")
    parser.add_argument("--band", help="band if not encoded in the file name (lowres/J/H/K)")

    def run(args) -> int:
        from charis_dpp.config.dataset import DatasetConfig

        cfg = DatasetConfig.from_info(args.info, band=args.band)
        out = args.out or args.info.with_suffix(".toml")
        cfg.save(out, write_info=False)
        print(f"wrote {out}")
        return 0

    return run


@register("show-config", help="Print the resolved configuration for a dataset", group="config")
def _show(parser):
    parser.add_argument("config", type=Path, help="dataset .toml (or legacy .info)")

    def run(args) -> int:
        from charis_dpp.config.dataset import DatasetConfig

        cfg = DatasetConfig.load(args.config)
        r = cfg.resolve()
        inst = r.instrument
        print(f"target        : {cfg.target.get('name')}")
        print(
            f"band          : {r.band} ({r.nlam} channels, {r.wavelengths_um[0]:.3f}-{r.wavelengths_um[-1]:.3f} um)"
        )
        print(
            f"MJD           : {r.mjd:.3f}  epoch={inst.epoch_label} (applied: {', '.join(inst.epochs_applied)})"
        )
        print(
            f"pixscale      : {inst.pixscale_arcsec} arcsec/px   Dtel {inst.dtel_m} m   north offset {inst.north_offset_deg} deg   parity {inst.parity}"
        )
        print(f"ref channel   : {r.reference_channel}   telluric {r.telluric_channels.tolist()}")
        print(
            f"sub-bands     : {{{', '.join(f'{k}: {v.tolist()}' for k, v in r.sub_band_channels.items())}}}"
        )
        print(
            f"frames        : science {cfg.frames('science')[:3]}... ({len(cfg.frames('science'))}), sky {cfg.frames('sky')}, ladder {cfg.frames('ladder')}"
        )
        print(f"lambda/D [px] : {r.fwhm_px.round(2).tolist()}")
        if cfg.data.get("instrument"):
            print(f"overrides     : {cfg.data['instrument']}")
        return 0

    return run


@register("constants", help="Print the instrument constants in force at a date/MJD", group="config")
def _constants(parser):
    parser.add_argument("when", help="MJD, YYYY-MM-DD or YYYYMMDD")

    def run(args) -> int:
        import json

        from charis_dpp.config.instrument import instrument_for

        inst = instrument_for(args.when)
        print(
            f"epoch {inst.epoch_label} (MJD {inst.mjd:.2f}); applied: {', '.join(inst.epochs_applied)}"
        )
        print(json.dumps(inst.values.to_dict(), indent=2, default=str))
        return 0

    return run
