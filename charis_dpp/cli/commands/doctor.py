"""``charis doctor`` — check dependencies, data tables and model libraries."""

from __future__ import annotations

import importlib

from charis_dpp.cli.registry import register


@register(
    "doctor",
    help="Check the installation: dependencies, instrument table, model libraries",
    group="misc",
)
def _add(parser):
    def run(args) -> int:
        ok = True
        for mod in ("numpy", "scipy", "astropy", "photutils", "astroquery", "matplotlib"):
            try:
                m = importlib.import_module(mod)
                print(f"  ok  {mod} {getattr(m, '__version__', '')}")
            except ImportError as exc:
                ok = False
                print(f" FAIL {mod}: {exc}")
        from charis_dpp.config.instrument import load_instrument

        try:
            t = load_instrument()
            print(
                f"  ok  instrument table {t.path} ({len(t.epochs)} epochs: {', '.join(t.epoch_labels())})"
            )
        except Exception as exc:
            ok = False
            print(f" FAIL instrument table: {exc}")
        from charis_dpp.io.wavelengths import wavelengths_um

        for band in ("lowres", "J", "H", "K"):
            try:
                w = wavelengths_um(band)
                print(f"  ok  {band}: {w.size} channels {w[0]:.3f}-{w[-1]:.3f} um")
            except Exception as exc:
                ok = False
                print(f" FAIL wavelength table {band}: {exc}")
        from charis_dpp.config.paths import btsettl_dir, models_dir

        try:
            print(f"  ok  model libraries at {models_dir()}")
        except FileNotFoundError as exc:
            print(f" warn {exc}")
        try:
            print(f"  ok  BT-Settl grid at {btsettl_dir()}")
        except FileNotFoundError:
            print(
                "  --  BT-Settl grid not configured ($CHARIS_DPP_BTSETTL); only needed for atmosfit"
            )
        print("CHARIS-DPP installation " + ("looks good" if ok else "has problems"))
        return 0 if ok else 1

    return run
