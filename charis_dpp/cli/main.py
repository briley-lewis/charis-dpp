"""Entry point for the ``charis`` command."""

from __future__ import annotations

import argparse
import importlib
import logging
import sys
from datetime import datetime
from pathlib import Path

from charis_dpp import __version__
from charis_dpp.cli.registry import all_subcommands

# Modules that register subcommands. Import failures (e.g. an optional
# dependency missing) must not take the whole CLI down, so each is guarded.
_COMMAND_MODULES = [
    "charis_dpp.cli.commands.doctor",
    "charis_dpp.cli.commands.config",
    "charis_dpp.cli.commands.prelims",
    "charis_dpp.cli.commands.skysub",
    "charis_dpp.cli.commands.reg",
    "charis_dpp.cli.commands.specphotcal",
    "charis_dpp.cli.commands.filt",
    "charis_dpp.cli.commands.psfsub",
    "charis_dpp.cli.commands.partialsub",
    "charis_dpp.cli.commands.pdi",
    "charis_dpp.cli.commands.tools",
    "charis_dpp.cli.commands.synthetic",
    "charis_dpp.cli.commands.autoreduce",
]


def _load_commands() -> list[str]:
    failed = []
    for mod in _COMMAND_MODULES:
        try:
            importlib.import_module(mod)
        except ModuleNotFoundError as exc:  # pragma: no cover - partial installs
            if exc.name and exc.name.startswith("charis_dpp"):
                continue  # not implemented yet
            failed.append(f"{mod}: {exc}")
    return failed


def _setup_logging(workdir: Path, level: str) -> None:
    fmt = "%(asctime)s %(levelname)s %(name)s: %(message)s"
    handlers: list[logging.Handler] = [logging.StreamHandler(sys.stderr)]
    logdir = workdir / "log"
    try:
        logdir.mkdir(exist_ok=True)
        handlers.append(logging.FileHandler(logdir / f"charis_{datetime.now():%Y%m%d}.log"))
    except OSError:
        pass
    logging.basicConfig(level=getattr(logging, level.upper()), format=fmt, handlers=handlers)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="charis",
        description="CHARIS Data Processing Pipeline",
    )
    parser.add_argument("--version", action="version", version=f"charis-dpp {__version__}")
    parser.add_argument("--workdir", type=Path, default=Path("."), help="dataset directory")
    parser.add_argument("--log-level", default="INFO", help="DEBUG, INFO, WARNING, ERROR")
    sub = parser.add_subparsers(dest="command", metavar="<command>")
    failed = _load_commands()
    for name, cmd in sorted(all_subcommands().items()):
        p = sub.add_parser(name, help=cmd.help, description=cmd.help)
        p.set_defaults(_run=cmd.add(p))
    if failed:
        parser.epilog = "Some commands failed to load:\n  " + "\n  ".join(failed)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command is None:
        parser.print_help()
        return 2
    _setup_logging(args.workdir, args.log_level)
    rc = args._run(args)
    return int(rc or 0)


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
