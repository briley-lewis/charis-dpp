"""Subcommand registry for the ``charis`` CLI.

A module registers a subcommand with::

    from charis_dpp.cli.registry import register

    @register("register", help="Register cubes to a common centre")
    def _add(parser):
        parser.add_argument("pfname")
        ...
        return run   # callable(args) -> int | None
"""

from __future__ import annotations

import argparse
from collections.abc import Callable
from dataclasses import dataclass

Runner = Callable[[argparse.Namespace], int | None]
Adder = Callable[[argparse.ArgumentParser], Runner]


@dataclass(frozen=True)
class Subcommand:
    name: str
    help: str
    add: Adder
    group: str = "pipeline"


_REGISTRY: dict[str, Subcommand] = {}


def register(name: str, *, help: str, group: str = "pipeline") -> Callable[[Adder], Adder]:
    def deco(add: Adder) -> Adder:
        if name in _REGISTRY:
            raise ValueError(f"subcommand {name!r} registered twice")
        _REGISTRY[name] = Subcommand(name, help, add, group)
        return add

    return deco


def all_subcommands() -> dict[str, Subcommand]:
    return dict(_REGISTRY)
