"""Locations of package data and of the large external model libraries."""

from __future__ import annotations

import os
from importlib import resources
from pathlib import Path

_PKG_DATA = resources.files("charis_dpp") / "data"


def package_data(relpath: str) -> Path:
    """Absolute path of a file shipped under ``charis_dpp/data``."""
    p = Path(str(_PKG_DATA / relpath))
    if not p.exists():
        raise FileNotFoundError(f"package data file not found: {relpath}")
    return p


def repo_root() -> Path | None:
    """Repository root when running from a source checkout, else ``None``."""
    here = Path(__file__).resolve()
    for parent in here.parents:
        if (parent / "idl").is_dir() and (parent / "pyproject.toml").exists():
            return parent
    return None


def models_dir() -> Path:
    """Directory holding ``planetmodels/`` and ``starmodels/``.

    ``$CHARIS_DPP_MODELS`` wins; otherwise ``<repo>/idl/models`` from a source
    checkout.  Raises a clear error when neither exists (``charis doctor``
    reports it).
    """
    env = os.environ.get("CHARIS_DPP_MODELS")
    if env:
        p = Path(env).expanduser()
        if p.is_dir():
            return p
        raise FileNotFoundError(f"$CHARIS_DPP_MODELS={env!r} is not a directory")
    root = repo_root()
    if root is not None and (root / "idl" / "models").is_dir():
        return root / "idl" / "models"
    raise FileNotFoundError(
        "model libraries not found: set $CHARIS_DPP_MODELS to a directory containing "
        "planetmodels/ and starmodels/ (the idl/models tree of the repository)"
    )


def btsettl_dir() -> Path:
    """BT-Settl model grid (``$CHARIS_DPP_BTSETTL``); not shipped."""
    env = os.environ.get("CHARIS_DPP_BTSETTL")
    if not env:
        raise FileNotFoundError("set $CHARIS_DPP_BTSETTL to the BT-Settl lte*.dat.txt directory")
    p = Path(env).expanduser()
    if not p.is_dir():
        raise FileNotFoundError(f"$CHARIS_DPP_BTSETTL={env!r} is not a directory")
    return p
