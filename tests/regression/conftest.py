"""Regression tests run only with ``CHARIS_DPP_TESTDATA=/path`` pointing at datasets with archived IDL outputs."""

from __future__ import annotations

import os
import tomllib
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent


def datasets() -> list[Path]:
    base = os.environ.get("CHARIS_DPP_TESTDATA")
    if not base:
        return []
    return sorted(p for p in Path(base).iterdir() if (p / "raw").is_dir())


@pytest.fixture(scope="session")
def tolerances() -> dict:
    with open(ROOT / "tolerances.toml", "rb") as f:
        return tomllib.load(f)


def pytest_generate_tests(metafunc):
    if "dataset_dir" in metafunc.fixturenames:
        ds = datasets()
        metafunc.parametrize("dataset_dir", ds, ids=[d.name for d in ds])
