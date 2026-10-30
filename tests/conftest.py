from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="session")
def repo_root() -> Path:
    return REPO


@pytest.fixture(scope="session")
def hr8799_info() -> Path:
    return REPO / "idl" / "prelims" / "HR8799_low.info"


@pytest.fixture
def rng() -> np.random.Generator:
    return np.random.default_rng(12345)


@pytest.fixture(scope="session")
def synthetic_dataset(tmp_path_factory) -> Path:
    """A small synthetic dataset (session-scoped; tests must not modify raw/)."""
    from charis_dpp.synthetic import SyntheticConfig, make_dataset

    out = tmp_path_factory.mktemp("synth")
    cfg = SyntheticConfig(nframes=8, shape=(101, 101), n_speckles=12, nsky=2, seed=3)
    make_dataset(out, cfg)
    return out


def pytest_collection_modifyitems(config, items):
    if not os.environ.get("CHARIS_DPP_TESTDATA"):
        skip = pytest.mark.skip(reason="set CHARIS_DPP_TESTDATA to run regression tests")
        for item in items:
            if "regression" in item.keywords:
                item.add_marker(skip)
