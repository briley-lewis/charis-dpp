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
    cfg = SyntheticConfig(
        nframes=8, shape=(101, 101), n_speckles=12, nsky=2, seed=3, astrogrid_lambda_d=11.2
    )
    make_dataset(out, cfg)
    return out


@pytest.fixture(scope="session")
def reduced_dataset(synthetic_dataset: Path, tmp_path_factory):
    """Synthetic dataset reduced through registration, calibration, filtering and A-LOCI."""
    import json
    import shutil

    from charis_dpp.filt import imrsub
    from charis_dpp.prelims.imprep import imprep
    from charis_dpp.psfsub import adialoci
    from charis_dpp.reg import register_cube
    from charis_dpp.skysub import subtract_sky
    from charis_dpp.specphotcal import make_empirical_psf, specphot_cal

    d = tmp_path_factory.mktemp("reduced")
    shutil.copytree(synthetic_dataset, d, dirs_exist_ok=True)
    cfg = next(d.glob("*.toml"))
    imprep(cfg)
    subtract_sky(cfg)
    register_cube(cfg)
    specphot_cal(cfg, calmethod=1, starlib="pickles")
    make_empirical_psf(cfg)
    imrsub(cfg, prad=True)
    adialoci(
        cfg,
        rsubval=1,
        nfwhm=0.75,
        drsub=5,
        na=100,
        svd=1e-6,
        nref=0,
        rmin=6,
        rmax=40,
        outfile="aloci.fits",
    )
    truth = json.loads((d / "truth.json").read_text())
    return d, cfg, truth


def pytest_collection_modifyitems(config, items):
    if not os.environ.get("CHARIS_DPP_TESTDATA"):
        skip = pytest.mark.skip(reason="set CHARIS_DPP_TESTDATA to run regression tests")
        for item in items:
            if "regression" in item.keywords:
                item.add_marker(skip)
