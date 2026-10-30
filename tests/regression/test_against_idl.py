"""Reduce an archived dataset with the Python pipeline and compare with the IDL products."""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest
from charis_dpp.config.dataset import DatasetConfig
from charis_dpp.testing.compare import compare_fits, compare_reduclog

pytestmark = pytest.mark.regression


def _config(d: Path) -> Path:
    tomls = sorted(d.glob("*.toml"))
    if tomls:
        return tomls[0]
    infos = sorted(d.glob("*.info"))
    if not infos:
        raise FileNotFoundError(f"{d}: no dataset .toml or .info")
    return DatasetConfig.from_info(infos[0]).save(d / (infos[0].stem + ".toml"))


def test_registration_and_calibration(dataset_dir: Path, tmp_path: Path, tolerances: dict):
    from charis_dpp.prelims.imprep import imprep
    from charis_dpp.reg import register_cube
    from charis_dpp.skysub import subtract_sky
    from charis_dpp.specphotcal import specphot_cal

    work = tmp_path / dataset_dir.name
    work.mkdir()
    shutil.copytree(dataset_dir / "raw", work / "data" / "raw")
    cfg_src = _config(dataset_dir)
    shutil.copy(cfg_src, work / cfg_src.name)
    cfg = work / cfg_src.name
    imprep(cfg)
    if DatasetConfig.load(cfg).frames("sky"):
        subtract_sky(cfg)
    register_cube(cfg)
    idl_log = dataset_dir / "idl" / "reduc.log"
    if idl_log.exists():
        d = compare_reduclog(work / "reduc.log", idl_log)
        t = tolerances["reduclog"]
        assert d["max_dx"] <= t["max_dx"] and d["max_dy"] <= t["max_dy"], d
        assert d["max_dpa"] <= t["max_dpa"], d
    specphot_cal(cfg)
    ti = tolerances["images"]
    for idl_file in sorted((dataset_dir / "idl" / "reduc" / "reg").glob("n*reg_cal.fits"))[:5]:
        r = compare_fits(work / "reduc" / "reg" / idl_file.name, idl_file, border=ti["border"])
        assert r.pearson >= ti["pearson"] and r.finite_agree >= ti["finite_agree"], str(r)


def test_final_products(dataset_dir: Path, tolerances: dict):
    """Compare archived Python products (reduc/proc) with the IDL ones when both are present."""
    py = dataset_dir / "python" / "reduc" / "proc"
    idl = dataset_dir / "idl" / "reduc" / "proc"
    if not (py.is_dir() and idl.is_dir()):
        pytest.skip("no archived final products")
    for f in sorted(idl.glob("*_collapsed.fits")):
        if (py / f.name).exists():
            r = compare_fits(py / f.name, f, border=tolerances["images"]["border"])
            assert r.pearson >= tolerances["final"]["pearson"], str(r)
