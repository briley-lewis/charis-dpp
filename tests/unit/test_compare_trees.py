"""``compare_trees`` on small synthetic reduction trees."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pytest
from charis_dpp.cli.main import main
from charis_dpp.io.fits import write_image
from charis_dpp.io.reduclog import ReducLog
from charis_dpp.testing.compare_trees import (
    compare_hwp_info,
    compare_trees,
    default_tolerances,
    pair_files,
    read_hwp_info,
)

SHAPE = (3, 21, 21)


def _cube(rng, seed_offset=0):
    r = np.random.default_rng(rng.integers(1 << 30) + seed_offset)
    return r.normal(1.0, 0.1, SHAPE)


def _tree(root: Path, cubes: dict[str, np.ndarray], *, idl_hwp: bool = False) -> Path:
    for rel, data in cubes.items():
        write_image(root / rel, data, dtype=np.float64)
    log = ReducLog()
    for n in (1, 2, 3):
        log.add(n, 100.0 + n * 0.1, 99.5, 30.0, -0.5 + n * 0.01, 10.0 * n)
    log.save(root / "reduc.log")
    if idl_hwp:
        lines = ["FILE, HWP_POS, HWP_ANG, CYCLE, "]
        lines += [f"{n} {(n - 1) % 4} {22.5 * ((n - 1) % 4):.1f} 0" for n in (1, 2, 3)]
    else:
        lines = ["# name hwp_pos hwp_ang hwp_cycle"]
        lines += [
            f"n{n:04d}leftreg_cal.fits {(n - 1) % 4} {22.5 * ((n - 1) % 4):.3f} 0"
            for n in (1, 2, 3)
        ]
    (root / "hwp_info.txt").write_text("\n".join(lines) + "\n")
    return root


@pytest.fixture
def trees(tmp_path, rng):
    cubes = {
        "prep/n0001.fits": _cube(rng),
        "reg/n0001reg_cal.fits": _cube(rng),
        "proc/final.fits": _cube(rng),
        "proc/final_collapsed.fits": _cube(rng)[0],
    }
    a = _tree(tmp_path / "a", cubes)
    b = _tree(tmp_path / "b", cubes, idl_hwp=True)
    return a, b, cubes


def test_identical_trees_pass(trees):
    a, b, cubes = trees
    res = compare_trees(a, b, "identical", default_tolerances())
    assert res.ok and res.n_missing == 0
    assert {f.relpath for f in res.files} == set(cubes) | {"reduc.log", "hwp_info.txt"}
    reg = next(f for f in res.files if f.relpath == "reg/n0001reg_cal.fits")
    assert reg.metrics["rel_rms"] == 0.0 and reg.metrics["finite_agree"] == 1.0
    assert len(reg.channels) == SHAPE[0]
    assert next(f for f in res.files if f.kind == "hwp_info").passed  # IDL vs Python layout


def test_perturbed_file_fails_identical_passes_port(trees, tmp_path):
    a, b, cubes = trees
    rel = "reg/n0001reg_cal.fits"
    pert = cubes[rel] * (1 + 1e-4)  # rel_rms 1e-4: above 1e-6, below the port 1e-3
    write_image(b / rel, pert, dtype=np.float64)
    strict = compare_trees(a, b, "identical", default_tolerances())
    assert not strict.ok and strict.n_fail == 1
    bad = next(f for f in strict.files if f.failed)
    assert bad.relpath == rel and bad.rule == "identical"
    assert 5e-5 < bad.metrics["rel_rms"] < 2e-4
    port = compare_trees(a, b, "port", default_tolerances())
    assert port.ok
    assert next(f for f in port.files if f.relpath == rel).rule == "port"
    # the same file on an allow-list is judged with the port thresholds in identical mode
    allow = tmp_path / "allow.txt"
    allow.write_text("# known difference\nreg/*reg_cal.fits\n")
    assert compare_trees(a, b, "identical", default_tolerances(), allow_list=allow).ok
    # a collapsed image uses [final] (correlation only) in port mode
    write_image(b / "proc/final_collapsed.fits", cubes["proc/final_collapsed.fits"] * 1.5)
    res = compare_trees(a, b, "port", default_tolerances())
    fin = next(f for f in res.files if f.relpath == "proc/final_collapsed.fits")
    assert fin.rule == "port/final" and fin.passed
    # a genuinely different image fails even in port mode
    write_image(b / "proc/final.fits", np.random.default_rng(1).normal(1, 0.1, SHAPE))
    assert not compare_trees(a, b, "port", default_tolerances()).ok


def test_missing_and_nan_files(trees):
    a, b, cubes = trees
    (b / "prep/n0001.fits").unlink()
    write_image(a / "proc/extra.fits", cubes["proc/final.fits"])
    nan = cubes["proc/final.fits"].copy()
    nan[1, 10, 10] = np.nan
    write_image(b / "proc/final.fits", nan, dtype=np.float64)
    res = compare_trees(a, b, "identical", default_tolerances())
    by = {f.relpath: f for f in res.files}
    assert (
        by["prep/n0001.fits"].status == "missing" and "only in A" in by["prep/n0001.fits"].message
    )
    assert (
        by["proc/extra.fits"].status == "missing" and "only in A" in by["proc/extra.fits"].message
    )
    assert res.n_missing == 2
    assert by["proc/final.fits"].status == "fail"  # finite_agree < 1
    assert by["proc/final.fits"].metrics["finite_agree"] < 1.0
    assert res.n_fail == 1  # missing files are not failures


def test_hwp_info_rows_and_layouts(tmp_path):
    py = tmp_path / "py.txt"
    py.write_text(
        "# name hwp_pos hwp_ang hwp_cycle\nn0001leftreg_cal.fits 0 0.000 0\nn0002leftreg_cal.fits 1 22.500 0\nn0003leftreg_cal.fits 2 45.000 1\n"
    )
    idl = tmp_path / "idl.txt"
    idl.write_text(
        "FILE, HWP_POS, HWP_ANG, CYCLE, \n       1       0      0.00000       0\n       2       1      22.5000       0\n       3       3      67.5000       1\n       4       0      0.00000       2\n"
    )
    assert read_hwp_info(py)[2] == (1, 22.5, 0)
    assert read_hwp_info(idl)[4] == (0, 0.0, 2)
    d = compare_hwp_info(py, idl)
    assert d["n"] == 3 and d["n_only_b"] == 1 and d["n_only_a"] == 0
    assert d["n_pos_diff"] == 1 and d["n_cycle_diff"] == 0 and d["max_dang"] == 22.5
    assert [x[0] for x in d["diffs"]] == [3]


def test_reduc_subdir_layout_and_reports(tmp_path, rng):
    cubes = {"reduc/reg/n0001reg_cal.fits": _cube(rng)}
    a = _tree(tmp_path / "a", cubes)
    b = _tree(tmp_path / "b", cubes)
    assert pair_files(a, b) == ["reduc/reg/n0001reg_cal.fits", "reduc.log", "hwp_info.txt"]
    out = tmp_path / "out"
    wd = ["--workdir", str(tmp_path)]  # keeps the CLI log out of the repository
    rc = main([*wd, "compare-trees", str(a), str(b), "--mode", "identical", "--out", str(out)])
    assert rc == 0
    rows = list(csv.DictReader(open(out / "compare.csv")))
    whole = [r for r in rows if r["channel"] == "all"]
    assert {r["relpath"] for r in whole} == set(pair_files(a, b))
    assert sum(r["relpath"] == "reduc/reg/n0001reg_cal.fits" for r in rows) == 1 + SHAPE[0]
    md = (out / "compare.md").read_text()
    assert "| `reduc.log` | pass |" in md and "3 file(s): 3 pass, 0 fail, 0 missing" in md
    # perturbation beyond the port tolerance -> exit code 1
    write_image(b / "reduc/reg/n0001reg_cal.fits", cubes["reduc/reg/n0001reg_cal.fits"] + 0.5)
    assert main([*wd, "compare-trees", str(a), str(b), "--mode", "port", "--out", str(out)]) == 1
    with pytest.raises(ValueError):
        compare_trees(a, b, "bogus")
