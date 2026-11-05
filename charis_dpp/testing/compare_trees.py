"""File-by-file comparison of two CHARIS-DPP reduction trees.

A *tree* is a dataset working directory (or its ``reduc/`` directory) produced
by any of the pipelines being validated against each other: the pre-port IDL
code, the patched IDL code or the Python port.  Files are paired by relative
path (``prep/*.fits``, ``reg/*.fits``, ``proc/*.fits``, ``reduc*.log``,
``hwp_info.txt``; the three product directories are also looked for under
``reduc/``) and compared with :func:`charis_dpp.testing.compare.compare_fits`
and :func:`charis_dpp.testing.compare.compare_reduclog`.  A file present on
one side only is reported as *missing*, never as an error.

Two modes decide what counts as a pass:

``identical``
    the two trees ran the same code: ``rel_rms <= 1e-6`` and every pixel
    finite on both sides or on neither (``finite_agree == 1``).  Relative
    paths matching a pattern in an *allow-list* file (``fnmatch`` patterns,
    one per line, ``#`` comments) are judged with the ``port`` thresholds
    instead, for products known to differ legitimately.
``port``
    one tree is the Python port: thresholds from
    ``tests/regression/tolerances.toml`` (``[images]`` for FITS products,
    ``[reduclog]`` for the registration logs, ``[final]`` for
    ``proc/*_collapsed.fits``, where speckle noise dominates and only the
    correlation is checked).

Used by ``charis compare-trees``.
"""

from __future__ import annotations

import csv
import fnmatch
import re
import tomllib
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from charis_dpp.io.fits import read_image
from charis_dpp.testing.compare import CompareResult, compare_arrays, compare_reduclog

__all__ = [
    "FileComparison",
    "TreeComparison",
    "compare_hwp_info",
    "compare_trees",
    "default_tolerances",
    "load_allow_list",
    "load_tolerances",
    "pair_files",
    "write_csv",
    "write_markdown",
]

MODES = ("identical", "port")

#: Relative-path patterns that make up a tree, in report order.
PATTERNS = (
    "prep/*.fits",
    "reg/*.fits",
    "proc/*.fits",
    "reduc/prep/*.fits",
    "reduc/reg/*.fits",
    "reduc/proc/*.fits",
    "reduc*.log",
    "reduc/reduc*.log",
    "hwp_info.txt",
    "reduc/hwp_info.txt",
)

IDENTICAL = {"rel_rms": 1e-6, "finite_agree": 1.0}

#: ``tests/regression/tolerances.toml`` as shipped, used when the file is not given.
DEFAULT_TOLERANCES: dict = {
    "reduclog": {"max_dx": 0.05, "max_dy": 0.05, "max_dpa": 0.01},
    "images": {"border": 2, "rel_rms": 1e-3, "pearson": 0.995, "finite_agree": 0.98},
    "final": {"pearson": 0.99},
}

_REPO_TOLERANCES = Path(__file__).resolve().parents[2] / "tests" / "regression" / "tolerances.toml"


@dataclass
class FileComparison:
    """Outcome for one relative path.

    ``status`` is ``pass``, ``fail``, ``missing`` (one side only) or ``error``
    (unreadable or incompatible files).  ``metrics`` holds the whole-file
    statistics; ``channels`` the per-channel :class:`CompareResult` of a cube.
    """

    relpath: str
    kind: str  # fits | reduclog | hwp_info
    status: str
    rule: str = ""  # which threshold set judged it
    metrics: dict[str, float] = field(default_factory=dict)
    channels: list[CompareResult] = field(default_factory=list)
    message: str = ""

    @property
    def passed(self) -> bool:
        return self.status == "pass"

    @property
    def failed(self) -> bool:
        return self.status in ("fail", "error")

    @property
    def worst_channel(self) -> CompareResult | None:
        if not self.channels:
            return None
        return max(self.channels, key=lambda c: (c.rel_rms, -c.pearson))

    def __str__(self) -> str:
        if self.message:
            return f"{self.relpath}: {self.status} ({self.message})"
        m = ", ".join(f"{k} {v:.4g}" for k, v in self.metrics.items())
        return f"{self.relpath}: {self.status} [{self.rule}] {m}"


@dataclass
class TreeComparison:
    a_dir: Path
    b_dir: Path
    mode: str
    files: list[FileComparison]

    @property
    def n_pass(self) -> int:
        return sum(f.passed for f in self.files)

    @property
    def n_fail(self) -> int:
        return sum(f.failed for f in self.files)

    @property
    def n_missing(self) -> int:
        return sum(f.status == "missing" for f in self.files)

    @property
    def ok(self) -> bool:
        return self.n_fail == 0

    def summary(self) -> str:
        return (
            f"{len(self.files)} file(s): {self.n_pass} pass, {self.n_fail} fail, "
            f"{self.n_missing} missing ({self.mode} mode)"
        )


# --------------------------------------------------------------------------- inputs


def load_tolerances(path: str | Path | None = None) -> dict:
    """``tolerances.toml`` as a dict; the repository copy when ``path`` is None and present."""
    p = Path(path) if path is not None else _REPO_TOLERANCES
    if path is None and not p.exists():
        return default_tolerances()
    with open(p, "rb") as f:
        return tomllib.load(f)


def default_tolerances() -> dict:
    return {k: dict(v) for k, v in DEFAULT_TOLERANCES.items()}


def load_allow_list(path: str | Path | None) -> list[str]:
    """``fnmatch`` patterns, one per line; blank lines and ``#`` comments ignored."""
    if path is None:
        return []
    out = []
    for line in Path(path).read_text().splitlines():
        s = line.split("#", 1)[0].strip()
        if s:
            out.append(s)
    return out


def _allowed(relpath: str, allow_list: list[str]) -> bool:
    name = Path(relpath).name
    return any(fnmatch.fnmatch(relpath, pat) or fnmatch.fnmatch(name, pat) for pat in allow_list)


def _kind(relpath: str) -> str:
    name = Path(relpath).name
    if name == "hwp_info.txt":
        return "hwp_info"
    if name.endswith(".log"):
        return "reduclog"
    return "fits"


def _is_final(relpath: str) -> bool:
    p = Path(relpath)
    return p.parent.name == "proc" and p.name.endswith("_collapsed.fits")


def pair_files(a_dir: str | Path, b_dir: str | Path) -> list[str]:
    """Sorted relative paths matching :data:`PATTERNS` in either tree."""
    found: set[str] = set()
    for root in (Path(a_dir), Path(b_dir)):
        for pat in PATTERNS:
            found.update(p.relative_to(root).as_posix() for p in root.glob(pat) if p.is_file())

    def key(rel: str):
        name = Path(rel).name
        for i, pat in enumerate(PATTERNS):
            if fnmatch.fnmatch(rel, pat):
                return (i, _natural(name))
        return (len(PATTERNS), _natural(name))  # pragma: no cover

    return sorted(found, key=key)


def _natural(s: str) -> tuple:
    return tuple(int(t) if t.isdigit() else t for t in re.split(r"(\d+)", s))


# --------------------------------------------------------------------------- hwp_info


def read_hwp_info(path: str | Path) -> dict[int, tuple[int, float, int]]:
    """``hwp_info.txt`` rows keyed by frame number: ``(hwp_pos, hwp_ang, hwp_cycle)``.

    Accepts the Python layout (``# name hwp_pos hwp_ang hwp_cycle``, names like
    ``n0012leftreg_cal.fits``) and the IDL one (``FILE, HWP_POS, HWP_ANG, CYCLE,``
    header line, then the bare frame number).
    """
    rows: dict[int, tuple[int, float, int]] = {}
    for line in Path(path).read_text().splitlines():
        s = line.strip().replace(",", " ")
        if not s or s.startswith("#"):
            continue
        parts = s.split()
        if len(parts) < 4:
            raise ValueError(f"{path}: malformed line {line!r}")
        if parts[0].upper() in ("FILE", "NAME"):
            continue  # header
        m = re.search(r"(\d+)", parts[0])
        if not m:
            raise ValueError(f"{path}: no frame number in {parts[0]!r}")
        rows[int(m.group(1))] = (int(float(parts[1])), float(parts[2]), int(float(parts[3])))
    return rows


def compare_hwp_info(a: str | Path, b: str | Path, *, ang_tol: float = 1e-3) -> dict:
    """Row-by-row comparison of position, angle and cycle for the common frames.

    Returns ``n`` (common frames), ``n_only_a``, ``n_only_b``, ``n_pos_diff``,
    ``n_cycle_diff``, ``max_dang`` and the list ``diffs`` of
    ``(frame, (pos, ang, cyc)_a, (pos, ang, cyc)_b)`` for rows that differ.
    """
    ra, rb = read_hwp_info(a), read_hwp_info(b)
    common = sorted(set(ra) & set(rb))
    diffs = []
    npos = ncyc = 0
    dang = 0.0
    for n in common:
        pa, aa, ca = ra[n]
        pb, ab, cb = rb[n]
        d = abs(aa - ab)
        dang = max(dang, d)
        bad = pa != pb or ca != cb or d > ang_tol
        npos += pa != pb
        ncyc += ca != cb
        if bad:
            diffs.append((n, ra[n], rb[n]))
    return {
        "n": len(common),
        "n_only_a": len(set(ra) - set(rb)),
        "n_only_b": len(set(rb) - set(ra)),
        "n_pos_diff": npos,
        "n_cycle_diff": ncyc,
        "max_dang": dang,
        "diffs": diffs,
    }


# --------------------------------------------------------------------------- per file


def _compare_fits_file(
    relpath: str, fa: Path, fb: Path, thresholds: dict, rule: str, *, border: int
) -> FileComparison:
    try:
        a, b = read_image(fa), read_image(fb)
    except Exception as exc:  # reported, not raised
        return FileComparison(relpath, "fits", "error", rule, message=f"cannot read: {exc}")
    if a.shape != b.shape:
        return FileComparison(
            relpath, "fits", "error", rule, message=f"shapes differ {a.shape} vs {b.shape}"
        )
    whole = compare_arrays(a, b, border=border, name=relpath)
    channels: list[CompareResult] = []
    if a.ndim == 3:
        channels = [
            compare_arrays(a[i], b[i], border=border, name=f"{relpath}[{i}]")
            for i in range(a.shape[0])
        ]
    metrics = {
        "rel_rms": whole.rel_rms,
        "pearson": whole.pearson,
        "finite_agree": whole.finite_agree,
        "max_abs": whole.max_abs,
        "n": float(whole.n),
    }
    ok = whole.ok(**thresholds)
    return FileComparison(relpath, "fits", "pass" if ok else "fail", rule, metrics, channels)


def _compare_log_file(relpath: str, fa: Path, fb: Path, tol: dict, rule: str) -> FileComparison:
    try:
        d = compare_reduclog(fa, fb)
    except Exception as exc:
        return FileComparison(relpath, "reduclog", "error", rule, message=str(exc))
    ok = all(d[k] <= tol[k] for k in ("max_dx", "max_dy", "max_dpa") if k in tol)
    return FileComparison(relpath, "reduclog", "pass" if ok else "fail", rule, d)


def _compare_hwp_file(relpath: str, fa: Path, fb: Path, rule: str) -> FileComparison:
    try:
        d = compare_hwp_info(fa, fb)
    except Exception as exc:
        return FileComparison(relpath, "hwp_info", "error", rule, message=str(exc))
    diffs = d.pop("diffs")
    ok = not diffs and d["n_only_a"] == 0 and d["n_only_b"] == 0 and d["n"] > 0
    msg = ""
    if diffs:
        shown = ", ".join(f"frame {n}: {x} vs {y}" for n, x, y in diffs[:5])
        more = f" (+{len(diffs) - 5} more)" if len(diffs) > 5 else ""
        msg = f"{len(diffs)} row(s) differ: {shown}{more}"
    elif d["n"] == 0:
        msg = "no common frames"
    return FileComparison(
        relpath,
        "hwp_info",
        "pass" if ok else "fail",
        rule,
        {k: float(v) for k, v in d.items()},
        [],
        msg,
    )


# --------------------------------------------------------------------------- driver


def compare_trees(
    a_dir: str | Path,
    b_dir: str | Path,
    mode: str = "identical",
    tolerances: dict | str | Path | None = None,
    allow_list: list[str] | str | Path | None = None,
) -> TreeComparison:
    """Compare every paired file of two reduction trees.

    ``tolerances`` is the ``tolerances.toml`` dict or path (default: the
    repository copy); ``allow_list`` a list of ``fnmatch`` patterns or the path
    of a file holding them.  Only ``identical`` mode consults the allow-list.
    """
    if mode not in MODES:
        raise ValueError(f"mode must be one of {MODES}, got {mode!r}")
    a_dir, b_dir = Path(a_dir), Path(b_dir)
    for d in (a_dir, b_dir):
        if not d.is_dir():
            raise FileNotFoundError(f"{d} is not a directory")
    tol = tolerances if isinstance(tolerances, dict) else load_tolerances(tolerances)
    allow = allow_list if isinstance(allow_list, list) else load_allow_list(allow_list)
    img, final, logtol = tol.get("images", {}), tol.get("final", {}), tol.get("reduclog", {})
    border = int(img.get("border", 2))
    port_img = {k: img[k] for k in ("rel_rms", "pearson", "finite_agree") if k in img}
    # [final]: speckle-dominated; only the keys it lists are checked.
    port_final = {k: final[k] for k in ("rel_rms", "pearson", "finite_agree") if k in final}
    port_final = {"rel_rms": np.inf, "pearson": -1.0, "finite_agree": 0.0} | port_final
    strict = {"rel_rms": IDENTICAL["rel_rms"], "pearson": -1.0, "finite_agree": 1.0}

    results: list[FileComparison] = []
    for rel in pair_files(a_dir, b_dir):
        fa, fb = a_dir / rel, b_dir / rel
        kind = _kind(rel)
        if not (fa.exists() and fb.exists()):
            only = "only in A" if fa.exists() else "only in B"
            results.append(FileComparison(rel, kind, "missing", message=only))
            continue
        use_port = mode == "port" or _allowed(rel, allow)
        if kind == "fits":
            if use_port:
                rule, thr = ("port/final", port_final) if _is_final(rel) else ("port", port_img)
            else:
                rule, thr = "identical", strict
            results.append(_compare_fits_file(rel, fa, fb, thr, rule, border=border))
        elif kind == "reduclog":
            if use_port:
                rule, thr = "port", logtol
            else:
                rule, thr = "identical", {"max_dx": 1e-6, "max_dy": 1e-6, "max_dpa": 1e-6}
            results.append(_compare_log_file(rel, fa, fb, thr, rule))
        else:
            results.append(_compare_hwp_file(rel, fa, fb, "exact"))
    return TreeComparison(a_dir, b_dir, mode, results)


# --------------------------------------------------------------------------- reports

_CSV_COLUMNS = (
    "relpath",
    "kind",
    "channel",
    "status",
    "rule",
    "rel_rms",
    "pearson",
    "finite_agree",
    "max_abs",
    "n",
    "max_dx",
    "max_dy",
    "max_dha",
    "max_dpa",
    "n_only_a",
    "n_only_b",
    "n_pos_diff",
    "n_cycle_diff",
    "max_dang",
    "message",
)


def _fmt(v) -> str:
    if v is None or v == "":
        return ""
    if isinstance(v, float):
        return f"{v:.6g}"
    return str(v)


def write_csv(result: TreeComparison, path: str | Path) -> Path:
    """One row per file (``channel = all``) plus one per channel of each cube."""
    path = Path(path)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=_CSV_COLUMNS)
        w.writeheader()
        for fc in result.files:
            row = {"relpath": fc.relpath, "kind": fc.kind, "channel": "all", "status": fc.status}
            row |= {"rule": fc.rule, "message": fc.message}
            row |= {k: _fmt(v) for k, v in fc.metrics.items() if k in _CSV_COLUMNS}
            w.writerow(row)
            for i, c in enumerate(fc.channels):
                w.writerow(
                    {
                        "relpath": fc.relpath,
                        "kind": fc.kind,
                        "channel": i,
                        "status": "",
                        "rule": fc.rule,
                        "rel_rms": _fmt(c.rel_rms),
                        "pearson": _fmt(c.pearson),
                        "finite_agree": _fmt(c.finite_agree),
                        "max_abs": _fmt(c.max_abs),
                        "n": c.n,
                    }
                )
    return path


def table_lines(result: TreeComparison) -> list[str]:
    """The summary table as plain-text lines (also the body of ``compare.md``)."""
    hdr = f"{'file':40s} {'status':8s} {'rule':10s} {'rel_rms':>10s} {'pearson':>9s} {'finite':>7s} {'max|diff|':>10s}  note"
    lines = [hdr, "-" * len(hdr)]
    for fc in result.files:
        m = fc.metrics
        if fc.kind == "fits" and m:
            cols = f"{m['rel_rms']:10.3e} {m['pearson']:9.5f} {m['finite_agree']:7.4f} {m['max_abs']:10.3e}"
            note = ""
            wc = fc.worst_channel
            if wc is not None:
                ch = wc.name.rsplit("[", 1)[-1].rstrip("]")
                note = f"worst ch {ch}: rel_rms {wc.rel_rms:.2e}, r {wc.pearson:.4f}"
        elif fc.kind == "reduclog" and m:
            cols = f"{'':10s} {'':9s} {'':7s} {'':10s}"
            note = f"n {int(m['n'])}, dx {m['max_dx']:.3g}, dy {m['max_dy']:.3g}, dpa {m['max_dpa']:.3g}"
        elif fc.kind == "hwp_info" and m:
            cols = f"{'':10s} {'':9s} {'':7s} {'':10s}"
            note = fc.message or (f"n {int(m['n'])}, max dang {m['max_dang']:.3g}")
        else:
            cols = f"{'':10s} {'':9s} {'':7s} {'':10s}"
            note = fc.message
        lines.append(f"{fc.relpath:40s} {fc.status:8s} {fc.rule:10s} {cols}  {note}")
    lines += ["-" * len(hdr), result.summary()]
    return lines


def write_markdown(result: TreeComparison, path: str | Path) -> Path:
    path = Path(path)
    rows = [
        f"# Tree comparison ({result.mode} mode)",
        "",
        f"* A: `{result.a_dir}`",
        f"* B: `{result.b_dir}`",
        f"* {result.summary()}",
        "",
        "| file | status | rule | rel_rms | pearson | finite agree | max abs diff | note |",
        "|---|---|---|---:|---:|---:|---:|---|",
    ]
    for fc in result.files:
        m = fc.metrics
        if fc.kind == "fits" and m:
            cells = [
                f"{m['rel_rms']:.3e}",
                f"{m['pearson']:.5f}",
                f"{m['finite_agree']:.4f}",
                f"{m['max_abs']:.3e}",
            ]
            wc = fc.worst_channel
            note = ""
            if wc is not None:
                ch = wc.name.rsplit("[", 1)[-1].rstrip("]")
                note = f"worst channel {ch}: rel_rms {wc.rel_rms:.2e}, r {wc.pearson:.4f}"
        elif fc.kind == "reduclog" and m:
            cells = ["", "", "", ""]
            note = f"n {int(m['n'])}, max dx {m['max_dx']:.3g}, dy {m['max_dy']:.3g}, dpa {m['max_dpa']:.3g}"
        elif fc.kind == "hwp_info" and m:
            cells = ["", "", "", ""]
            note = fc.message or f"n {int(m['n'])}, max dang {m['max_dang']:.3g}"
        else:
            cells = ["", "", "", ""]
            note = fc.message
        note = note.replace("|", "\\|")
        rows.append(
            f"| `{fc.relpath}` | {fc.status} | {fc.rule} | " + " | ".join(cells) + f" | {note} |"
        )
    path.write_text("\n".join(rows) + "\n")
    return path
