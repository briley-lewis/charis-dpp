"""Guard against the hard-coded constants that plagued the IDL tree."""

from __future__ import annotations

import io
import re
import tokenize
from pathlib import Path

import pytest

PKG = Path(__file__).resolve().parents[2] / "charis_dpp"

FORBIDDEN = [
    (r"(?<![\w.])201(?![\w.])", "frame size 201"),
    (r"0\.0162\b|0\.01615\b|0\.0164\b", "pixel scale literal"),
    (r"\b7\.9[02]?\b", "telescope diameter literal"),
    (r"\b57995\b|\b60933\b", "MJD epoch literal"),
    (r"\b1\.475\b", "read time literal"),
    (r"\[\s*5\s*,\s*6\s*,\s*14\s*,\s*(15|21)\s*\]", "telluric channel list"),
    (r"\bdialog_pickfile\b|\binput\(", "interactive prompt"),
    (r"\[\s*0\s*:\s*21\s*\]", "22-channel slice"),
]
ALLOWLIST = {"data", "synthetic.py", "plotting.py"}


def _files():
    for p in PKG.rglob("*.py"):
        if any(part in ALLOWLIST for part in p.relative_to(PKG).parts):
            continue
        yield p


def _code_lines(path: Path) -> dict[int, str]:
    """Source lines with string literals and comments blanked out."""
    src = path.read_text()
    lines = src.splitlines()
    keep = {i + 1: line for i, line in enumerate(lines)}
    for tok in tokenize.generate_tokens(io.StringIO(src).readline):
        if tok.type in (tokenize.STRING, tokenize.COMMENT):
            (r0, c0), (r1, c1) = tok.start, tok.end
            for r in range(r0, r1 + 1):
                line = keep[r]
                a = c0 if r == r0 else 0
                b = c1 if r == r1 else len(line)
                keep[r] = line[:a] + " " * (b - a) + line[b:]
    return keep


@pytest.mark.parametrize("pattern,label", FORBIDDEN, ids=[f[1] for f in FORBIDDEN])
def test_no_hardcoded_constants(pattern, label):
    rx = re.compile(pattern)
    hits = []
    for p in _files():
        for i, code in _code_lines(p).items():
            if rx.search(code):
                hits.append(f"{p.relative_to(PKG)}:{i}: {code.strip()}")
    assert not hits, f"{label} hard-coded:\n" + "\n".join(hits)
