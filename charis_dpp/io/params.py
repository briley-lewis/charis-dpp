"""Legacy ``.info`` parameter files (FITS-card style text).

Replaces ``subroutines/param.pro``, ``sxpar_param.pro`` and
``sxaddpar_dl.pro``.  Differences from the IDL implementation:

* reading never rewrites the file (``param.pro`` rewrote it on every ``/get``
  and created an empty one if missing);
* the ``****`` separator cards and comments are preserved verbatim;
* keys are case-insensitive and truncated to 8 characters on write, exactly as
  IDL did, so ``klipnzone`` is stored as ``KLIPNZON``.
"""

from __future__ import annotations

import re
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path

_NUM_RE = re.compile(r"^[-+]?(\d+\.?\d*|\.\d+)([eEdD][-+]?\d+)?$")


@dataclass
class Card:
    key: str
    value: str | float | int | bool | None
    comment: str = ""
    raw: str | None = None  # original line for separators / unknown lines

    def format(self) -> str:
        if self.raw is not None and self.value is None:
            return self.raw
        key = self.key.upper()[:8].ljust(8)
        v = self.value
        if isinstance(v, bool):
            vs = "T" if v else "F"
            vs = vs.rjust(20)
        elif isinstance(v, int):
            vs = f"{v:>20d}"
        elif isinstance(v, float):
            vs = f"{v:>20.10g}" if abs(v) < 1e15 else f"{v:>20.6e}"
        elif v is None:
            vs = ""
        else:
            s = str(v).replace("'", "''")
            vs = f"'{s}'".ljust(20)
        line = f"{key}= {vs}"
        if self.comment:
            line = f"{line} /{self.comment}"
        return line[:80]


def _parse_value(text: str) -> tuple[str | float | int | bool | None, str]:
    text = text.strip()
    if not text:
        return None, ""
    if text.startswith("'"):
        # string: find closing quote (doubled quotes are escapes)
        i = 1
        buf = []
        while i < len(text):
            ch = text[i]
            if ch == "'":
                if i + 1 < len(text) and text[i + 1] == "'":
                    buf.append("'")
                    i += 2
                    continue
                break
            buf.append(ch)
            i += 1
        rest = text[i + 1 :].strip()
        comment = rest[1:] if rest.startswith("/") else ""
        return "".join(buf).rstrip(), comment
    if "/" in text:
        val, comment = text.split("/", 1)
    else:
        val, comment = text, ""
    val = val.strip()
    if val in ("T", "F"):
        return val == "T", comment
    if _NUM_RE.match(val):
        norm = val.replace("D", "e").replace("d", "e")
        if re.match(r"^[-+]?\d+$", norm):
            return int(norm), comment
        return float(norm), comment
    return val, comment


@dataclass
class ParamFile:
    """An ordered collection of cards backed by a ``.info`` file."""

    path: Path | None = None
    cards: list[Card] = field(default_factory=list)

    # --- I/O ---------------------------------------------------------------------
    @classmethod
    def load(cls, path: str | Path) -> ParamFile:
        path = Path(path)
        pf = cls(path=path)
        for line in path.read_text().splitlines():
            if not line.strip():
                continue
            key = line[:8].strip()
            if key == "END":
                continue
            if len(line) > 8 and line[8:10] == "= ":
                value, comment = _parse_value(line[10:])
                pf.cards.append(Card(key.upper(), value, comment.strip()))
            else:
                pf.cards.append(Card(key.upper(), None, "", raw=line))
        return pf

    def save(self, path: str | Path | None = None) -> Path:
        path = Path(path or self.path)
        lines = [c.format() for c in self.cards] + ["END"]
        tmp = path.with_suffix(path.suffix + ".tmp")
        tmp.write_text("\n".join(lines) + "\n")
        tmp.replace(path)
        self.path = path
        return path

    # --- access ------------------------------------------------------------------
    @staticmethod
    def _norm(key: str) -> str:
        return key.upper()[:8]

    def _find(self, key: str) -> Card | None:
        k = self._norm(key)
        for c in self.cards:
            if c.key == k and c.raw is None:
                return c
        return None

    def __contains__(self, key: str) -> bool:
        return self._find(key) is not None

    def get(self, key: str, default=None):
        c = self._find(key)
        return default if c is None else c.value

    def __getitem__(self, key: str):
        c = self._find(key)
        if c is None:
            raise KeyError(key)
        return c.value

    def comment(self, key: str) -> str:
        c = self._find(key)
        return "" if c is None else c.comment

    def set(self, key: str, value, comment: str | None = None) -> None:
        c = self._find(key)
        if c is None:
            self.cards.append(Card(self._norm(key), value, comment or ""))
        else:
            c.value = value
            if comment is not None:
                c.comment = comment

    def set_separator(self, label: str) -> None:
        """A ``****  = '******label'`` card as written by ``charis_newobs``."""
        stars = "*" * (4 + sum(1 for c in self.cards if c.key.startswith("****")))
        self.cards.append(Card(stars[:8], f"******{label}", ""))

    def delete(self, key: str) -> None:
        k = self._norm(key)
        self.cards = [c for c in self.cards if not (c.key == k and c.raw is None)]

    def items(self) -> Iterator[tuple[str, object]]:
        for c in self.cards:
            if c.raw is None and not c.key.startswith("*"):
                yield c.key, c.value

    def to_dict(self) -> dict[str, object]:
        return dict(self.items())

    # --- frame lists -------------------------------------------------------------
    def frames(self, key: str) -> list[int]:
        from charis_dpp.io.filelists import parse_range

        v = self.get(key, "")
        return parse_range(str(v) if v is not None else "")

    def set_frames(self, key: str, numbers, comment: str | None = None) -> None:
        from charis_dpp.io.filelists import format_range

        self.set(key, format_range(numbers), comment)
