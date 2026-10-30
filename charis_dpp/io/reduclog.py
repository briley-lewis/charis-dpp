"""``reduc.log`` — registration log consumed by every later step.

Columns (whitespace separated, same as IDL): ``filenum xc yc rsat ha pa``.
Entries are keyed by frame number so later steps match PA/rsat by number,
not by row order (BUGS.md #59).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

DEFAULT_NAME = "reduc.log"


@dataclass
class ReducEntry:
    filenum: int
    xc: float
    yc: float
    rsat: float
    ha: float
    pa: float


class ReducLog:
    def __init__(self, entries: dict[int, ReducEntry] | None = None, path: Path | None = None):
        self.entries: dict[int, ReducEntry] = dict(entries or {})
        self.path = path

    @classmethod
    def load(cls, path: str | Path = DEFAULT_NAME) -> ReducLog:
        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(f"{path} not found: run registration first")
        entries: dict[int, ReducEntry] = {}
        for line in path.read_text().splitlines():
            s = line.strip()
            if not s or s.startswith("#") or s[0].isalpha():
                continue
            parts = s.split()
            if len(parts) < 6:
                raise ValueError(f"{path}: malformed line {line!r}")
            n = int(float(parts[0]))
            entries[n] = ReducEntry(n, *(float(p) for p in parts[1:6]))
        return cls(entries, path)

    def save(self, path: str | Path | None = None) -> Path:
        path = Path(path or self.path or DEFAULT_NAME)
        lines = [
            f"{e.filenum:5d} {e.xc:10.4f} {e.yc:10.4f} {e.rsat:8.3f} {e.ha:12.6f} {e.pa:12.6f}"
            for e in sorted(self.entries.values(), key=lambda e: e.filenum)
        ]
        path.write_text("\n".join(lines) + "\n")
        self.path = path
        return path

    def add(self, filenum: int, xc: float, yc: float, rsat: float, ha: float, pa: float) -> None:
        self.entries[int(filenum)] = ReducEntry(int(filenum), xc, yc, rsat, ha, pa)

    def __contains__(self, filenum: int) -> bool:
        return int(filenum) in self.entries

    def __getitem__(self, filenum: int) -> ReducEntry:
        try:
            return self.entries[int(filenum)]
        except KeyError:
            raise KeyError(f"frame {filenum} is not in {self.path or DEFAULT_NAME}") from None

    def __len__(self) -> int:
        return len(self.entries)

    def require(self, filenums: list[int]) -> None:
        missing = [n for n in filenums if n not in self.entries]
        if missing:
            raise ValueError(
                f"frames {missing} are in the science frame list but not in {self.path or DEFAULT_NAME}; "
                "re-run registration or fix observation.frames.science"
            )

    def pa(self, filenums: list[int]) -> np.ndarray:
        self.require(filenums)
        return np.array([self.entries[n].pa for n in filenums], dtype=float)

    def ha(self, filenums: list[int]) -> np.ndarray:
        self.require(filenums)
        return np.array([self.entries[n].ha for n in filenums], dtype=float)

    def rsat(self, filenums: list[int]) -> np.ndarray:
        self.require(filenums)
        return np.array([self.entries[n].rsat for n in filenums], dtype=float)

    def centers(self, filenums: list[int]) -> np.ndarray:
        self.require(filenums)
        return np.array([(self.entries[n].xc, self.entries[n].yc) for n in filenums], dtype=float)
