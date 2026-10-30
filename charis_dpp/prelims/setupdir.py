"""Dataset directory layout (``charis_setupdir.pro``)."""

from __future__ import annotations

from pathlib import Path

SUBDIRS = (
    "data/raw",
    "reduc/prep",
    "reduc/reg",
    "reduc/rsub",
    "reduc/proc",
    "psfmodel",
    "log",
)


class Dirs:
    """Resolved dataset directories."""

    def __init__(self, root: str | Path = "."):
        self.root = Path(root)

    @property
    def raw(self) -> Path:
        return self.root / "data" / "raw"

    @property
    def prep(self) -> Path:
        return self.root / "reduc" / "prep"

    @property
    def reg(self) -> Path:
        return self.root / "reduc" / "reg"

    @property
    def rsub(self) -> Path:
        return self.root / "reduc" / "rsub"

    @property
    def proc(self) -> Path:
        return self.root / "reduc" / "proc"

    @property
    def psfmodel(self) -> Path:
        return self.root / "psfmodel"

    @property
    def ladder(self) -> Path:
        return self.root / "ladder"

    @property
    def log(self) -> Path:
        return self.root / "log"

    @property
    def reduclog(self) -> Path:
        return self.root / "reduc.log"

    @property
    def obslog(self) -> Path:
        return self.root / "obslog.txt"

    def ensure(self) -> Dirs:
        for s in SUBDIRS:
            (self.root / s).mkdir(parents=True, exist_ok=True)
        return self


def setup_dirs(root: str | Path = ".") -> Dirs:
    return Dirs(root).ensure()
