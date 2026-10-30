"""``obslog.txt`` written by ``imprep`` and read back to pick science/sky frames.

Columns: ``file_in file_out filter exptime parang health``.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

DEFAULT_NAME = "obslog.txt"
HEADER = "File_in File_out Filter ExpTime ParAngle Health_Flag"


@dataclass
class ObsEntry:
    file_in: str
    file_out: str
    filter: str
    exptime: float
    parang: float
    health: str

    @property
    def filenum(self) -> int:
        from charis_dpp.io.filelists import frame_number

        return frame_number(self.file_out)


class ObsLog:
    def __init__(self, entries: list[ObsEntry] | None = None, path: Path | None = None):
        self.entries = list(entries or [])
        self.path = path

    @classmethod
    def load(cls, path: str | Path = DEFAULT_NAME) -> ObsLog:
        path = Path(path)
        entries = []
        for line in path.read_text().splitlines():
            s = line.strip()
            if not s or s.startswith("File_in") or s.startswith("#"):
                continue
            p = s.split()
            if len(p) < 6:
                raise ValueError(f"{path}: malformed line {line!r}")
            entries.append(ObsEntry(p[0], p[1], p[2], float(p[3]), float(p[4]), p[5]))
        return cls(entries, path)

    def save(self, path: str | Path | None = None) -> Path:
        path = Path(path or self.path or DEFAULT_NAME)
        lines = [HEADER] + [
            f"{e.file_in} {e.file_out} {e.filter} {e.exptime:.2f} {e.parang:.2f} {e.health}"
            for e in self.entries
        ]
        path.write_text("\n".join(lines) + "\n")
        self.path = path
        return path

    def frames(self, health: str) -> list[int]:
        return sorted(e.filenum for e in self.entries if e.health == health)
