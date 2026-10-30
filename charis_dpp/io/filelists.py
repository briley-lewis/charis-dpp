"""Frame-number ranges and file naming.

Replaces ``nbrlist.pro``, ``filelist.pro``, ``nbr2txt.pro`` and the range
formatter inside ``charis_set_files.pro``.  Frame numbers are integers; file
names are ``<prefix><NNNN><suffix>.fits`` (e.g. ``n0004reg_cal.fits``).
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from pathlib import Path

DEFAULT_PREFIX = "n"
DEFAULT_NDIGITS = 4


def parse_range(spec: str | Sequence[int] | None) -> list[int]:
    """``"1-5,8,10-12"`` → ``[1, 2, 3, 4, 5, 8, 10, 11, 12]`` (sorted, unique).

    An empty or ``None`` spec gives an empty list.  Whitespace is ignored.
    """
    if spec is None:
        return []
    if not isinstance(spec, str):
        return sorted({int(x) for x in spec})
    out: set[int] = set()
    for part in spec.replace(" ", "").split(","):
        if not part:
            continue
        if "-" in part:
            a, b = part.split("-", 1)
            lo, hi = int(a), int(b)
            if hi < lo:
                raise ValueError(f"descending range {part!r} in {spec!r}")
            out.update(range(lo, hi + 1))
        else:
            out.add(int(part))
    return sorted(out)


def format_range(numbers: Iterable[int]) -> str:
    """Inverse of :func:`parse_range`: ``[1,2,3,5,7,8] → "1-3,5,7-8"``.

    Fixes ``charis_set_files.pro`` which crashed on a single frame and assumed
    sorted, unique input.
    """
    nums = sorted({int(n) for n in numbers})
    if not nums:
        return ""
    runs: list[list[int]] = [[nums[0], nums[0]]]
    for n in nums[1:]:
        if n == runs[-1][1] + 1:
            runs[-1][1] = n
        else:
            runs.append([n, n])
    return ",".join(f"{a}-{b}" if a != b else f"{a}" for a, b in runs)


def frame_filename(
    number: int,
    suffix: str = "",
    prefix: str = DEFAULT_PREFIX,
    ndigits: int = DEFAULT_NDIGITS,
    ext: str = ".fits",
) -> str:
    """``frame_filename(4, "reg_cal") == "n0004reg_cal.fits"``."""
    if number < 0:
        raise ValueError("frame number must be non-negative")
    return f"{prefix}{number:0{ndigits}d}{suffix}{ext}"


def frame_paths(
    numbers: Iterable[int] | str,
    directory: str | Path = ".",
    suffix: str = "",
    prefix: str = DEFAULT_PREFIX,
    ndigits: int = DEFAULT_NDIGITS,
    ext: str = ".fits",
) -> list[Path]:
    """Paths for a range spec or list of numbers (no existence check)."""
    nums = parse_range(numbers) if isinstance(numbers, str) else list(numbers)
    d = Path(directory)
    return [d / frame_filename(n, suffix, prefix, ndigits, ext) for n in nums]


def frame_number(path: str | Path, prefix: str = DEFAULT_PREFIX) -> int:
    """Recover the frame number from a pipeline file name (``n0004reg.fits`` → 4)."""
    name = Path(path).name
    if not name.startswith(prefix):
        raise ValueError(f"{name!r} does not start with prefix {prefix!r}")
    digits = ""
    for ch in name[len(prefix) :]:
        if ch.isdigit():
            digits += ch
        else:
            break
    if not digits:
        raise ValueError(f"no frame number in {name!r}")
    return int(digits)


def existing(paths: Iterable[Path]) -> list[Path]:
    return [p for p in paths if p.exists()]
