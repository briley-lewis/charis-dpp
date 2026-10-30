"""Match every frame to a half-wave-plate angle and group frames into HWP cycles
(``charis_pdi_hwp_match.pro``; the header path absorbs ``python/hdr_hwp_match.py``).

Sources: ``'log'`` (a VAMPIRES ``vampHWPLog_*.csv``: angle, UT start, UT end,
target), ``'angles'`` (a ``hwp_angles*.csv`` with ``file`` and ``hwp_ang``
columns) or ``'header'`` (``RET-ANG1``, copied by ``imprep`` from the raw
cube into the prep header).  The result is written as ``hwp_info.txt``
(``name hwp_pos hwp_ang hwp_cycle``) and returned as :class:`HwpTable`.

Fixes vs. IDL: exact and after-the-exposure log times are handled (#222),
the log time format is parsed per row and a UT midnight roll-over is
unwrapped (#223), the exposure reference time is an option (#224), HWP
positions are restricted to the cycle angles (#225).
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from charis_dpp.astro.time_coords import ten
from charis_dpp.config.dataset import DatasetConfig
from charis_dpp.io.filelists import frame_filename
from charis_dpp.io.fits import header_get, read_cube
from charis_dpp.io.tables import read_csv_columns
from charis_dpp.prelims.setupdir import Dirs

log = logging.getLogger(__name__)

__all__ = ["HwpTable", "hwp_match", "assign_cycles", "read_hwp_log"]


@dataclass
class HwpTable:
    names: list[str]
    filenums: list[int]
    hwp_pos: np.ndarray  # index into the cycle angles, -1 if unknown
    hwp_ang: np.ndarray
    hwp_cycle: np.ndarray  # -1 if not part of a complete cycle
    angles_deg: list[float] = field(default_factory=lambda: [0.0, 22.5, 45.0, 67.5])

    @property
    def ncycles(self) -> int:
        c = self.hwp_cycle[self.hwp_cycle >= 0]
        return int(np.unique(c).size)

    def save(self, path: str | Path) -> Path:
        path = Path(path)
        lines = ["# name hwp_pos hwp_ang hwp_cycle"]
        for n, p, a, c in zip(self.names, self.hwp_pos, self.hwp_ang, self.hwp_cycle, strict=True):
            lines.append(f"{n} {int(p)} {float(a):.3f} {int(c)}")
        path.write_text("\n".join(lines) + "\n")
        return path

    @classmethod
    def load(cls, path: str | Path, angles_deg=None) -> HwpTable:
        names, nums, pos, ang, cyc = [], [], [], [], []
        for line in Path(path).read_text().splitlines():
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            n, p, a, c = line.split()[:4]
            names.append(n)
            m = re.search(r"(\d+)", n)
            nums.append(int(m.group(1)) if m else -1)
            pos.append(int(p))
            ang.append(float(a))
            cyc.append(int(c))
        t = cls(names, nums, np.array(pos), np.array(ang), np.array(cyc))
        if angles_deg is not None:
            t.angles_deg = list(angles_deg)
        return t


def _parse_ut(s: str) -> float:
    """UT hours from ``HH:MM:SS(.s)``, ``YYYY-MM-DDTHHMMSS`` or ``...THH:MM:SS`` strings."""
    s = s.strip()
    if "T" in s:
        s = s.split("T", 1)[1]
    if ":" in s:
        return ten(s)
    digits = re.sub(r"\D", "", s)
    return ten(
        f"{digits[0:2]}:{digits[2:4]}:{digits[4:6]}{('.' + digits[6:]) if len(digits) > 6 else ''}"
    )


def read_hwp_log(path: str | Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """VAMPIRES HWP log → ``(angles, ut_start_hours, ut_end_hours)``; end is NaN when absent."""
    angs, t0, t1 = [], [], []
    for line in Path(path).read_text().splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        parts = [p.strip() for p in re.split(r"[,\s]+", s) if p.strip()]
        try:
            a = float(parts[0])
        except ValueError:
            continue  # header line
        angs.append(a)
        t0.append(_parse_ut(parts[1]))
        t1.append(_parse_ut(parts[2]) if len(parts) > 2 and re.search(r"\d", parts[2]) else np.nan)
    if not angs:
        raise ValueError(f"no HWP entries read from {path}")
    return np.array(angs), np.array(t0), np.array(t1)


def _unwrap_hours(t: np.ndarray, ref: float) -> np.ndarray:
    """Hours on the same night as ``ref``: a UT midnight roll-over in either direction
    (log before midnight, exposures after it, or the reverse) is unwrapped (#223)."""
    t = np.asarray(t, dtype=float).copy()
    t[t < ref - 12.0] += 24.0
    t[t > ref + 12.0] -= 24.0
    return t


def assign_cycles(hwp_pos: np.ndarray, times: np.ndarray, nangs: int) -> np.ndarray:
    """Group frames into complete cycles (one frame per HWP position), nearest in time
    to the least-populated position; forward and reverse passes, keep the tighter one."""
    pos = np.asarray(hwp_pos)
    times = np.asarray(times, dtype=float)
    n = pos.size
    counts = np.array([(pos == i).sum() for i in range(nangs)])
    if counts.min() == 0:
        return np.full(n, -1)
    ncyc = int(counts.min())
    anchor = int(np.argmin(counts))
    anchors = np.flatnonzero(pos == anchor)

    def one_pass(order):
        cyc = np.full(n, -1)
        used = np.zeros(n, dtype=bool)
        for i in order:
            a = anchors[i]
            used[a] = True
            cyc[a] = i
            dts = np.abs(times - times[a])
            for j in range(nangs):
                if j == anchor:
                    continue
                cand = np.flatnonzero((pos == j) & ~used)
                if cand.size:
                    k = cand[np.argmin(dts[cand])]
                    used[k] = True
                    cyc[k] = i
        return cyc

    c1 = one_pass(range(ncyc))
    c2 = one_pass(range(ncyc - 1, -1, -1))

    def spread(c):
        return np.mean([np.ptp(times[c == i]) for i in range(ncyc) if (c == i).any()])

    cyc = c1 if spread(c1) <= spread(c2) else c2
    for i in range(ncyc):  # only complete cycles survive
        if (cyc == i).sum() != nangs:
            cyc[cyc == i] = -1
    return cyc


def hwp_match(
    cfg: DatasetConfig | str | Path,
    *,
    source: str | None = None,
    hwp_log: str | Path | None = None,
    time_ref: str | None = None,
    prefix: str = "n",
    suffix: str | None = None,
    outfile: str = "hwp_info.txt",
) -> HwpTable:
    cfg = cfg if isinstance(cfg, DatasetConfig) else DatasetConfig.load(cfg)
    res = cfg.resolve()
    dirs = Dirs(cfg.workdir)
    pd = cfg.pdi
    source = source or pd.get("hwp_source", "log")
    hwp_log = hwp_log or (pd.get("hwp_log") or None)
    time_ref = time_ref or pd.get("hwp_time_ref", "mid")
    if time_ref not in ("start", "mid", "end"):
        raise ValueError("time_ref must be 'start', 'mid' or 'end'")
    angles = [float(a) for a in res.instrument.pdi.hwp_angles_deg]
    nums = cfg.frames("science")
    if suffix is None:
        for cand in ("leftreg_cal", "leftreg", "e_flat_skysub", "e_skysub", "e_flat", "e"):
            d = dirs.reg if "reg" in cand else dirs.prep
            if (d / frame_filename(nums[0], cand, prefix)).exists():
                suffix = cand
                break
        else:
            raise FileNotFoundError("no prepared or registered cubes found for the science frames")
    datadir = dirs.reg if "reg" in suffix else dirs.prep
    names = [frame_filename(n, suffix, prefix) for n in nums]
    if source == "log" and hwp_log is None:
        cands = sorted(dirs.root.glob("vampHWPLog_*.csv")) + sorted(
            dirs.root.glob("hwp_angles*.csv")
        )
        if len(cands) == 1:
            hwp_log = cands[0]
            source = "angles" if "hwp_angles" in cands[0].name else "log"
        elif len(cands) > 1:
            raise ValueError(f"several HWP logs found, give hwp_log=: {[c.name for c in cands]}")
        else:
            source = "header"
            log.info("no HWP log found; using RET-ANG1 from the headers")
    if hwp_log is not None:
        hwp_log = Path(hwp_log) if Path(hwp_log).is_absolute() else dirs.root / hwp_log
        if source == "log" and "hwp_angles" in hwp_log.name:
            source = "angles"

    utc = np.zeros(len(nums))
    t_start = np.zeros(len(nums))
    t_end = np.zeros(len(nums))
    hwp_ang = np.full(len(nums), np.nan)
    orig = []
    for i, nm in enumerate(names):
        c = read_cube(datadir / nm)
        h = c.h0
        t = _parse_ut(str(header_get(h, "UTC-TIME", "00:00:00")))
        exptime = float(header_get(h, "EXPTIME", 0.0))
        mjd0, mjd1 = header_get(h, "MJD-STR"), header_get(h, "MJD-END")
        dit = (
            (float(mjd1) - float(mjd0)) * 86400.0
            if mjd0 is not None and mjd1 is not None
            else exptime
        )
        # UTC-TIME is the end of the read-out in CHARIS cubes: 'start' = t − dit, 'mid' = t − dit/2
        t_end[i] = t
        t_start[i] = t - dit / 3600.0
        if time_ref == "start":
            t -= dit / 3600.0
        elif time_ref == "mid":
            t -= 0.5 * dit / 3600.0
        utc[i] = t
        orig.append(str(header_get(h, "ORIGNAME", "")).strip())
        if source == "header":
            v = header_get(h, "RET-ANG1", None)
            hwp_ang[i] = np.nan if v is None or float(v) < 0 else float(v)
    ref_t = utc[0]
    utc = _unwrap_hours(utc, ref_t)
    t_start = _unwrap_hours(t_start, ref_t)
    t_end = _unwrap_hours(t_end, ref_t)
    if source == "log":
        angs, t0, _ = read_hwp_log(hwp_log)
        t0 = _unwrap_hours(t0, ref_t)
        order = np.argsort(t0)
        angs, t0 = angs[order], t0[order]
        span = max(t0.max(), t_end.max()) - min(t0.min(), t_start.min())
        if span > 16.0:
            log.warning(
                "HWP log and exposures span %.1f h: check that the log is from the same night", span
            )
        for i in range(len(nums)):
            k = (
                int(np.searchsorted(t0, utc[i], side="right")) - 1
            )  # last move at or before the exposure
            if k < 0:
                log.warning(
                    "%s: exposure (UT %.4f h) before the first HWP log entry (UT %.4f h); no angle",
                    names[i],
                    utc[i],
                    t0[0],
                )
                continue
            hwp_ang[i] = angs[k]
            inside = (t0 > t_start[i]) & (t0 < t_end[i])
            if inside.any():  # a move inside the exposure: the frame mixes two HWP angles
                log.warning(
                    "%s: HWP moved during the exposure (log entries %s); angle %.1f assumed",
                    names[i],
                    np.flatnonzero(inside).tolist(),
                    hwp_ang[i],
                )
    elif source == "angles":
        cols = read_csv_columns(hwp_log)
        lc = {k.lower(): v for k, v in cols.items()}
        if "file" not in lc or "hwp_ang" not in lc:
            raise ValueError(f"{hwp_log}: needs 'file' and 'hwp_ang' columns")
        table = {
            re.sub(r"(CRSA|_cube|\.fits)", "", str(f)).strip(): float(a)
            for f, a in zip(lc["file"], lc["hwp_ang"], strict=True)
        }
        for i, o in enumerate(orig):
            hwp_ang[i] = table.get(o, table.get(o.lstrip("0"), np.nan))
    elif source != "header":
        raise ValueError("source must be 'log', 'angles' or 'header'")
    hwp_pos = np.full(len(nums), -1)
    ok = np.isfinite(hwp_ang)
    diff = np.abs(hwp_ang[ok, None] - np.array(angles)[None, :])
    hwp_pos[ok] = np.argmin(diff, axis=1)
    hwp_pos[ok & (np.min(diff, axis=1) > 2.0)] = -1  # not one of the cycle angles
    cyc = assign_cycles(np.where(hwp_pos >= 0, hwp_pos, len(angles)), utc, len(angles))
    table = HwpTable(names, list(nums), hwp_pos, hwp_ang, cyc, angles)
    out = table.save(dirs.root / outfile)
    for i, nm in enumerate(names):
        log.info(
            "  %s  UT %8.4f h  angle %6s  pos %2d  cycle %2d",
            nm,
            utc[i],
            "-" if not np.isfinite(hwp_ang[i]) else f"{hwp_ang[i]:.1f}",
            hwp_pos[i],
            cyc[i],
        )
    nbad = int((hwp_pos < 0).sum())
    if nbad:
        log.warning(
            "%d frame(s) have no HWP angle (not one of %s); excluded from cycles", nbad, angles
        )
    log.info(
        "hwp_match (%s): %d frames, %d complete cycles -> %s",
        source,
        len(nums),
        table.ncycles,
        out,
    )
    if table.ncycles == 0:
        log.warning("no complete HWP cycles found")
    return table
