from __future__ import annotations

import numpy as np
import pytest
from astropy.io.fits import Header
from charis_dpp.io.filelists import (
    format_range,
    frame_filename,
    frame_number,
    frame_paths,
    parse_range,
)
from charis_dpp.io.fits import Cube, get_sat_positions, read_cube, set_sat_positions, write_cube
from charis_dpp.io.obslog import ObsEntry, ObsLog
from charis_dpp.io.params import ParamFile
from charis_dpp.io.reduclog import ReducLog
from charis_dpp.io.wavelengths import (
    band_channels,
    reference_channel,
    telluric_channels,
    wavelengths_um,
)


@pytest.mark.parametrize(
    "spec,expect",
    [
        ("1-5,8,10-12", [1, 2, 3, 4, 5, 8, 10, 11, 12]),
        ("7", [7]),
        ("", []),
        (None, []),
        ("3,1,2,2", [1, 2, 3]),
        (" 1 - 3 ", [1, 2, 3]),
    ],
)
def test_parse_range(spec, expect):
    assert parse_range(spec) == expect


def test_format_range_roundtrip(rng):
    for _ in range(50):
        nums = sorted(set(rng.integers(1, 60, size=rng.integers(1, 25)).tolist()))
        assert parse_range(format_range(nums)) == nums
    assert format_range([1, 2, 3, 5, 7, 8]) == "1-3,5,7-8"
    assert format_range([4]) == "4"
    assert format_range([]) == ""


def test_frame_names():
    assert frame_filename(4, "reg_cal") == "n0004reg_cal.fits"
    assert frame_number("n0004reg_cal.fits") == 4
    assert [p.name for p in frame_paths("1-2", "reduc/reg", "reg")] == [
        "n0001reg.fits",
        "n0002reg.fits",
    ]


def test_paramfile_roundtrip(hr8799_info, tmp_path):
    pf = ParamFile.load(hr8799_info)
    assert pf["NAME"] == "HR 8799"
    assert pf["RSAT"] == -1.0 and pf["RMIN"] == 5 and pf["ALOCISVD"] == -6
    assert pf.get("missing", "d") == "d"
    pf.set("klipnzone", 3, "zones")
    pf.set_frames("fnum_sat", [1, 2, 3, 9])
    out = pf.save(tmp_path / "x.info")
    text = out.read_text()
    assert "KLIPNZON=" in text and text.rstrip().endswith("END")
    assert "****    = '******PSF Sub Range'" in text
    back = ParamFile.load(out)
    assert back["klipnzone"] == 3 and back.frames("fnum_sat") == [1, 2, 3, 9]
    assert back.comment("rmin") == "Min radius for PSF sub"
    # loading must never rewrite the file
    mtime = out.stat().st_mtime_ns
    ParamFile.load(out)
    assert out.stat().st_mtime_ns == mtime


def test_cube_roundtrip(tmp_path):
    data = np.arange(2 * 4 * 5, dtype=float).reshape(2, 4, 5)
    h0 = Header({"OBJECT": "x"})
    h1 = Header({"PARANG": 12.5})
    p = write_cube(tmp_path / "c.fits", data, h0, h1)
    c = read_cube(p)
    assert c.data.shape == (2, 4, 5) and c.nlam == 2 and c.shape == (4, 5)
    assert np.allclose(c.data, data)
    assert c["OBJECT"] == "x" and c["PARANG"] == 12.5
    sats = np.array([[[1.5, 2.5]] * 4] * 2)
    set_sat_positions(c.h1, sats)
    assert np.allclose(get_sat_positions(c.h1, 2), sats)
    assert isinstance(c.copy(), Cube)


def test_reduclog_keyed_by_number(tmp_path):
    log = ReducLog()
    log.add(3, 50.1, 50.2, -1, 0.5, 10.0)
    log.add(1, 50.0, 50.0, -1, 0.1, 5.0)
    p = log.save(tmp_path / "reduc.log")
    back = ReducLog.load(p)
    assert back[3].pa == 10.0 and back[1].xc == 50.0
    assert back.pa([3, 1]).tolist() == [10.0, 5.0]
    with pytest.raises(ValueError):
        back.require([1, 2])


def test_obslog(tmp_path):
    log = ObsLog(
        [
            ObsEntry("CRSA1.fits", "n0001e.fits", "Broadband", 29.5, 1.0, "science"),
            ObsEntry("CRSA2.fits", "n0002e.fits", "Broadband", 29.5, 1.0, "sky"),
        ]
    )
    p = log.save(tmp_path / "obslog.txt")
    back = ObsLog.load(p)
    assert back.frames("sky") == [2] and back.frames("science") == [1]
    assert p.read_text().splitlines()[0].split() == [
        "File_in",
        "File_out",
        "Filter",
        "ExpTime",
        "ParAngle",
        "Health_Flag",
    ]


def test_wavelengths():
    for band, n in (("lowres", 22), ("J", 15), ("H", 20), ("K", 17)):
        w = wavelengths_um(band)
        assert w.size == n and 1.1 < w[0] < w[-1] < 2.4
    w = wavelengths_um("lowres")
    assert reference_channel(w) == 10
    assert telluric_channels(w).tolist() == [5, 6, 14, 15]
    assert band_channels(w)["K"].tolist() == [16, 17, 18, 19, 20]
    assert telluric_channels(wavelengths_um("J"), band="J").size == 0
    h = Header({"CAL_BAND": "Broadband"})
    assert wavelengths_um(h).size == 22
