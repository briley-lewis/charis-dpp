"""The IDL epoch table in idl/setup/charis_get_constant.pro must agree with instrument.toml."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pytest
from charis_dpp.config.instrument import load_instrument

REPO = Path(__file__).resolve().parents[2]
PRO = REPO / "idl" / "setup" / "charis_get_constant.pro"


def parse_idl_epochs() -> list[dict]:
    epochs = []
    for line in PRO.read_text().splitlines():
        m = re.match(r"\s*epochs\.add,\s*\{(.*)\}\s*$", line)
        if not m:
            continue
        d = {}
        for item in re.split(r",\s*(?=[A-Za-z_0-9]+:)", m.group(1)):
            k, v = item.split(":", 1)
            v = v.strip()
            if v.startswith("'"):
                d[k.strip().lower()] = v.strip("'")
            else:
                d[k.strip().lower()] = float(
                    re.sub(r"[dD]0?$|d(?=[-+]?\d)", "", v.replace("d0", ""))
                )
        epochs.append(d)
    return epochs


def test_idl_epoch_table_matches_instrument_toml():
    idl = parse_idl_epochs()
    assert len(idl) >= 5
    table = load_instrument()
    py_mjds = sorted(e.valid_from_mjd for e in table.epochs)
    assert [e["mjd"] for e in idl] == pytest.approx(py_mjds, abs=1.0)
    merged: dict = {}
    for ep in idl:
        merged.update(ep)
        inst = table.for_mjd(ep["mjd"] + 0.5)
        assert inst.epoch_label == ep["label"]
        assert merged["dtel"] == pytest.approx(inst.dtel_m)
        assert merged["pixscale"] == pytest.approx(inst.pixscale_arcsec)
        assert merged["angoffset"] == pytest.approx(inst.north_offset_deg)
        assert merged["lat"] == pytest.approx(inst.latitude_deg)
        assert merged["lng"] == pytest.approx(inst.longitude_deg)
        assert int(merged["parity"]) == inst.parity
        assert merged["readtime"] == pytest.approx(inst.read_time_s)
        cm = inst.sat_spots.contrast_mag
        assert merged["satcon25"] == pytest.approx(float(cm["25"]))
        assert merged["satcon50"] == pytest.approx(float(cm["50"]))
        assert merged["satconerr"] == pytest.approx(float(inst.sat_spots.contrast_err_mag["25"]))
        m3 = inst.get("pdi.mueller.m3")
        for k in ("m1", "b1", "m2", "b2"):
            assert merged[f"m3_{k}"] == pytest.approx(float(m3[k]))
        rotzp = inst.get("nbs.rotzp_deg", None)
        if rotzp is not None:
            assert merged["rotzp"] == pytest.approx(float(rotzp))
    # the Python table also stores the spot guesses the IDL has as absolute 201-px positions
    text = PRO.read_text()
    for band, key in (("lowres", "lowspots"), ("J", "jspots"), ("H", "hspots"), ("K", "kspots")):
        xs = re.search(rf"'{key}_x': value = \[(.*?)\]", text).group(1)
        ys = re.search(rf"'{key}_y': value = \[(.*?)\]", text).group(1)
        idl_xy = np.column_stack(
            [[float(v) for v in xs.split(",")], [float(v) for v in ys.split(",")]]
        )
        inst = table.for_mjd(58000.0)
        py_xy = inst.spot_guesses_px(band) + inst.detector.cube_size_px // 2
        assert np.allclose(idl_xy, py_xy, atol=0.05), band
