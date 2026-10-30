from __future__ import annotations

import numpy as np
import pytest
from charis_dpp.config.dataset import DatasetConfig, canonical_band
from charis_dpp.config.instrument import load_instrument, to_mjd


def test_epochs_in_order_and_inherit():
    t = load_instrument()
    mjds = [e.valid_from_mjd for e in t.epochs]
    assert mjds == sorted(mjds)
    early = t.for_mjd(57900)
    assert early.epoch_label == "commissioning"
    assert early.pixscale_arcsec == 0.0162 and early.dtel_m == 7.9 and early.north_offset_deg == 2.2
    assert dict(early.sat_spots.contrast_mag)["25"] == pytest.approx(7.5037)
    after = t.for_mjd(57995.0)
    assert after.epoch_label == "satspot-recal-2017"
    assert dict(after.sat_spots.contrast_mag)["25"] == pytest.approx(6.412)
    assert after.pixscale_arcsec == 0.0162  # unchanged by that epoch
    chen = t.for_mjd(to_mjd("2024-10-18"))
    assert chen.pixscale_arcsec == 0.01615 and chen.north_offset_deg == 2.03 and chen.dtel_m == 7.92
    nbs = t.for_mjd(60933.0)
    assert nbs.parity == -1 and nbs.get("pdi.mueller.parity_applied_upstream") is True
    assert t.for_mjd(60932.9).parity == 1
    recoat = t.for_mjd(to_mjd("2026-04-16"))
    assert recoat.get("pdi.mueller.m3.m1") == pytest.approx(1.768)
    assert t.for_mjd(to_mjd("2026-04-15")).get("pdi.mueller.m3.m1") == pytest.approx(2.104)


def test_overrides_and_header_cards():
    inst = load_instrument().for_mjd(58000, {"pixscale_arcsec": 0.0164})
    assert inst.pixscale_arcsec == 0.0164
    cards = dict((k, v) for k, v, _ in inst.to_header_cards())
    assert (
        cards["PIXSCALE"] == 0.0164
        and cards["INSTEPOC"] == "satspot-recal-2017"
        and cards["INSTPAR"] == 1
    )


def test_before_first_epoch_raises():
    with pytest.raises(ValueError):
        load_instrument().for_mjd(50000)


def test_spot_guesses_scale_and_parity():
    t = load_instrument()
    g = t.for_mjd(58000).spot_guesses_px("lowres")
    assert g.shape == (4, 2)
    g2 = t.for_mjd(58000).spot_guesses_px("lowres", astrogrid_lambda_d=11.2)
    assert np.allclose(g2, g * 11.2 / 15.9)
    gp = t.for_mjd(61000).spot_guesses_px("lowres")
    assert np.allclose(gp[:, 0], -g[:, 0]) and np.allclose(gp[:, 1], g[:, 1])


def test_to_mjd_formats():
    assert to_mjd("20170814") == to_mjd("2017-08-14")
    assert to_mjd(57979.5) == 57979.5


@pytest.mark.parametrize(
    "name,expect",
    [
        ("Broadband", "lowres"),
        ("low", "lowres"),
        ("j", "J"),
        ("Ks", "K"),
        ("H", "H"),
        ("bogus", None),
    ],
)
def test_canonical_band(name, expect):
    assert canonical_band(name) == expect


def test_info_roundtrip(hr8799_info, tmp_path):
    cfg = DatasetConfig.from_info(hr8799_info)
    assert cfg.target["name"] == "HR 8799"
    assert cfg.frames("science") == list(range(1, 101))
    assert cfg.frames("sky") == [101, 102, 103]
    assert cfg.reduction["loci"]["nref"] == 50 and cfg.reduction["klip"]["meansub"] is True
    p = cfg.save(tmp_path / "HR8799_lowres.toml")
    back = DatasetConfig.load(p)
    assert back.target == cfg.target
    assert back.frames("science") == cfg.frames("science")
    info = (tmp_path / "HR8799_lowres.info").read_text()
    assert "FNUM_SAT= '1-100'" in info and "KLIPNZON=" in info and "ALOCINRE=" in info
    again = DatasetConfig.from_info(tmp_path / "HR8799_lowres.info")
    assert again.reduction["loci"]["svd_log10"] == cfg.reduction["loci"]["svd_log10"]


def test_resolve_channel_sets(hr8799_info):
    r = DatasetConfig.from_info(hr8799_info).resolve()
    assert r.band == "lowres" and r.nlam == 22
    assert r.reference_channel == 10
    assert r.telluric_channels.tolist() == [5, 6, 14, 15]
    assert r.sub_band_channels["J"].tolist() == [0, 1, 2, 3, 4]
    assert r.sub_band_channels["H"].tolist() == list(range(7, 14))
    assert r.sub_band_channels["K"].tolist() == list(range(16, 21))
    assert r.good_channels().size == 18
    assert r.instrument.epoch_label == "commissioning"  # 2017-08-14 < 2017-08-30


def test_unknown_section_rejected(tmp_path):
    p = tmp_path / "x.toml"
    p.write_text('[observation]\nband="lowres"\ndate="2020-01-01"\n[bogus]\na=1\n')
    with pytest.raises(ValueError):
        DatasetConfig.load(p)
