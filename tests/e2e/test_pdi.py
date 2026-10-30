"""PDI chain on a synthetic two-beam dataset with a tangentially polarised ring."""

from __future__ import annotations

import json

import numpy as np
import pytest
from charis_dpp.geometry.grids import image_center, radial_grid
from charis_dpp.io.fits import read_cube
from charis_dpp.pdi.hwp_match import hwp_match
from charis_dpp.pdi.pdi import pdi
from charis_dpp.pdi.register import pdi_register_cube
from charis_dpp.pdi.specphot_cal import pdi_specphot_cal
from charis_dpp.pdi.subtract_sky import pdi_subtract_sky
from charis_dpp.prelims.imprep import imprep
from charis_dpp.synthetic import SyntheticConfig, make_dataset

pytestmark = pytest.mark.slow


@pytest.fixture(scope="module")
def pdi_dataset(tmp_path_factory):
    d = tmp_path_factory.mktemp("pdi")
    cfg = SyntheticConfig(
        nframes=8,
        shape=(161, 161),  # beams at ±(25, 13) px keep the two beams' sat spots apart
        n_speckles=8,
        nsky=2,
        seed=3,
        astrogrid_lambda_d=11.2,
        pdi=True,
        planet=None,
        pa_step_deg=3.0,
    )
    cfg_path = make_dataset(d, cfg)
    truth = json.loads((d / "truth.json").read_text())
    return d, cfg_path, truth


def test_pdi_chain(pdi_dataset):
    d, cfg, truth = pdi_dataset
    obslog = imprep(cfg, pdi=True)
    assert len(obslog.frames("science")) == 8 and len(obslog.frames("sky")) == 2
    pdi_subtract_sky(cfg)
    rl_l, rl_r = pdi_register_cube(cfg)
    assert len(rl_l) == 8 and len(rl_r) == 8
    off = np.array(truth["pdi"]["polcent_offset_px"])
    c0 = np.array(truth["frames"][0]["center"])
    assert np.allclose(rl_r.centers([1])[0], c0 + off, atol=0.3)
    assert np.allclose(rl_l.centers([1])[0], c0 - off, atol=0.3)
    assert (d / "reduc_left.log").exists() and (d / "reduc.log").exists()
    pdi_specphot_cal(cfg, calmethod=1, starlib="pickles")
    assert (d / "reduc" / "reg" / "n0001leftreg_cal.fits").exists()
    # HWP matching from the VAMPIRES log and from the headers agree
    t_log = hwp_match(cfg, source="log")
    t_hdr = hwp_match(cfg, source="header", outfile="hwp_info_hdr.txt")
    assert t_log.ncycles == 2 and t_hdr.ncycles == 2
    assert np.allclose(t_log.hwp_ang, t_hdr.hwp_ang) and list(t_log.hwp_cycle) == list(
        t_hdr.hwp_cycle
    )
    # ideal synthetic optics (HWP not synchronised with the field rotation): no Mueller
    # corrections; the ring must come out in Qphi, not Uphi
    r = pdi(cfg, mm_ip_corr=False, mm_crosstalk_corr=False, adisync=False, bg_sub=True, plot=False)
    assert r.ncycles == 2
    shape = r.qphi.shape[-2:]
    rad = radial_grid(shape, image_center(shape))
    on = np.abs(rad - truth["pdi"]["ring_r_px"]) < 1.5
    off_ring = (rad > 26) & (rad < 40)
    qphi = read_cube(d / "reduc" / "proc" / "pdi_qphi_collapsed.fits").data
    uphi = read_cube(d / "reduc" / "proc" / "pdi_uphi_collapsed.fits").data
    qphi = qphi[0] if qphi.ndim == 3 else qphi
    uphi = uphi[0] if uphi.ndim == 3 else uphi
    sig = np.nanstd(qphi[off_ring])
    assert np.nanmedian(qphi[on]) > 5 * sig
    assert abs(np.nanmedian(uphi[on])) < 0.15 * np.nanmedian(qphi[on])
    # degree of polarisation on the ring: Qphi / (ring + star halo) ~ pol_fraction * ring/total
    cal = read_cube(d / "reduc" / "reg" / "n0001leftreg_cal.fits")
    assert "FSCALE0" in cal.h1
    # the full Mueller-matrix path runs and still shows the ring
    r2 = pdi(cfg, adisync=False, unres_pol_sub=True, outfile="pdimm", plot=True)
    assert np.isfinite(r2.qphi).any() and r2.unres_q.shape == (truth["nlam"], 2)
    # the real-instrument model rotates the axes of these ideal data, so test the
    # polarised intensity rather than Qphi
    pi2 = np.hypot(np.nanmedian(r2.q, axis=0), np.nanmedian(r2.u, axis=0))
    assert np.nanmedian(pi2[on]) > 3 * np.nanmedian(pi2[off_ring])
    assert (d / "reduc" / "proc" / "pdimm_unres_pol.png").exists()
    assert (d / "reduc" / "proc" / "coverage_map.fits").exists()


def test_pdi_xcorr_registration(pdi_dataset):
    d, cfg, truth = pdi_dataset
    if not (d / "reduc" / "prep" / "n0001e.fits").exists():
        imprep(cfg, pdi=True)
    rl_l, rl_r = pdi_register_cube(cfg, method="xcorr", nomask=True, suffix="e")
    assert len(rl_l) == 8 and len(rl_r) == 8
    # frames are aligned to the first one: relative shifts follow the dithers of the truth
    c = np.array([f["center"] for f in truth["frames"]])
    rel_truth = c[1:] - c[0]
    rel_meas = rl_r.centers(list(range(2, 9))) - rl_r.centers([1])[0]
    assert np.allclose(rel_meas, rel_truth, atol=0.3)
    assert (d / "reduc" / "reg" / "n0001rightreg.fits").exists()
