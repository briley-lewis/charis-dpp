"""Forward modelling, spectral extraction, detection and contrast on the reduced synthetic dataset."""

from __future__ import annotations

import numpy as np
import pytest
from charis_dpp.io.fits import read_cube
from charis_dpp.partialsub.attenmap import attenuation_map
from charis_dpp.partialsub.fwdmod_planet import fwdmod_planet
from charis_dpp.partialsub.inject import Planet
from charis_dpp.specphotcal.extract_spectrum import extract_spectrum
from charis_dpp.tools.contrast import final_contrast, raw_contrast
from charis_dpp.tools.detection import find_detections

pytestmark = pytest.mark.slow


def test_detection_finds_planet(reduced_dataset):
    d, cfg, truth = reduced_dataset
    dets = find_detections(cfg, "aloci.fits", outfile=d / "reduc" / "proc" / "aloci_det.targ")
    assert dets, "no detections"
    p = truth["planet"]
    best = dets[0]
    assert abs(best.x - p["x_northup"]) < 1.0 and abs(best.y - p["y_northup"]) < 1.0
    assert best.snr > 5
    assert abs(best.log_contrast - p["log_contrast"]) < 0.6  # before throughput correction


def test_fwdmod_throughput_and_spectrum(reduced_dataset):
    d, cfg, truth = reduced_dataset
    p = truth["planet"]
    thr = fwdmod_planet(
        cfg,
        reducname="aloci.fits",
        planets=[Planet(p["x_northup"], p["y_northup"], p["log_contrast"])],
        psf="empirical",
    )
    assert thr.shape == (1, truth["nlam"])
    good = np.isfinite(thr[0])
    assert good.sum() > truth["nlam"] // 2
    assert np.all((thr[0][good] > 0.05) & (thr[0][good] < 1.2)), thr
    assert (d / "reduc" / "proc" / "synth_throughput0.txt").exists()
    # throughput-corrected spectrum recovers the injected contrast per channel
    spec = extract_spectrum(
        cfg, "aloci.fits", (p["x_northup"], p["y_northup"]), throughput="synth_throughput0.txt"
    )
    cal = read_cube(d / "reduc" / "reg" / "n0001reg_cal.fits")
    fstar = np.array([float(cal.h1[f"FSTAR_{l}"]) for l in range(truth["nlam"])])
    recovered = spec.flux / fstar
    ratio = recovered[good] / np.array(p["contrast"])[good]
    assert 0.6 < np.nanmedian(ratio) < 1.4, ratio
    assert spec.snr_collapsed > 5


def test_attenuation_map_and_contrast(reduced_dataset):
    d, cfg, truth = reduced_dataset
    pl, pc = attenuation_map(
        cfg,
        reducname="aloci.fits",
        ntc=1,
        lrmin=10,
        lrmax=30,
        dlrstep=10,
        contrast=-4.5,
        psf="gaussian",
    )
    att = read_cube(pl).data
    assert att.shape[0] == truth["nlam"]
    assert np.nanmedian(att) > 0.05
    c = final_contrast(cfg, "aloci.fits", attenmap=pc, attenmap_lambda=pl, plot=True)
    assert c.broadband.size > 10 and np.nanmin(c.broadband) > 0
    assert set(c.bands) == {"J", "H", "K"}
    r = raw_contrast(cfg)
    assert np.nanmedian(r.broadband) > np.nanmedian(
        c.broadband
    )  # PSF subtraction improves contrast
