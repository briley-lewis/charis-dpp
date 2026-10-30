"""End-to-end: synthetic raw cubes -> imprep -> sky subtraction -> registration -> reduc.log."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import numpy as np
import pytest
from charis_dpp.config.dataset import DatasetConfig
from charis_dpp.io.fits import get_psf_center, get_sat_positions, read_cube
from charis_dpp.io.obslog import ObsLog
from charis_dpp.io.reduclog import ReducLog
from charis_dpp.prelims.imprep import imprep
from charis_dpp.reg import check_psf_quality, register_cube
from charis_dpp.skysub import subtract_sky


@pytest.fixture(scope="module")
def workdir(synthetic_dataset: Path, tmp_path_factory) -> Path:
    d = tmp_path_factory.mktemp("work")
    shutil.copytree(synthetic_dataset, d, dirs_exist_ok=True)
    return d


@pytest.fixture(scope="module")
def truth(workdir) -> dict:
    return json.loads((workdir / "truth.json").read_text())


@pytest.fixture(scope="module")
def cfg_path(workdir) -> Path:
    return next(workdir.glob("*.toml"))


def test_imprep(workdir, cfg_path, truth):
    log = imprep(cfg_path)
    assert len(log.entries) == len(truth["frames"]) + len(truth["sky_frames"])
    assert log.frames("science") == [f["filenum"] for f in truth["frames"]]
    assert log.frames("sky") == [f["filenum"] for f in truth["sky_frames"]]
    back = ObsLog.load(workdir / "obslog.txt")
    assert back.frames("sky") == log.frames("sky")
    cfg = DatasetConfig.load(cfg_path)
    assert cfg.frames("science") == log.frames("science")
    c = read_cube(workdir / "reduc" / "prep" / "n0001e.fits")
    for k in (
        "HA",
        "LST",
        "EXP1TIME",
        "COADDS",
        "LAT",
        "LNG",
        "PA",
        "PARITY",
        "INSTEPOC",
        "D_IMRANG",
    ):
        assert k in c.h0, k
    assert c.h0["COADDS"] == truth["config"]["reads"]
    assert abs(c.h0["PA"] - truth["frames"][0]["parang"]) < 1e-6
    assert c.h0["PARITY"] == 0
    assert c.h1["CUNIT3"] == "microns"


def test_subtract_sky(workdir, cfg_path, truth):
    out = subtract_sky(cfg_path)
    assert len(out) == len(truth["frames"])
    c = read_cube(out[0])
    cfgd = truth["config"]
    orig = read_cube(workdir / "reduc" / "prep" / "n0001e.fits")
    # what was removed in the corner is the sky level (the halo stays)
    removed = np.nanmedian(orig.data[:, :10, :10] - c.data[:, :10, :10])
    assert abs(removed - cfgd["sky_level"]) < 1.0
    assert c.h1["SKYSUB"] is True
    # scaled modes run (the synthetic halo dominates the sky, so no level check)
    out2 = subtract_sky(cfg_path, scalesky=True)
    assert "SKYSCALE" in read_cube(out2[0]).h1
    out3 = subtract_sky(cfg_path, scaleindchan=True, skylim=1.0)
    assert "SKYSC_0" in read_cube(out3[0]).h1
    subtract_sky(cfg_path)  # leave plain subtraction for the registration tests


def test_register_sats(workdir, cfg_path, truth):
    log = register_cube(cfg_path)
    assert isinstance(log, ReducLog) and len(log) == len(truth["frames"])
    for fr in truth["frames"]:
        e = log[fr["filenum"]]
        # gcntrd marginal fits are biased by nearby speckles at the 0.1 px level (as in IDL)
        assert abs(e.xc - fr["center"][0]) < 0.2, (e, fr["center"])
        assert abs(e.yc - fr["center"][1]) < 0.2
        assert abs(e.pa - fr["parang"]) < 1e-6
    reg = read_cube(workdir / "reduc" / "reg" / "n0001reg.fits")
    xc, yc = get_psf_center(reg.h1)
    assert (xc, yc) == (50.0, 50.0)
    sats = get_sat_positions(reg.h1, reg.nlam)
    assert np.isfinite(sats).all()
    # after registration the spots should sit at truth offsets about the centre
    fr = truth["frames"][0]
    tr = np.array(fr["spots"]) - np.array(fr["center"])
    assert np.nanmax(np.abs((sats - np.array([50.0, 50.0])) - tr)) < 0.4
    # the star is now at the centre in every channel
    for il in (0, reg.nlam // 2, reg.nlam - 1):
        iy, ix = np.unravel_index(np.nanargmax(reg.data[il][40:61, 40:61]), (21, 21))
        assert abs(ix + 40 - 50) <= 1 and abs(iy + 40 - 50) <= 1
    assert "REFSLICE" in reg.h1 and "SHIFT_0" in reg.h1
    cfg = DatasetConfig.load(cfg_path)
    assert cfg.geometry["star_center"] == [50.0, 50.0]


def test_register_fitgauss_is_less_biased(workdir, cfg_path, truth):
    log = register_cube(cfg_path, fitgauss=True, out_suffix="reggauss")
    errs = [abs(log[fr["filenum"]].xc - fr["center"][0]) for fr in truth["frames"]]
    assert np.median(errs) < 0.12


def test_register_psf_method(workdir, cfg_path, truth):
    log = register_cube(cfg_path, method="psf", out_suffix="regpsf")
    for fr in truth["frames"][:3]:
        e = log[fr["filenum"]]
        assert abs(e.xc - fr["center"][0]) < 0.15 and abs(e.yc - fr["center"][1]) < 0.15


def test_psf_quality(workdir, cfg_path):
    r = check_psf_quality(cfg_path, plot=True)
    assert np.isfinite(r).all() and (r > 0).all()
    assert (workdir / "reduc" / "reg" / "PSFquality.txt").exists()


def test_specphot_cal(workdir, cfg_path, truth):
    from charis_dpp.geometry.apertures import aper_idl
    from charis_dpp.specphotcal import specphot_cal
    from charis_dpp.specphotcal.specphot_cal import read_fscale

    out = specphot_cal(cfg_path, calmethod=0, starlib="pickles")
    assert len(out) == len(truth["frames"])
    c = read_cube(out[0])
    fscale, cerr, fstar, units = read_fscale(c.h1)
    assert units == "mJy" and np.all(np.isfinite(fscale)) and np.all(fscale > 0)
    assert np.all(cerr < 0.2)
    # the calibrated stellar aperture flux must equal FSTAR (spots are exact PSF copies)
    xc, yc = get_psf_center(c.h1)
    for il in (2, 10, 18):
        r = float(c.h1[f"R_AP{il}"])
        flux = aper_idl(c.data[il], xc, yc, r, setskyval=0.0)[0]
        assert abs(flux / fstar[il] - 1) < 0.12, (il, flux, fstar[il])
    # the median-cube method gives the same scale to a few per cent
    out1 = specphot_cal(cfg_path, calmethod=1, starlib="pickles", out_suffix="reg_cal1")
    f1 = read_fscale(read_cube(out1[0]).h1)[0]
    assert np.nanmedian(np.abs(f1 / fscale - 1)) < 0.05


def test_psf_models(workdir, cfg_path, truth):
    from charis_dpp.specphotcal import make_empirical_psf, make_synthetic_psf

    p = make_empirical_psf(cfg_path)
    psf = read_cube(p)
    assert psf.data.shape == (truth["nlam"], 21, 21)
    assert np.allclose(np.nansum(psf.data, axis=(1, 2)), 1.0, atol=1e-6)
    peak = np.array([np.unravel_index(np.nanargmax(s), s.shape) for s in psf.data])
    assert np.all(peak == 10)
    assert "WAVEL0" in psf.h1
    # empirical PSF width tracks the injected FWHM
    from charis_dpp.image.centroid import fit_gaussian2d

    f = fit_gaussian2d(psf.data[10], 10, 10, truth["fwhm_px"][10])
    assert abs(f["fwhm_x"] - truth["fwhm_px"][10] * truth["config"]["fwhm_factor"]) < 0.3
    q = make_synthetic_psf(cfg_path, kind="gaussian")
    assert read_cube(q).data.shape == (truth["nlam"], 21, 21)
