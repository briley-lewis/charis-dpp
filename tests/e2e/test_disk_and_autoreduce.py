"""Disk forward modelling, spectral covariance and the automatic chain on the synthetic dataset."""

from __future__ import annotations

import shutil

import numpy as np
import pytest
from astropy.io import fits
from charis_dpp.geometry.grids import radial_grid
from charis_dpp.io.fits import read_cube
from charis_dpp.partialsub.disk_fwdmod import disk_fwdmod
from charis_dpp.partialsub.negspec_disk import negative_disk
from charis_dpp.partialsub.recalc_disk_chisq import recalc_disk_chisq
from charis_dpp.specphotcal.spec_covar import spectral_covariance

pytestmark = pytest.mark.slow


def _ring_model(path, shape, center, r0=18.0, width=3.0, peak=1.0):
    r = radial_grid(shape, center)
    im = peak * np.exp(-0.5 * ((r - r0) / width) ** 2)
    hdr = fits.Header()
    hdr["R0"] = r0
    hdr["ITILT"] = 0.0
    hdr["KSI0"] = width
    fits.HDUList([fits.PrimaryHDU(im.astype(np.float32), hdr)]).writeto(path, overwrite=True)
    return im


def test_disk_fwdmod_attenuates_ring(reduced_dataset):
    d, cfg, truth = reduced_dataset
    c = read_cube(d / "reduc" / "proc" / "aloci.fits")
    center = (c.shape[1] // 2, c.shape[0] // 2)
    mpath = d / "ring_model.fits"
    ring = _ring_model(mpath, c.shape, center)
    res = disk_fwdmod(
        cfg,
        [mpath],
        reducname="aloci.fits",
        psf="empirical",
        absflux=True,
        manscale=True,
        lrmin=8,
        lrmax=35,
        calcattenfact=True,
    )
    assert len(res) == 1
    r = res[0]
    assert r.scale == 1.0 and np.isfinite(r.chisq_binned) and r.nbins > 0
    assert r.params["R0"] == 18.0
    assert r.psfsub_col.exists() and r.input_col.exists() and r.relatten_cube.exists()
    out = read_cube(r.psfsub_col).data
    inp = read_cube(r.input_col).data
    out = out[0] if out.ndim == 3 else out
    inp = inp[0] if inp.ndim == 3 else inp
    rad = radial_grid(c.shape, center)
    on_ring = (np.abs(rad - 18.0) < 2) & np.isfinite(out)
    # an azimuthally symmetric ring is heavily self-subtracted by ADI: processed < input
    assert np.nanmedian(out[on_ring]) < np.nanmedian(inp[on_ring])
    assert np.nanmedian(inp[on_ring]) > 0.3 * ring[on_ring].mean()
    hdr = fits.getheader(r.psfsub_col, 1)
    assert hdr["SFAC"] == 1.0 and "CSQ_COL" in hdr and "SFAC_0" in fits.getheader(r.psfsub_cube, 1)
    assert list((d / "reduc" / "proc" / "model_stat").glob("fit_outcomes_*.dat"))
    # recomputing with a different ROI runs on the stored outputs
    again = recalc_disk_chisq(
        cfg, reducname="aloci.fits", lrmin=10, lrmax=30, roirange=(60, 20, 45)
    )
    assert again and again[0][0] == "ring_model" and np.isfinite(again[0][1])
    # negative injection of the same model re-reduces without error
    p = negative_disk(cfg, mpath, reducname="aloci.fits", outfile="ring_np.fits")
    assert p.exists() and np.isfinite(read_cube(p).data).any()


def test_spectral_covariance_fit(reduced_dataset):
    d, cfg, truth = reduced_dataset
    r = spectral_covariance(cfg, "aloci.fits", rho=5.0, width=1.0, plot=True)
    assert r.psi.shape == (truth["nlam"], truth["nlam"])
    assert np.allclose(np.diag(r.psi), 1.0, atol=1e-6)
    assert 0 <= r.a_rho <= 1 and 0 <= r.a_lam <= 1 and -0.05 <= r.a_del <= 1.0
    assert r.params_file.exists() and (d / "reduc" / "proc" / "spec_covar5lam_D.png").exists()


def test_autoreduce_chain(synthetic_dataset, tmp_path):
    from charis_dpp.cli.commands.autoreduce import autoreduce

    d = tmp_path / "auto"
    shutil.copytree(synthetic_dataset, d)
    for p in d.glob("*.toml"):
        p.unlink()  # autoreduce creates its own dataset config from the raw headers
    from charis_dpp.synthetic import DEFAULT_STAR

    # no network in CI: pass the star (the autoreduce CLI would query SIMBAD and,
    # failing that, leave target.mags empty, which specphot-cal refuses)
    cfg_path = autoreduce(d, stop="detect", starlib="pickles", ntc=1, star=DEFAULT_STAR)
    assert cfg_path.exists()
    assert (d / "reduc" / "proc" / "aloci.fits").exists()
    assert (d / "reduc" / "proc" / "aloci_det.targ").exists()
