"""End-to-end PSF subtraction on the synthetic dataset: the injected planet must survive ADI."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import numpy as np
import pytest
from charis_dpp.filt import imrsub
from charis_dpp.geometry.apertures import aper_idl
from charis_dpp.image.snr import snr_map
from charis_dpp.io.fits import read_cube
from charis_dpp.prelims.imprep import imprep
from charis_dpp.psfsub import adialoci, adiklip
from charis_dpp.reg import register_cube
from charis_dpp.skysub import subtract_sky
from charis_dpp.specphotcal import specphot_cal

pytestmark = pytest.mark.slow


@pytest.fixture(scope="module")
def reduced(synthetic_dataset: Path, tmp_path_factory):
    d = tmp_path_factory.mktemp("psfsub")
    shutil.copytree(synthetic_dataset, d, dirs_exist_ok=True)
    cfg = next(d.glob("*.toml"))
    imprep(cfg)
    subtract_sky(cfg)
    register_cube(cfg)
    specphot_cal(cfg, calmethod=1, starlib="pickles")
    imrsub(cfg, prad=True)
    truth = json.loads((d / "truth.json").read_text())
    return d, cfg, truth


def _planet_snr(path: Path, truth: dict, fwhm: float) -> tuple[float, tuple[float, float]]:
    col = read_cube(path).data
    p = truth["planet"]
    r = snr_map(col, fwhm, 45, coord=(p["x_northup"], p["y_northup"]))
    return r.snr, r.position


def test_classical_adi_products(reduced):
    d, cfg, truth = reduced
    col = read_cube(d / "reduc" / "proc" / "imcubefinal_derot_collapsed.fits").data
    assert col.shape == (101, 101)
    assert (d / "reduc" / "rsub" / "n0001rsub.fits").exists()


def test_adialoci_recovers_planet(reduced):
    d, cfg, truth = reduced
    final, col = adialoci(
        cfg,
        rsubval=1,
        nfwhm=0.75,
        drsub=5,
        na=100,
        svd=1e-6,
        nref=0,
        rmin=6,
        rmax=40,
        outfile="aloci.fits",
    )
    assert final.exists() and col.exists()
    c = read_cube(final)
    assert c.data.shape == (truth["nlam"], 101, 101)
    for k in ("LOCI_NFW", "LOCI_NA", "SVD", "ROTANG", "PSFALGO"):
        assert k in c.h1, k
    snr, pos = _planet_snr(col, truth, truth["fwhm_px"][10])
    p = truth["planet"]
    assert snr > 5, snr
    assert abs(pos[0] - p["x_northup"]) < 1.0 and abs(pos[1] - p["y_northup"]) < 1.0
    # residuals inside the subtraction region should be much smaller than the raw halo
    raw = read_cube(d / "reduc" / "proc" / "imcubefinal_derot_collapsed.fits").data
    y, x = np.mgrid[0:101, 0:101]
    rr = np.hypot(x - 50, y - 50)
    ring = (rr > 10) & (rr < 30)
    assert np.nanstd(read_cube(col).data[ring]) < 0.5 * np.nanstd(raw[ring])


def test_adiklip_recovers_planet(reduced):
    d, cfg, truth = reduced
    final, col = adiklip(
        cfg, rsubval=1, npca=4, nzone=2, drsub=5, rmin=6, rmax=40, outfile="klip.fits"
    )
    c = read_cube(final)
    assert "KLIP_PCA" in c.h1 and c.h1["KLIP_PCA"] == 4
    snr, pos = _planet_snr(col, truth, truth["fwhm_px"][10])
    p = truth["planet"]
    assert snr > 5, snr
    assert abs(pos[0] - p["x_northup"]) < 1.0 and abs(pos[1] - p["y_northup"]) < 1.0


def test_sdialoci_and_ssdi_run(reduced):
    from charis_dpp.psfsub.sdialoci import sdialoci
    from charis_dpp.psfsub.ssdi import ssdi

    d, cfg, truth = reduced
    final, col = sdialoci(
        cfg, rsubval=1, drsub=5, na=50, rmin=6, rmax=40, nfwhm=0.75, outfile="sdi.fits"
    )
    assert read_cube(final).data.shape[0] == truth["nlam"]
    out = ssdi(cfg, source="aloci")
    assert out.exists()


def test_planet_flux_after_fwdmod_is_consistent(reduced):
    """Aperture flux of the planet in the A-LOCI result is a fraction (throughput < 1) of the injected flux."""
    d, cfg, truth = reduced
    c = read_cube(d / "reduc" / "proc" / "aloci.fits")
    cal = read_cube(d / "reduc" / "reg" / "n0001reg_cal.fits")
    p = truth["planet"]
    il = 10
    r_ap = float(cal.h1[f"R_AP{il}"])
    flux = aper_idl(c.data[il], p["x_northup"], p["y_northup"], r_ap, setskyval=0.0)[0]
    fstar = float(cal.h1[f"FSTAR_{il}"])
    injected = fstar * p["contrast"][il]  # planet flux in the same (aperture-calibrated) units
    assert 0.2 < flux / injected < 1.1, (flux, injected)
