CHARIS Data Processing Pipeline
===============================

Reduction and analysis of SCExAO/CHARIS IFS data cubes. Two implementations live here:

* `charis_dpp/` — the **Python package** (no IDL licence needed). See `docs/WORKFLOW.md`
  and the tutorial `docs/tutorial.ipynb`.
* `idl/` — the original IDL pipeline, with the fixes listed in `BUGS.md`.

Python quick start
------------------

```bash
pip install -e ".[dev]"
charis doctor
```

In a directory with extracted cubes (`CRSA*_cube.fits`):

```bash
charis newobs 20230730 "HD 1160" lowres     # writes HD1160_lowres.toml
charis imprep HD1160_lowres.toml
charis subtract-sky HD1160_lowres.toml
charis register HD1160_lowres.toml          # --method xcorr if there are no satellite spots
charis specphot-cal HD1160_lowres.toml
charis make-emp-psf HD1160_lowres.toml
charis imrsub HD1160_lowres.toml --prad
charis adialoci HD1160_lowres.toml
charis detect HD1160_lowres.toml
charis fwdmod-planet HD1160_lowres.toml --targfile reduc/proc/final_det.targ
charis extract-spectrum HD1160_lowres.toml --coords X Y --throughput synth_throughput0.txt
```

`charis autoreduce` runs the chain; `charis <step> --help` lists the options. Parameters come
from the dataset `.toml` (edit it) and, for instrument constants, from the dated table in
`charis_dpp/data/instrument.toml`.

Tests: `pytest -m "not slow"` (fast), `pytest -m slow` (synthetic end-to-end chains).

---

IDL pipeline (original README)
==============================

CHARIS Data Processing Pipeline [documentation under construction!]
===============================================================

The CHARIS Data Processing Pipeline is a full end-to-end, supported pipeline for processing and analyzing SCExAO/CHARIS integral field spectrograph data cubes.   

Requirements: 
-------------------------------
- the distributed CHARIS-DPP package (available here)
- an IDL license (code is not pre-compiled).   Note that a future release (~late 2021) will be available in Python 3.
- the IDL Astronomer User's Library.  http://idlastro.gsfc.nasa.gov/
- SAOImage/DS9 or similar image display program (STRONGLY recommended but not strictly required)

  Note that the CHARIS DPP requires rectified CHARIS data cubes as input.   It does not create the cubes themselves from raw CHARIS reads.   To create CHARIS data cubes, see Tim Brandt's cube extraction pipeline, here: https://github.com/PrincetonUniversity/charis-dep 
  
**Installation**: 
-------------------------------
To set up this package you need to do the following ...

A. With Git

1. cd to your preferred directory
2. git clone https://github.com/thaynecurrie/charis-dpp.git

B. Manually

1. Download and copy entire zip'd package to your preferred directory path: (e.g. cp [package_version].zip [your/preferred/path/])
2. Unzip the package (i.e. cd to [your/preferred/path/] and then "unzip [package_version].zip" from standard Terminal command line)

THEN ...
3. Add that directory path to your IDL PATH, which is usually defined in a .bash_profile or .cshrc file.  

For example, if your preferred directory path is charis_dpp and you are using a **tcsh** shell, then the following line should be added to your .cshrc file:

      setenv IDL_PATH +$IDL_DIR/lib:+$HOME/[your/preferred/path/charis_dpp_v7]
      
If your preferred directory path is charis_dpp and you are using a **bash** shell, then the following line should be added to your .bash_profile file:

      export IDL_PATH=+$IDL_DIR/lib:+HOME/[your/preferred/path/charis_dpp_v7]

4. "source" your .cshrc or .bash_profile file.

5. Navigate to your CHARIS DPP package directory.   Under the "setup" subdirectory, edit line 13 of charis_path.pro, to change the directory path to your full charis PATH

     (e.g. if you unpack this package as charis_dpp, then edit line 13 as follows: charispath='[path to CHARIS-DPP]/charis_dpp/'   )


To test whether your installation was (likely) successful, start IDL prompt from anywhere on your computer and run a program called charis_test:
`IDL> charis_test`

You should see a print out of the path to charis_test and the following message:
`CHARIS-DPP Installation likely successful`


Documentation:
-----------------

My plan is to provide full documentation on readthedocs at some point in the near future.   The full documentation will also included tutorials to cover standard cases for data reduction.   Stay tuned.

Overview:
-----------
-__Pipeline Workflow__: The typical pipeline execution has the following key steps:
 1. __Preliminaries__: creation of a parameter file ([string].info) populated by basic information about your target and default pipeline settings, fits header manipulation
      - charis_newobs, charis_imprep
     
 2. __Sky Subtraction__ : removes thermal emission from each data cube
      - charis_subtract_sky
 3. __Data Cube Registration__: recenters each slice of each cube to coordinate [100,100] (indexed from 0). Creates a log file of parallactic angle.
      - charis_register_cube
 4. __Model PSF Creation__: creates an empirical PSF from the satellite spots
      - charis_makeemppsf
 5. __Spectrophotometric Calibration__: provides absolute flux calibration for each cube using either the satellite spots or an unsaturated, unocculted data cube
      - charis_specphot_cal
 6. __Spatial Filtering__: subtracts off either a radial profile from each slice or a moving-box median filter.   Also generates a sequenced-combined classical ADI cube
      - charis_imrsub
 7. __PSF Subtraction__: subtracts off the stellar PSF using A-LOCI or KLIP employed with different combinations of angular and/or spectral differential imaging or reference star differential imaging
      - charis_adialoci, charis_sdialoci, charis_adiklip, charis_rdiklip
 8. __Forward-Modeling: Throughput Corrections and Attenuation Maps__: uses forward-modeling to determine the signal loss per channel of a planet at a given position or the azimuthal average signal loss for planets over a range of separations (i.e. an attenuation map).
      - charis_aloci_fwmod_planet, charis_klip_fwdmod_planet, charis_aloci_attenmap, charis_klip_attenmap

 9. __Spectral Extraction__: Extracts the spectrum of a planet from a PSF-subtracted CHARIS data cube at a given location with or without a throughput correction
       - charis_extract_1d_spectrum

10. __Contrast Curves__: Calculations the 5-sigma contrast curve over the CHARIS bandpass and in individual channels
      - charis_calc_final_contrast

For __Disk sources__, steps 8--9 are replaced by a disk forward-modeling module
      - charis_aloci_disk_fwdmod, charis_klip_disk_fwdmod



Quick Start Guide:
-----------
Coming Soon.

A note on satellite spot contrasts:
-----------------
The CHARIS satellite spot contrast measurement has changed over the years. The pipeline currently does not have a decision tree that leads you to more recent values -- you may have to update them yourself, whether here or in [pyklip-charis](https://pyklip.readthedocs.io/en/latest/instruments/CHARIS.html). A guide to the best of my knowledge is below, and can be cited as T. Currie priv. comm.

Using whatever the standard modulation frequency was:

Before Aug 30 2017: 9.97e-04

Aug 30 2017 to Dec 31 2021: 2.72e-03 (a Lozi+2018 calibration differs by ~8%, 2.94e-3)

Dec 31 2021 to May 20 2022: 1.36e-03

After May 20 2022: 1.25e-03

If using 500 Hz modulation frequency:

Before Dec 31 2021: ???????? (presumably no data)

Dec 31 2021 to May 20 2022: 1.80e-03

May 20 2022 to Oct 20 2022: ???????? (presumably no data)

After Oct 20, 2022: 2.56e-03

The dated epoch table in `charis_dpp/data/instrument.toml` and `idl/setup/charis_get_constant.pro` currently encodes only the 2017-08-30 recalibration; the later values above are not yet in the table (see BUGS.md).
