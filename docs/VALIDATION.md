# Validation: old IDL → new IDL → Python on real 2023-07-30 CHARIS PDI data

Status (2026-11-06, 01:30 HST): **HD 236894 and Europa three-way runs complete;
every port defect they exposed is fixed on `integrate-port`**. Every number below comes
from `~/charis-val/compare/*/compare.csv` on mueller or from the aperture
photometry script quoted in the Results section.

## Why

The Python port (`charis_dpp`) was validated only on synthetic data, and the
65 IDL patches on this branch had never been compiled. The only prior real
reduction (Europa, 2023-07-30, June 2024) was a hybrid: IDL flat-field and sky
subtraction, notebook-Python registration with hand-set offsets, notebook-Python
spectrophotometric calibration (whose sky subtraction never took effect), then
IDL HWP matching and `charis_pdi`. Its intermediate products are gone. It is a
sanity check, not a reference. The reference is produced here by running the
pre-port IDL code (tag `old-idl-2025-04`) end to end.

Two questions, in order:

1. **Does new IDL match old IDL?** Same `.info`, same steps, tag
   `old-idl-2025-04` vs branch `integrate-port`. Expected: bit-identical except
   where a BUGS.md row or a deliberate change explains the difference.
2. **Does Python match new IDL?** Same data via `charis migrate-info`. Expected:
   agreement at `tests/regression/tolerances.toml`, with the documented
   deliberate differences.

Tool: `charis compare-trees A B --mode identical|port [--allow-list f] [--out d]`
(`charis_dpp/testing/compare_trees.py`), which pairs `prep/`, `reg/`, `proc/`,
`reduc*.log` and `hwp_info.txt` by relative path.

## Layout on mueller (as built)

```
~/charis-val/
  old-idl/      git worktree at tag old-idl-2025-04   (IDL 8.9, ~/idl/idl89; charis_path.pro edited, uncommitted)
  new-idl/      git worktree at branch integrate-port (local branch integrate-port-val tracking it)
  datasets/
    frames-2023-07-30.md, frames-2023-07-30-headers.csv    frame lists (copies in docs/validation/)
    calib_2023-07-30/           CHARIS-DEP 2.0.1 buildcal on CRSA00123405 (1550 nm) with darks 123417-420,
                                oversampled; recovered spot shift (-5.377487, 3.110374, -0.000860) equals the
                                CAL_DX/DY/DPHI of the 2023 cubes to six digits; charis_2023-07-30.ini = sample.ini + calibdir
    HD236894_bb/
      data/raw/                 symlinks to the 32 cubes CRSA00123552-583 extracted in 2023
      pdi_flats/                CRSA00123412-416 extracted here with the calibration above
      HD236894_low.info         hand-built (see below); byte-identical copy in old/, new/, py/
      old/ new/ py/             one working directory per run (reduc/, reduc*.log, hwp_info.txt, logs)
      old_nomm/ new_nomm/       charis_pdi with /no_mm_ip_corr,/no_mm_crosstalk_corr on the same reg/ (symlink)
      py_regtest/               Python registration run on the new-IDL flat-fielded frames
      py_pditest_mm/, _nomm/    Python charis pdi run on the new-IDL registered frames
    Europa_bb/
      data/raw/                 CRSA00123628-855 being extracted (220 science + 8 sky)
      pdi_flats/ -> HD236894_bb/pdi_flats cubes;  Europa_low.info in old/ new/ py/
  compare/                      compare-trees outputs: HD236894_<A>_vs_<B>/{compare.csv,compare.md,stdout.txt}
```

IDL is started per run with `IDL_PATH="+<worktree>/idl:+$HOME/IDLAstro/pro:<IDL_DEFAULT>"`
(`~/IDLAstro` was cloned on mueller and on the laptop: no astrolib existed on
either, and `robust_poly_fit` is called from `charis_pdi` and
`charis_register_cube`) and, for the new tree, `CHARIS_DPP_PATH=~/charis-val/new-idl/idl`.
Python: conda env `charisenv` (3.11), `pip install -e ~/charis-val/new-idl`.

Compile gate (`charis_compile_check`): expected output on IDL 9.1 (laptop) and
8.9 (mueller) is 329 files, **1 failed** (`subroutines/craigm/xwindow.pro`,
third-party, needs Coyote's `pswindow`; fails on the old tag too; unused), 4
skipped (`@include` fragments) and two `DUP` lines (`charis_getcenter` in
`reg/` and `subroutines/`, `gbfunc` in `gbfunc.pro` and inside
`charis_calc_spec_covar.pro`; BUGS.md #325). Both machines give exactly this.

## Inputs and decisions

| item | value | why |
|---|---|---|
| `.info` | built with the old tree's `param,/set` calls copied from `charis_newobs` (`docs/validation/HD236894_low.info`); SIMBAD values from the Harvard mirror (RA 28.05713, Dec +58.43509, plx 0.3334 mas, O8V, J/H/K 8.936/8.877/8.909) | CDS Strasbourg was unreachable from both machines, so `charis_newobs` died in `charis_querysimbadf`; HD 236894 is not in Hipparcos, so `newobs` would have written sptype -999 and zero parallax |
| `fnum_sat`, `fnum_sky` | `9-32`, `''` | **sequence indices**, not CRSA numbers: `imprep` renumbers the raw cubes in sorted order (`n0001e` = 123552 … `n0032e` = 123583), so 123560-583 is 9-32. Both `imprep`s rewrite the lists (IDL: all 32 as sky; Python: all 32 as science), so they are set again after `imprep` |
| flats | `pdi_flats/` = CRSA00123412-416, 5.9 s, no lamp-off frames, no `darkflatdir` | only PDI flats of the night |
| sky subtraction | skipped on all sides | no sky frames for the standard |
| registration | `method=1` (`psf`, gcntrd on the unsaturated PSF), `polcent_offset` default (31.758, 16.542), `fitgauss` off | no satellite spots anywhere that night (`X_GRDST=OFF`) |
| specphot | **skipped on all sides** | the old tree's `/unsat` path dies on broadband data (BUGS.md #327); P and AoLP do not depend on it |
| HWP angles | `/use_headers` (IDL reads `RET-ANG1` from the raw cube's 3rd extension) / `--source header` (Python `imprep` copies it to the primary header) | frames-2023-07-30.md; the cubes' 3rd extension is the raw primary header |
| `charis_pdi` | `/meanadd, /save_derot_seq`, and a second run with `/no_mm_ip_corr, /no_mm_crosstalk_corr`; `angoffset` 2.2 via `charis_set_epoch` (MJD 60155) | isolates the Mueller term |
| `flip_matrix` | false on both sides | fork behaviour; BUGS.md #308 |

## Results: HD 236894

### Old IDL vs new IDL (identical mode)

| step | files | result | attribution |
|---|---|---|---|
| `imprep` (`prep/n####e.fits`) | 32 | **identical** (rel_rms 0) | |
| `pdi_flatfield` (`master_pdiflat`, `n####e_flat`) | 24 | **identical** | |
| `pdi_register_cube` | 24×3 | old tree **crashes** at `charis_pdi_register_cube.pro:135` (`LSHIFT` out of range) after registering both beams; its per-slice centroids and `reduc_left/right.log` are identical to the new tree's | BUGS.md #336: the old wrapper needs `SHIFT_n` keywords that only the sat-spot branch writes; the fork writes them on every branch. The old native `leftreg/rightreg` (pre-NaN re-shift) are kept in `old/reduc/reg_oldnative/`; the old downstream steps ran on the new tree's `reg/` |
| `pdi_hwp_match` (`hwp_info.txt`) | 24 rows | **identical** (6 cycles) | only the integer field width differs (nbrlist returns long, #285) |
| `charis_pdi`, Mueller off (`old_nomm` vs `new_nomm`) | 132 | 106 pass, **2 fail at rel_rms 1.6e-6** (`pdi_aolp.fits`, `pdi_c0005_u_derot.fits`; max diff 0.012 deg / 0.013) | floating-point noise just above the 1e-6 threshold; allow-list candidates |
| `charis_pdi`, Mueller on (`old` vs `new`) | 24 | **24 fail**: `pdi_q` 9 %, `pdi_u` 41 %, `pdi_aolp` 57 % rel_rms | **deliberate**: the fork's McIntosh+2026 Mueller model (derotator retardances from `mcintosh_2026_spie_wavelength_fits.csv` instead of the Hart+2021 mirror model, HWP thicknesses 1.653/1.291 mm instead of 1.623/1.268, Wollaston efficiency from the CSV instead of 1, `dtheta_hwp` −0.04 instead of −0.002 deg, `dtheta_der` −0.02 instead of −0.0118 deg; `m3_model` constants unchanged for 2023). Commit `dfb9331`. |

Verdict: **the core PDI arithmetic is unchanged between the old tag and this
branch**; every difference in the corrected products is the Mueller-model
update, and the old tree could not have run this dataset's registration or
spectrophotometric calibration at all.

### New IDL vs Python (port mode, `tests/regression/tolerances.toml`)

Step-isolated comparisons (each Python step run on the new-IDL products of the
previous step):

| step | files | result | notes |
|---|---|---|---|
| `imprep` | 32 | identical | |
| `pdi_flatfield` (`master_pdiflat`, `n####e_flat`) | 25 | **identical** (rel_rms 4e-8) | after #326, #328, #334 the ROI mask is the IDL one pixel for pixel; the #216 near-zero guard is now opt-in (`[pdi] flat_floor`) because it changed 109 edge pixels per channel |
| `pdi_register_cube` on the IDL flats (`py_regtest`) | 128 | **128 pass, 0 fail**: `leftreg` rel_rms 1.7e-6, `reg` 9e-5; per-slice `SHIFT_n` agree to ≤ 0.001 px | after #329 (gcntrd) and #333 (robust_poly_fit) |
| `charis pdi` Mueller off on the IDL `reg/` (`py_pditest_nomm`) | 145 | **132 pass, 0 fail**; pre-derotation `c####_{q,u,iq,iu}` agree to 2.5e-8, derotated to 4e-6 | derotation kernel (`rotate_image` = `rotat` port) |
| `charis pdi` Mueller on (`py_pditest_mm`) | 177 | **107 pass, 1 fail** (`pdi_aolp.fits` cube, rel_rms 9e-3: the angle wraps at ±90 deg where Q and U are ~0) | after #330, #331, #332, #335; Mueller rows agree to 3e-4 before #335, better after |

Full Python chain (`py`, flats → registration → HWP → pdi) vs new IDL:
**146 pass, 18 fail** of 177 (13 "missing" are the Python-only
`_derot_uncorr` cubes). `prep/` and `reg/` all pass (`leftreg` 1.7e-6,
`rightreg` 1.8e-4: the right-beam centres differ by ≤ 0.001 px on a few
frames). The 18 failures are the per-cycle `q/u_derot` cubes (2-9e-3),
`pdi_q`/`pdi_qphi` (6e-3), `pdi_pi` (7e-3), `pdi_u` (4e-2) and the `pdi_aolp`
cube (angle wrap): Q and U are ~1 % of I, so a 1e-4 registration residual in I
becomes 1e-2 in Q/U per pixel. Every `_collapsed` product passes, and the
aperture-integrated P and AoLP agree to 0.001 % and 0.01 deg (table below).

### Recovered polarisation of the standard

Aperture sums (r < 12 px) of the collapsed Q, U and of the mean of the I_Q/I_U
double sums; P in per cent, AoLP in degrees (pipeline frame, north-up, before
any catalogue comparison):

| run | P % | AoLP | q % | u % | P % per channel 0 / 5 / 10 / 15 / 21 |
|---|---|---|---|---|---|
| old IDL, Mueller on | 1.287 | −89.11 | −1.286 | −0.040 | 1.99 1.73 1.09 0.93 0.18 |
| new IDL, Mueller on | 1.360 | −85.63 | −1.344 | −0.206 | 2.00 1.74 1.18 0.98 0.26 |
| Python, Mueller on (full chain) | 1.360 | −85.62 | −1.344 | −0.207 | 2.01 1.74 1.18 0.98 0.26 |
| Python on IDL reg/, Mueller on | 1.360 | −85.63 | −1.344 | −0.206 | 2.00 1.74 1.18 0.98 0.26 |
| old IDL, Mueller off | 1.861 | 54.70 | −0.618 | 1.755 | 2.66 2.48 1.67 1.44 0.53 |
| new IDL, Mueller off | 1.861 | 54.70 | −0.618 | 1.755 | 2.66 2.48 1.67 1.44 0.53 |
| Python, Mueller off | 1.861 | 54.70 | −0.618 | 1.755 | 2.66 2.48 1.67 1.44 0.53 |

Python and new IDL agree to 0.001 % in P and 0.01 deg in AoLP. The old-vs-new
difference (1.29 % / −89 deg vs 1.36 % / −86 deg) is the Mueller-model update.
Catalogue: the observing target list gives P = 3.92 % (Heiles 2000, optical)
and no position angle for HD 236894 (H = 8.877); no near-infrared standard value
was found in the literature (searched 2026-11-06). For interstellar polarisation
with a Serkowski curve, P(H)/P(V) ≈ 0.3, so ≈ 1.2 % is expected around 1.6 μm,
falling towards K: the measured 2.0 % (1.16 μm) → 1.18 % (1.63 μm) → 0.26 %
(2.37 μm) is of the right size and shape. A proper check needs the Heiles 2000
position angle and the pipeline's AoLP convention; left for the user.

## Port defects found on this dataset (all fixed on `integrate-port`)

| BUGS.md | what | effect before the fix |
|---|---|---|
| #326 | `roirange` widths taken as half-widths | master flat all ones |
| #328 | FOV mask rotated the wrong way | flat normalisation off by 1-6 % |
| #329 | `gcntrd` was a non-linear fit, not astrolib's | centres 0.05-0.1 px off IDL |
| #330 | `D_IMRANG` read from the wrong header, silently 0 | corrected P/AoLP 1.56 % / −71.5 deg instead of 1.36 % / −85.6 deg |
| #331 | HWP Sellmeier instead of the precomputed table; `np.interp` clamps where `interpol` extrapolates | Mueller terms off at 2e-4 |
| #332 | `getrot` port ignored the CD determinant sign | every north-up product rotated by 2× the field angle (6.6 deg) |
| #333 | `robust_poly_fit` weights and stopping rule | fitted centres 0.01 px off |
| #334 | analytic ROI vs IDL's `rotat`-based one | flat normalisation off by 0.06-0.35 % |
| #335 | (IDL) single-precision literals in `hwp_model.pro`; the port mirrors them | HWP retardance 0.04 deg |

Old-tree defects recorded: #327 (`charis_specphot_cal_unsat` has no
`'Broadband'` case), #336 (PDI registration wrapper needs `SHIFT_n` keywords
the old `psf` branch never writes). Both were already fixed by the fork.

## Deliberate differences (allow-list candidates)

| product | reason |
|---|---|
| `proc/pdi_*` with Mueller on, old vs new | McIntosh+2026 Mueller model (commit `dfb9331`) |
| `proc/pdi_aolp.fits`, `proc/pdi_c0005_u_derot.fits`, old vs new, Mueller off | 1.6e-6 rel_rms floating-point noise |
| `proc/pdi_c####_{q,u}_derot.fits`, `pdi_q`, `pdi_u`, `pdi_pi`, `pdi_qphi` cubes, IDL vs Python | registration residual (≤ 0.001 px) amplified by I/Q ≈ 100; collapsed products and aperture P/AoLP agree |
| `proc/pdi_aolp.fits` cube, IDL vs Python | angle wrap where Q, U ≈ 0 |
| `reduc_left/right.log` XC, YC columns | IDL logs the re-measured post-shift centre (≈ 99.98), Python the pre-shift centre (≈ 66, 82); HA/PA agree |

## Results: Europa (2023-07-30, 220 science + 8 sky frames, 13.3 s broadband PDI)

Inputs: `Europa_low.info` (`docs/validation/`; `fnum_sat` 1-220, `fnum_sky`
221-228, sptype −999, coordinates of frame 123740), the 228 cubes extracted on
mueller with the same CHARIS-DEP calibration as the flats, lamp flats as above,
`pdi_subtract_sky /scalesky` with the 8 sky frames.

### Stage 1: imprep, flat-field, sky subtraction

| pair | mode | files | result |
|---|---|---|---|
| old IDL vs new IDL | identical | 676 | **676 pass** |
| new IDL vs Python | port | 676 | **676 pass** |

### Stage 2: registration, HWP matching, `charis_pdi`

Neither IDL tree can register an extended source without satellite spots
(methods are `sats`, `psf`, `calsat`), exactly as in 2023, so registration is
Python-only: `charis pdi-register --method xcorr` with per-beam centroid boxes
(`[pdi.xcorr] left_box = [43, 93, 58, 108]`, `right_box = [107, 157, 92, 142]`),
i.e. the 2023 notebook method. Its `reg/` feeds all three `charis_pdi` runs.

| pair | mode | files | result | attribution |
|---|---|---|---|---|
| HWP table, old vs new IDL (before #337/#338) | — | 220 rows | identical, 55 cycles, **reproduces the 2023 c0019-c0036 pairing anomaly exactly** (cycle 19 = 123704+123709, …, 123705 in cycle 36) | BUGS.md #338: the mid-exposure seconds go negative for frames whose `UTC-TIME` seconds field is < dit/2 and `ten()` returns a negative hour; #337: frame CRSA00123812 has no `RET-ANG1` and was taken as 0 deg, making a bogus 55th cycle |
| HWP table, new IDL (patched) vs Python | — | 220 rows | **identical**: 54 complete cycles, 123704 paired with 123705, 123812 (and its 45 deg partner) left out with a warning | |
| `charis_pdi` Mueller off, old vs new | identical | 892 | 756 pass, 132 fail (per-cycle cubes of the re-paired cycles ≥ 19; collapsed `pdi_q` 5.5e-3, `pdi_u` 4.3e-3) | #337/#338 only: before the patch the same comparison gave 817 pass and 75 fails at ≤ 5.6e-6 |
| `charis_pdi` Mueller off, new vs Python | port | 997 | **888 pass, 0 fail** (`pdi_q_collapsed` 9e-5, `pdi_u_collapsed` 6e-6) | |
| `charis_pdi` Mueller on, new vs Python | port | 1673 | **1564 pass, 0 fail** (`pdi_q_collapsed` 2e-5, `pdi_aolp_collapsed` 3e-4) | 109 "missing" are the Python-only `_derot_uncorr` cubes |
| `charis_pdi` Mueller on, old vs new | identical | 1568 | 1446 pass, 122 fail (`pdi_q` 17 %, `pdi_u` 49 %, `pdi_aolp` 76 %) | Mueller-model update (as for HD 236894) plus #337/#338 |

### Stage 3: disk-integrated polarisation, trends, 2023 comparison

Disk aperture r < 25 px on the collapsed products; per-cycle values from the
`c####_{q,u,iq,iu}_derot` cubes (`docs/validation/per_cycle_qu.csv`,
`europa_qu_vs_time.png`); bootstrap over cycles (2000 draws).

| run | cycles | q % | u % | P % | AoLP | dq/dt, du/dt per hour (significance) |
|---|---|---|---|---|---|---|
| new IDL, Mueller on | 54 | −0.067 ± 0.004 | −0.097 ± 0.003 | 0.118 ± 0.003 | −62.2 ± 0.9 | +0.039 (3.3σ), +0.026 (2.6σ) |
| Python, Mueller on | 54 | −0.067 ± 0.004 | −0.097 ± 0.003 | 0.118 ± 0.004 | −62.2 ± 0.9 | +0.039 (3.3σ), +0.026 (2.6σ) |
| old IDL, Mueller on (Hart+2021 model, old pairing) | 55 | −0.096 | −0.043 | 0.105 | −78.0 | +0.041 (3.6σ), +0.012 (1.3σ) |
| new IDL / Python, Mueller off | 54 | +0.543 ± 0.005 | +0.722 ± 0.007 | 0.903 ± 0.004 | 26.5 ± 0.3 | +0.094 (9.5σ), −0.154 (13σ) |

Reading: the instrumental term is ≈ 0.9 % and rotates strongly with time (10-13σ
trends when the Mueller correction is off); the McIntosh+2026 model removes
it to 0.12 % with residual trends of ~3σ (0.04 %/h), which is Kellen's "model
or non-sidereal bookkeeping slightly off" regime and worth a closer look with
the 2025 nights. Disk-median Q is negative in the corrected products (−16
counts; the 2023 products had −30 with a different flux scale), positive
before correction.

Against the 2023 hybrid products (same frames, notebook registration with hand
offsets, old pairing, Python specphot): the total-intensity images correlate
at 0.989 over the disk after a 1.3/2.1 px centre shift, Q at 0.90, U at 0.82
and PI at 0.79, with a flux ratio of 2.06 (the 2023 products were
spectrophotometrically scaled). The 27 cycles that pair differently are
cycles 19-36 (#338) and 46-54 (#337).

## Overall verdict

1. **New IDL = old IDL** wherever the old code runs: `imprep`, flats, sky
   subtraction, HWP tables and the Mueller-off `charis_pdi` are identical to
   ≤ 6e-6. The only differences are deliberate (McIntosh+2026 Mueller model) or
   fixes of old defects found here (#327, #336, #337, #338). The old code could
   not register either dataset.
2. **Python = new IDL**: identical `imprep`/flats/sky/HWP tables; registration
   to ≤ 0.001 px; `charis_pdi` products within port tolerance on Europa
   (0 failures of 2452 files) and, on the star, everywhere except the per-pixel
   Q/U cubes whose residual is the 1e-4 registration difference amplified by
   I/Q; disk-integrated P and AoLP agree to 0.001 % and 0.1 deg on both targets.
   Nine port defects (#326, #328-#335 mirror) and one matcher bug were found and
   fixed on the way; none was visible on synthetic data.
