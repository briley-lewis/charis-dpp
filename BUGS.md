# CHARIS-DPP bug inventory

This file lists every defect found in a line-by-line audit of the IDL pipeline
(`idl/`, excluding the vendored `subroutines/craigm` and `subroutines/textoidl`).
The audit was done by reading the code; no IDL interpreter was available, so IDL
patches listed as applied have been syntax-checked by eye only and should be
exercised in IDL before being relied on.

Line numbers were taken from the pre-sync tree (`briley-lewis/charis-dpp@main`,
commit 4af0c1e); the current `idl/` is the autoformatted
`thomasmcintoshh/charis-dpp-mm-updates` tree, so lines shift by a few in most
files. Section 14 was taken from the fork tree directly.

Columns: **Sev** — CRASH (won't compile / always errors), WRONG (runs, wrong
output), GEOM (hard-coded 201×201 / centre 100 / 22-channel / pixel-scale
assumption), DEAD (unused keyword, unreachable code, debug leftovers).
**IDL patch** — Y: patched in `idl/` on this branch (the patch carries a
`; [BUGS.md #n]` comment); N: documented only, fixed in the Python port.
Every row is fixed in the Python port.

Status of IDL patches is tracked in the "IDL patch" column: `Y` = applied on
this branch (65 rows; grep `BUGS.md #<n>` in `idl/` to find each
patch), `Y (pending)` = a local IDL fix is intended but not yet applied
(220 rows), `N` = documented only.  Every row is fixed in the Python
port regardless.  Patches were syntax-checked by eye (no IDL available) and
should be exercised in IDL before being relied on.

The IDL constants were also restructured: `setup/charis_get_constant.pro`
now holds the dated epoch table (`mjd=` keyword; `tests/unit/test_idl_parity.py`
checks it against `charis_dpp/data/instrument.toml`), `setup/charis_path.pro`
reads `$CHARIS_DPP_PATH`/`$CHARIS_DPP_BTSETTL`, `cals/get_charis_sat_contrast.pro`
uses the epoch table and accepts a per-channel manual vector, and the three
routines that were never shipped (`get_xycind`, `calc_centroid_spots4`,
`medfitsme`) now exist under `subroutines/`.

## Defect inventory

Severity: **CRASH** = won't compile or always errors on that path · **WRONG** = runs but gives wrong output · **GEOM** = hard-coded 201×201 / centre 100 / 22-channel / pixscale assumption (breaks other modes or sizes) · **DEAD** = unused keyword, unreachable code, leftover debug output. Line numbers are from the older tree; the fork is autoformatted so they shift. "IDL patch": Y = patched in `idl/`, N = documented only (fixed in Python). Every row is fixed in Python.

### 1. prelims/ (newobs, imprep, health, set/add/del_files, sxpar_charis)

| # | File:line | Sev | What is wrong | Fix | IDL patch |
|---|---|---|---|---|---|
| 1 | charis_newobs.pro:135,139,150,159,168 | WRONG | `if (plxset gt 0)` — `where` returns index ≥0, so a tag at index 0 is treated as missing (plx/pm/sptype set to 0) | `ge 0` | Y |
| 2 | charis_newobs.pro:152,161 | WRONG | `epmra`/`epmde` read without checking the tag exists → crash if Hipparcos row lacks them | guard with `epmraset ge 0` | Y |
| 3 | charis_newobs.pro:122 vs 148 | WRONG | missing-distance sentinel is 0 in one branch, −999 in the other | use −999 consistently | Y (pending) |
| 4 | charis_newobs.pro:194-196 | CRASH | 2MASS result not checked to be a struct; `s._r` crashes when VizieR returns nothing | `size(s,/type) eq 8` guard, write −999 mags | Y |
| 5 | charis_newobs.pro:89 | DEAD | `stop` when star not found | `message` | Y |
| 6 | charis_newobs.pro:95 | DEAD | first 2MASS query result `s2mass` never used (queried twice) | remove | Y (pending) |
| 7 | charis_newobs.pro:216 | WRONG | writes `procdir='./proc/'` but pipeline uses `./reduc/proc/` | fix string | Y |
| 8 | charis_newobs.pro:209-211 | WRONG | canned `fnum_sat='1-100'`, `fnum_sky='101-103'` survive if imprep finds no matching frames | imprep writes empty string when none | Y (pending) |
| 9 | charis_newobs.pro | DEAD | `cfa` keyword, `obj_main` unused; comment "INCOMPLETE" | document | N |
| 10 | charis_imprep.pro:16 | CRASH | `if strlen(possiblenames) gt 4` fails with >1 `.info` file (array in `if`) | `n_elements(...) eq 1` | Y |
| 11 | charis_imprep.pro:153-199 | CRASH | `ra`,`dec` only set in the "no LST in header" branch; `rap=ra` at :194 is undefined when header has ST/LST | always read RA/DEC | Y (pending) |
| 12 | charis_imprep.pro:199 | WRONG | header LST is an `hh:mm:ss` string but never passed through `ten()` → `ha=lst-rap/15` wrong or errors | `ten()` | Y (pending) |
| 13 | charis_imprep.pro:159-163 vs 199 | WRONG | if header has HA, that value is kept and the freshly computed `ha` discarded | compute once, write once | Y (pending) |
| 14 | charis_imprep.pro:241 | WRONG | output numbering is `i+1` in `file_search` order; any extra FITS file in `raw/` shifts every number | number from `ORIGNAME`/frame number; fallback sorted | Y (pending) |
| 15 | charis_imprep.pro:127 vs 245 | WRONG | obslog header has 5 labels, rows have 6 columns | add FILTER label | Y (pending) |
| 16 | charis_imprep.pro:62 | WRONG | `file_test('CRSA*fits')` doesn't expand wildcards → raw files never moved | `file_search` | Y |
| 17 | charis_imprep.pro:264-265 | CRASH | `writefits,/compress` then `/append` — appending to gzip unsupported | write then gzip | Y (pending) |
| 18 | charis_imprep.pro:139,147 | WRONG | `exp1time=1.475` hard-coded; `coadds=lastrd-firstrd+1` may be off by one for CDS readout | document; Python reads from header with override | N |
| 19 | charis_imprep.pro:42,46,57 | DEAD | `prefname`,`reducdir` not keywords; `./data/prep/` created but never used; lines 73-92 dead after `goto,skipsort` | remove | Y (pending) |
| 20 | charis_imprep.pro:211,213 | DEAD | `coadds` written with `f9.5`; `a10` ignored for strings | `i4` | Y (pending) |
| 21 | charis_check_health.pro:19-20 | DEAD | `'good'` never returned (always overwritten) | remove branch | Y (pending) |
| 22 | charis_check_health.pro:30,38-39 | WRONG | comments say median & 10×, code uses mean & 3×; empty `inner`/`outer` gives `avgimage[-1]` | fix comment, add count guard | Y (pending) |
| 23 | charis_set_files.pro:29,46 | CRASH | with one frame, `fnumber[nframes-2]` = `[-1]` | handle `nframes eq 1` | Y |
| 24 | charis_set_files.pro + charis_add_files.pro:12 | WRONG | input assumed sorted/unique; `add_files` passes unsorted with possible duplicates | `sort(uniq())` | Y |
| 25 | charis_set_files.pro:15-18 | WRONG | when nothing matches, parameter keeps its old/default value | write `''` | Y |
| 26 | charis_del_files.pro:11 | CRASH | uses `nbrlist(filetoadd)` instead of `filestodel` → `/stringinput` broken | rename | Y |
| 27 | charis_del_files.pro:14-17 | CRASH | no match → `todelete=-1` → `remove` errors | guard | Y |
| 28 | sxpar_charis.pro:190 | CRASH | `/justfirst` on a vector search with no match sets `nfound=-1, matches=1` → `hdr[-1]` | guard | Y |
| 29 | HR8799_low.info, HR8799_low_example.info | DEAD | identical; lack `FNUM_LAD`, `PLANMOD` that newobs writes | regenerate example | Y (pending) |

### 2. skysub/

| # | File:line | Sev | What is wrong | Fix | IDL patch |
|---|---|---|---|---|---|
| 30 | charis_subtract_sky.pro:141 | GEOM | `dist_circle,gg,201` hard-coded, centred at 100 | use image size / PSFCENT | Y (pending) |
| 31 | charis_subtract_sky.pro:146,152 | WRONG | empty `good` → `median(targslice[-1])` | count guard | Y (pending) |
| 32 | charis_subtract_sky.pro:147,154-155 | WRONG | `/even` used for sky median but not target median | consistent | Y (pending) |
| 33 | charis_subtract_sky.pro:77-85,102-104 | DEAD | `skytime` computed, exposure scaling commented out; `medbox`, `fwhm`, `dim`, `d0` unused | document | N |
| 34 | charis_subtract_sky.pro:173 | DEAD | `mastersky.fits` rewritten every loop iteration | move out of loop | Y (pending) |

### 3. reg/ (register_cube, getcenter, getcenter_pdi, check_psfquality)

| # | File:line | Sev | What is wrong | Fix | IDL patch |
|---|---|---|---|---|---|
| 35 | charis_register_cube.pro:112 | CRASH | checks `method eq 'calsats'` but help documents `'calsat'`; string falls through `case` | accept both | Y |
| 36 | charis_register_cube.pro:169 | WRONG | `n_elements(test) gt 1` — a single `*skysub.fits` file is ignored, `e` files used instead | `count=` | Y (pending) |
| 37 | charis_register_cube.pro:121 | GEOM | `refslice=10` hard-wired ("ONLY GOOD FOR BROADBAND"); J has 15 channels | derive from wavelengths | Y (pending) |
| 38 | charis_register_cube.pro:389 | WRONG | `/smask` loop runs over `n_elements(x)` with `x` undefined → 0 iterations, mask never applied | `xsource` | Y |
| 39 | charis_register_cube.pro:538-541,745-748 | WRONG | on centroid failure `deltax=-9999` and the cube is still shifted by ~9999 px and written | skip/flag frame | Y (pending) |
| 40 | charis_register_cube.pro:564,942 | WRONG | NaN fallback rescales slice-0 geometry with `wvl[refslice]/wvl[s]` (ratio relative to refslice, not slice 0) | `wvl[0]/wvl[s]` | Y (pending) |
| 41 | charis_register_cube.pro:994 | WRONG | y guess uses `PSFcensref[0,s]` (x) instead of `[1,s]` | fix index | Y |
| 42 | charis_register_cube.pro:903,916 | CRASH | `stop` statements in method 2 (`:916` unconditional) → method 2 always halts | remove | Y |
| 43 | charis_register_cube.pro:695-700 | WRONG | `/revise,/keepprev` uses `PSFcens` that was zeroed at :663 → forces centroid to (0,0) | keep previous values | Y (pending) |
| 44 | charis_register_cube.pro:674 | WRONG | PSF method loops over `refsz[3]` instead of `sz[3]` | fix | Y (pending) |
| 45 | charis_register_cube.pro:528,733,918,1015 | DEAD | plot `yrange` uses `guessoffsets[0]`/`[1]` as lower/upper bound | cosmetic | N |
| 46 | charis_register_cube.pro:985-988,1115 | GEOM | slice 10 and 201 hard-coded in method 2 and checkquality | derive | Y (pending) |
| 47 | charis_register_cube.pro:1079 | WRONG | `/checkquality` only processes slice 0 (`for il=0L,0`) though comment says all | loop all | Y (pending) |
| 48 | charis_register_cube.pro:1126,1132 | DEAD | `set_plot,'ps'` never reset; `bsize` not a keyword | fix | Y (pending) |
| 49 | charis_register_cube.pro:308,671,978 | WRONG | PA wrap is a heuristic (`PA>180 and dec<lat`), comment "need to fix" | Python computes PA properly | N |
| 50 | charis_register_cube.pro:198-203,219-220,547-548 | DEAD | reads `spang*`,`cenmask`,`spmask`,`RA` unused; `xc0/yc0` not keywords; `cenxref` unused | document | N |
| 51 | charis_register_cube.pro (fork) `fitgauss` | CRASH | calls `calc_centroid_spots4`, defined nowhere in either tree | add wrapper on `mpfit2dpeak` | Y |
| 52 | charis_getcenter.pro:10, _pdi:103,117 | DEAD | write `im.fits`/`mask.fits` on every call | remove | Y (pending) |
| 53 | charis_getcenter.pro:62 / _pdi:167 | DEAD | `spang` undefined in dead `inter_spider` path; `median(cor,corlim)` uses 3 as filter width | remove path | Y (pending) |
| 54 | charis_getcenter_pdi.pro:154 | CRASH | writes `imz` (assignment commented out); tests `sy[i]` not `sy[j]` | fix | Y (pending) |
| 55 | charis_getcenter_pdi.pro:140 | WRONG | `im1[where(im1 lt 0)]=NaN` with no negatives → last pixel set NaN | count guard | Y (pending) |
| 56 | reg/charis_getcenter.pro vs subroutines/charis_getcenter.pro | DEAD | defined twice; path order decides which runs | delete `subroutines/` copy | Y (pending) |
| 57 | charis_check_psfquality.pro:54,90 | WRONG/GEOM | same slice-0-only loop and 201 as #46/#47; ladder list built then unused | fix | Y (pending) |

### 4. filt/

| # | File:line | Sev | What is wrong | Fix | IDL patch |
|---|---|---|---|---|---|
| 58 | charis_imrsub.pro:61 | WRONG | `n_elements(test) gt 0` always true (`file_search` returns `''`) → suffix always `reg_cal`; fails if specphot_cal not run | `count=` | Y (pending) |
| 59 | charis_imrsub.pro:80,119 | WRONG | `readcol` overwrites `filenum`; PA/rsat matched by row index, not file number → wrong PAs if `FNUM_SAT` edited after registration | match by number | Y (pending) |
| 60 | charis_imrsub.pro:133-135,159-160 | WRONG | empty `bad` → `[-1]` last pixel set NaN | count guard | Y (pending) |
| 61 | charis_imrsub.pro:248-250 | WRONG | `$` continuation only covers `angoffset`; `charis_northup` runs even with `/nonorthup`, with undefined `angoffset` | fix block | Y (pending) |
| 62 | charis_imrsub.pro:150 | DEAD | `lowpass` not a keyword → block dead | document | N |
| 63 | charis_imrsub.pro:97-98,199,218-219,231 | DEAD | `xc0/yc0` not keywords; `imt` re-read unused; `roi`,`good` unused; file written every iteration | tidy | Y (pending) |
| 64 | charis_imrsub.pro | DEAD | `imcubefinal_*.fits`, `ralignsdi.fits` written to CWD not `reduc/proc` | Python: `reduc/proc` | N |

### 5. specphotcal/

| # | File:line | Sev | What is wrong | Fix | IDL patch |
|---|---|---|---|---|---|
| 65 | charis_specphot_cal.pro:559, satspot_cal:176, specphot_cal_unsat:388 | WRONG | keyword declared `filtername=fname` but tested `keyword_set(filtername)` → user's filter override never applies | test `fname` | Y (pending) |
| 66 | charis_specphot_cal.pro:459 | CRASH | typo `cerro0` for `cerr0` (calmethod 5, ncal>1) | fix | Y |
| 67 | charis_specphot_cal.pro:719, satspot_cal:281, raw_contrast:311 | WRONG | `sat2noise` uses `cens[1,2,il]` (spot 2's x) instead of `cens[1,1,il]` (spot 1's y) → wrong noise for spot 1 | fix index | Y |
| 68 | charis_specphot_cal.pro:254,280 + get_charis_sat_contrast:165 | WRONG | calmethod 3 reads a per-channel attenuation vector but `get_charis_sat_contrast` uses only `manual0[0]` → per-channel ladder calibration lost | vector support | Y |
| 69 | charis_specphot_cal.pro:449-451 | WRONG | calmethod 5 reads every `ical` row from the last cube's header | per-cube | Y (pending) |
| 70 | charis_specphot_cal.pro:592-605 | WRONG | spot positions taken from the first cube's header for all cubes | per-cube | Y (pending) |
| 71 | charis_specphot_cal.pro:567-587 | CRASH | `case filtname0` has no `else` → unexpected CAL_BAND leaves `starmagval` undefined | `else: message` | Y (pending) |
| 72 | charis_specphot_cal.pro:115-120, _unsat:98, photometric_calibration:61-64 | WRONG | three different stellar-library numbering schemes (comment 0/1/2, prompt 1/2/3, code 1=Kurucz 2=Pickles 3=emp); `starlib=0` unsettable | Python: string enum; IDL: fix comments/prompt | Y (pending) |
| 73 | charis_specphot_cal.pro:235,244 | CRASH | compares numeric `modsize1` with `'unavailable'`; multiplies string by 1d3 if only one header is unavailable | type check | Y |
| 74 | charis_specphot_cal.pro:656-673 | GEOM | `fitbckgd` hard-codes 201 | derive | Y (pending) |
| 75 | charis_specphot_cal.pro:693-702 | WRONG | `aper` called with `[0,0]` in the `skyrad` slot (argument order) | fix order | Y (pending) |
| 76 | charis_specphot_cal.pro:797 | CRASH | `satflux[*,19]` hard-coded → fails for J band with `/verbose` | last channel | Y |
| 77 | charis_specphot_cal.pro:869-881 | WRONG | calmethod 0/3 assumes calibrated files == calibration cubes; with `/pick`/`datacube` count differs → index error | guard | Y (pending) |
| 78 | charis_specphot_cal.pro:917-938 | DEAD | `time_ratio` always 1 → `notimeratio`, `tint_calf` do nothing | document | N |
| 79 | charis_specphot_cal.pro:659-660, extract_1d_spectrum:171-275 | CRASH | call `get_xycind`, defined nowhere in repo (GPI routine) | add `subroutines/get_xycind.pro` | Y |
| 80 | charis_photometric_calibration_calculation.pro:76,109 | CRASH | `return,error(...)` — `error()` undefined | `message` | Y |
| 81 | charis_photometric_calibration_calculation.pro:23-39 | CRASH | `case filtname` no `else` (e.g. `'low'` crashes) | `else` | Y (pending) |
| 82 | charis_photometric_calibration_calculation.pro:178 | WRONG | smoothing kernel width `lambda[0]+dloglam` adds µm to a dimensionless log step → kernel far too narrow | `λ·ln(10)·dloglam` | Y (pending) |
| 83 | charis_photometric_calibration_calculation.pro:207 | WRONG | Vega convolved with a kernel built for the model's wavelength grid | per-grid kernel | Y (pending) |
| 84 | charis_photometric_calibration_calculation.pro:274-298 | WRONG | printed unit labels shifted by one (Jy labelled 'mJy', etc.); factors correct | fix labels | Y (pending) |
| 85 | charis_photometric_calibration_calculation.pro:180,332 | DEAD | `fwhm=0` → `Npixel=0`; dead `stop` after `return` | guard/remove | Y (pending) |
| 86 | charis_specphot_cal_unsat.pro:453 | CRASH | `cens` undefined with `/verbose` | fix | Y (pending) |
| 87 | charis_specphot_cal_unsat.pro:476 | WRONG | `if n_elements(ncal) gt 1` always false → `snratio_sats` stays 2-D, `CERR` indexed wrongly for several cubes | `ncal gt 1` | Y (pending) |
| 88 | charis_specphot_cal_unsat.pro:451 vs 606-607 | WRONG | header says sky annulus 2/6, code uses `[4,10]` | consistent | Y (pending) |
| 89 | charis_specphot_cal_unsat.pro:456,468,171-175 | DEAD | `setskyannulus` not a keyword; `eflux` with `phpadu=1` meaningless; `mjd_in` unused | document | N |
| 90 | charis_satspot_cal.pro:272-275 | WRONG | `dist_circle,mask1,dimcal[0],[x,y]` passes the centre as one array → ignored → spot masks centred on image centre | pass `x,y` | Y |
| 91 | charis_satspot_cal.pro:313-315,321,327,332 | GEOM | assumes ≥3 cubes, channels 0..21, `p0=[9]`/`[10]` "near 1.55 µm" (lowres only) | derive from wavelengths | Y (pending) |
| 92 | charis_satspot_cal.pro:333 | WRONG | `good` can be empty | guard | Y (pending) |
| 93 | charis_makeemppsf.pro:45, makesynthpsf:276 | WRONG | same `file_search` `n_elements gt 0` bug as #58 | `count=` | Y (pending) |
| 94 | charis_makeemppsf.pro:132 | WRONG | "missing SATS" test uses `(size(hdrval))[1]` (type code for scalars) → test never works | `size(/type)` | Y (pending) |
| 95 | charis_makeemppsf.pro:134,141 | WRONG | `goto` leaves previous spot's stamp in `psf_spot` (stale data) | reset | Y (pending) |
| 96 | charis_makeemppsf.pro:151 vs 90-91 | WRONG | `subarray_dist` sized `fsize` but stamp sized `dpsf`; mismatch if `psfsize≠fsize` | same size | Y (pending) |
| 97 | charis_makeemppsf.pro:158 | WRONG | empty `where` sets last pixel to 0 | guard | Y (pending) |
| 98 | charis_makeemppsf.pro:29,103-110 | DEAD | `outname` ignored; whole cube re-read per (λ,file); output lacks wavelength headers | honour `outname`; write wavelengths | Y (pending) |
| 99 | charis_makesynthpsf.pro:270,336,355,370 | WRONG | uses `pick` (not a keyword); "fwhm<1 means arcsec" guessed; empty `where` → `[-1]`; `outname` with path breaks `strsplit` | fix | Y (pending) |
| 100 | charis_extract_1d_spectrum.pro:95-112 | DEAD | code after `goto` with typo `extractube` and keywords specphot_cal doesn't accept → `pick`,`r_ap`,`calcube`,`fname`,`starname`,`mag`,`test`,`fluxunits` do nothing | remove; document | Y (pending) |
| 101 | charis_extract_1d_spectrum.pro:154,161 | WRONG | `fitrange` defaulted before tested → `~keyword_set(fitrange)` branch unreachable | reorder | Y (pending) |
| 102 | charis_extract_1d_spectrum.pro:236 | WRONG | `/breakout` only catches NaN but `cntrd/gcntrd` signal failure with −1 → never fires (also breaks `est_astrom_grid`) | test `lt 0` | Y (pending) |
| 103 | charis_extract_1d_spectrum.pro:388-391 | WRONG | `file_search` `''` → `ntest=1`, "0 files" case unreachable, `readcol,''`; >1 file prompts GUI (breaks autoreduce) | `count=`, explicit `throughput_file` | Y (pending) |
| 104 | charis_extract_1d_spectrum.pro:171-177,274-294,519,522-523 | GEOM | 201, 201/2, 0.0162 hard-coded | derive | Y (pending) |
| 105 | charis_extract_1d_spectrum.pro:519,471-507 | WRONG | `(x-100)*pixscale` labelled "E" (sign reversed for north-up/east-left); y-axis always "mJy" | fix labels | Y (pending) |
| 106 | charis_extract_1d_spectrum_disk.pro:45-59 | CRASH | live specphot_cal call with invalid keywords and `extractube` typo → crashes on uncalibrated cubes | remove/fix | Y (pending) |
| 107 | charis_extract_1d_spectrum_disk.pro:144-147 vs 166 | WRONG | aperture flux discarded, replaced by single pixel `ff[xposf,yposf]` (float indices truncated) | use aperture flux | Y (pending) |
| 108 | charis_extract_1d_spectrum_disk.pro:171 | WRONG | error divides by `pixscale`, flux by `pixscale^2` | consistent | Y (pending) |
| 109 | charis_extract_1d_spectrum_disk.pro:155,226,186 | WRONG | slice 11 hard-coded; `throughput` undefined with `/clipoutliers` alone; same `file_search` bug | fix | Y (pending) |
| 110 | charis_calc_spec_covar.pro:429 + gbfunc.pro | CRASH | inner `GBFUNC` uses undefined `rho` and shadows `gbfunc.pro` once compiled | delete inner copy | Y (pending) |
| 111 | charis_calc_spec_covar.pro:284 | WRONG | `1-(a+b)+a+b` always 1 (should be `A_δ+A_ρ+A_λ`) | fix | Y (pending) |
| 112 | charis_calc_spec_covar.pro:199-205,263 | WRONG | spline array has `num+1` elements vs `num`; `where` combines different lengths (off-by-one) | fix | Y (pending) |
| 113 | charis_calc_spec_covar.pro:209 | WRONG | bins `±step` wide overlap | half-width | Y (pending) |
| 114 | charis_calc_spec_covar.pro:324-330,361-367,299,336 | GEOM | J/H/K channel ranges `[0:8]`,`[9:14]`,`[16:20]` and `15*lam_ratio` hard-coded (lowres only) | derive from wavelengths / `rho` | Y (pending) |
| 115 | charis_calc_spec_covar.pro:187 | WRONG | sorted index applied to unsorted `mod_speccor` (result unused) | remove | Y (pending) |
| 116 | charis_est_bbphot.pro:109-131,213 | CRASH | only units 0 and 3 handled; units 1,2,4 (documented) crash with `c` undefined | implement | Y (pending) |
| 117 | charis_est_bbphot.pro:224 vs get_phot_zeropoint | WRONG | zero points `[1560,1040,670]` Jy disagree with `get_phot_zeropoint` (H=1050, K=667) | single table | Y (pending) |
| 118 | charis_est_bbphot.pro:54,71,89,91,141,99 | GEOM/DEAD | channel pairs hard-coded; errors added linearly; `diagnostic` not a keyword | derive; quadrature | Y (pending) |
| 119 | charis_est_astrom_grid.pro:253-256,292 | WRONG | defaults `xrange=[99,102]`,`yrange=[61,64]`,`fitrange=[0:21]` are one target's values; `noexten` ignores `nothroughputcor`; relies on broken `/breakout` (#102) | require args | Y (pending) |
| 120 | charis_get_phot_zeropoint.pro | CRASH | `case` no `else` | `else` | Y (pending) |

### 6. cals/ and setup/

| # | File:line | Sev | What is wrong | Fix | IDL patch |
|---|---|---|---|---|---|
| 121 | get_charis_wvlh.pro:26,39-57 | CRASH | `strtrim` without flag 2; unknown band (`'low'`,`'broadband'` — which specphot_cal accepts) leaves `wavecalfile` undefined | normalise names, `else` | Y (pending) |
| 122 | get_charis_sat_contrast.pro:101,165,168,175 | WRONG | modifies caller's `modulation_size`; only `manual0[0]` used (#68); manual error always 0.05; `/help` falls off end without `return` | fix | Y |
| 123 | charis_path.pro:13,21 | WRONG | hard-coded `~/idl_tools/ADI_dl/charis_dpp_new/` and `/Users/thaynecurrie/.../btsettl/`; `case` no `else` | env var / `else` | Y |
| 124 | charis_get_constant.pro | DEAD | no `else`; local `Dtel`,`pixscale` unused | `else` | Y |
| 125 | charis_test.pro:9, charis_which.pro:74,102 | CRASH | `output` unset when found in `.`; `fetch=level` undefined | fix | Y (pending) |

### 7. psfsub/ (adialoci, sdialoci, adiklip, rdiklip/_2, ssdi)

| # | File:line | Sev | What is wrong | Fix | IDL patch |
|---|---|---|---|---|---|
| 126 | adialoci:148-158, sdialoci:139, adiklip:115, rdiklip | WRONG | `file_search` `''` → `suffname` always `reg_cal` (#58) | `count=` | Y (pending) |
| 127 | adialoci:328, adiklip:264, rdiklip | WRONG | `if n_elements(drsub0 eq 1)` misplaced paren → always true; variable-width `drsub` branch unreachable, 4-element `drsub0` goes into scalar math | `n_elements(drsub0) eq 1` | Y |
| 128 | adialoci:1006, sdialoci:980, adiklip:947, rdiklip:959 | WRONG | `if ~keyword_set(mean)` — `mean` never defined; keyword is `meanadd` → combination always median | `meanadd` | Y |
| 129 | adialoci:765-793 | WRONG | with `usecoeff`, loaded `c` overwritten by a fresh inversion (`goto,skipmatrixinversion` commented out) → `usecoeff` has no effect; `c_use*` never reset on no-match (stale) | restore goto; reset | Y (pending) |
| 130 | adialoci:827-831 | WRONG | fwdmod `refpert` filled at `refpert[igood]` where `igood` indexes `isub`, not `iopt` → self-subtraction term on wrong pixel set | consistent index set | Y (pending) |
| 131 | adialoci:838 vs 658 | WRONG | `pert1` built from `annuli` only while covariance `aa` includes model → inconsistent system | consistent | Y (pending) |
| 132 | adialoci:703-710 | WRONG | 7σ outlier mask `igoodf` applied only when `/postsdi`; discarded otherwise | apply always | Y (pending) |
| 133 | adialoci:746-747 | CRASH | `nref=0` (ALOCINRE missing) → `[0:-1]` | guard | Y |
| 134 | adialoci:240-241 | WRONG | ALOCISVD missing → `svd=10^0=1` truncates nearly everything | sane default | Y (pending) |
| 135 | adialoci:971,993-997 | WRONG | `channel=k` copies the single-channel image into every slice of output | fill only k | Y (pending) |
| 136 | adialoci:719,879 | WRONG | `append=(ir+it gt 0)` assumes zone (0,0) always written; skipped zone → first write appends to stale file | Python: no scratch files | N |
| 137 | adialoci:437-443,526-552,1017,1022-1039,1046 | DEAD | `norm=0` block, `if 0` block, unreachable `outfile` branch, re-read/re-write same file, `free_lun,1` unguarded | tidy | Y (pending) |
| 138 | sdialoci:830,840,851 | WRONG | fwdmod `refpert=fltarr(n_elements(ioptf))` indexed with `ioptf` (indices into `ia`, up to `n_elements(ia)`) → IDL clips → SDI forward model silently corrupted | `annuli_fwd[ioptf,indim[k]]` | Y (pending) |
| 139 | sdialoci:704,719 | WRONG | `replicate(1,n_elements(iopt))` sized by `iopt` but `z` has `n_elements(ioptf)` rows → silent truncation; row mask applied to wrong array | size by `ioptf` | Y (pending) |
| 140 | sdialoci:633,638 | WRONG | empty-sector guard commented out → `iopt=-1` → `iaopt[-1]` last element | restore guard | Y (pending) |
| 141 | sdialoci:891,922 | WRONG | per-channel temp files never deleted; skipped zone appends to previous frame's file | Python: in-memory | N |
| 142 | sdialoci:748-753,209-211,309 | DEAD | `nref` correlation selection commented out → `nref` ignored; `na` no default; only scalar `drsub` | document | N |
| 143 | sdialoci:362-415 | DEAD | ring index files recomputed for every frame and channel (invariant) | Python: once | N |
| 144 | adiklip:539 vs 584 | WRONG | `indices_images.dat` written before `if cgood lt 1 then continue` → skipped zone misaligns every later zone's pixels; temp files never deleted | write after; Python: masks | Y (pending) |
| 145 | adiklip:710,770 (+rdiklip_2) | WRONG | `ref` filled only at `igoodf` pixels → outlier pixels left as raw data instead of NaN | set NaN | Y (pending) |
| 146 | adiklip:660-662 | WRONG | `eigenvect_full` reversed only `if npcaf gt 1` but eigenvalues always reversed → `npca=1` fwdmod indirect term mismatched | reverse always | Y (pending) |
| 147 | adiklip:735 | WRONG | sum over modes divides by `(λ_k−λ_j)` and takes `√λ_j`; rank-deficient → Inf/NaN | guard | Y (pending) |
| 148 | adiklip:371,376 | WRONG | `writeu,im[ia]` without `float()` but read back as float → BITPIX −64 cubes misread | `float()` | Y (pending) |
| 149 | adiklip:635 | WRONG | `usecoeff` reforms with current `c1×npcaf`; different `c1` crashes; stale reuse | guard | Y (pending) |
| 150 | rdiklip:741-743 | WRONG | `savecoeff` format `(i3,i3,i3,...)` no separators → fields merge at `nf≥100`, `usecoeff` misparses | add `1x` | Y (pending) |
| 151 | rdiklip:165 | WRONG | without `refpath` opens `dialog_pickfile` | require `refpath` | Y (pending) |
| 152 | rdiklip:692,722 | WRONG | `npca=npca<c1` overwrites caller's keyword; `npca=1` eigenvector transposed (inconsistent with adiklip) | local var | Y (pending) |
| 153 | rdiklip:437-508,798 | DEAD | copy-pasted zone geometry with undefined `geom`; `goto,skipme` dead block | tidy | Y (pending) |
| 154 | ssdi:118 | CRASH | `medfitsme` undefined | `charis_combinefits` | Y |
| 155 | ssdi:120-125 | WRONG | `/northup` result `imcomb0` never written; `ssdi.fits` gets un-rotated image | write it | Y (pending) |
| 156 | ssdi:129 | GEOM | channels 0..21 hard-coded | derive | Y (pending) |
| 157 | ssdi:89 | WRONG | reference median includes slice `il` itself | exclude | Y (pending) |
| 158 | ssdi:7 | DEAD | help documents `pca=` but keyword is `klip` | fix | Y (pending) |
| 159 | all four psfsub | DEAD | help text versions wrong; history records `theta=-allpa` though rotation is `charis_northup`+`angoffset`; `spang`,`rsat`,`northpa`,`lng`,`el`,`az`,`noise_im`,`dtmean` unused | document | N |

### 8. partialsub/ (fwdmod, fakeplanet, attenmap, negspec, disk)

| # | File:line | Sev | What is wrong | Fix | IDL patch |
|---|---|---|---|---|---|
| 160 | klip_fakeplanet:60 | CRASH | `goto,skiptotheend` but no such label; `help` not declared → file does not compile | add label/keyword | Y |
| 161 | aloci_fakeplanet:643, klip_fakeplanet:520 | CRASH | `profrad_tc,imslice,1,1,rmax` with `rmax` undefined | read from header | Y |
| 162 | aloci_fwdmod_planet:101-102 | CRASH | `adipsdi`,`sdipadi` used but not declared keywords (help advertises them) | declare | Y |
| 163 | aloci_fwdmod_planet:589,606-608 (+attenmap) | WRONG | chain 2/3 intermediate psfsub calls have no `outfile` → overwrite user's `reduc/proc/final_fc.fits` | unique outfiles | Y (pending) |
| 164 | aloci_fwdmod_planet:180-246 | WRONG | `case method` mixes string/int; help says targfile order "log(contrast),x,y" but `readcol` reads x,y,contrast | fix help, case | Y (pending) |
| 165 | aloci/klip_fwdmod, fakeplanet, attenmap, negspec | WRONG | PSF normalisation differs (aperrad vs 0.5·fwhm) across the family; moot because `insert_planet_into_cube` renormalises to `R_AP` | Python: one normalisation | N |
| 166 | klip_fwdmod:520, disk_fwdmod:599, negspec_disk:314 | WRONG | `where(finite eq 0)` → `[-1]` last pixel zeroed | guard | Y (pending) |
| 167 | klip_fwdmod:203,397-400 | DEAD | `asdi`/`sdi_reducname` not keywords (block dead); `cv_coord` on arrays gives wrong shape for `nplanet>1` | fix | Y (pending) |
| 168 | aloci_attenmap:615, klip_attenmap:482 | WRONG | `openw,5` inside the `irc` loop → each `atten*.dat` keeps only the last planet | move out | Y (pending) |
| 169 | aloci_attenmap:625, klip_attenmap:492 | WRONG | `throughput_map_lambda[xp0,yp0,il]=res_flux[ilambda,irc]` drops the `itc` index (always offset 0); float subscripts truncated | add `itc`, `round()` | Y (pending) |
| 170 | klip_attenmap:431 | CRASH | passes `refpath=refpath` but `refpath` not a keyword → RDI always falls to GUI | declare | Y (pending) |
| 171 | aloci_attenmap:84,470,140/413 | DEAD | `sdi`,`nonorthup` undeclared; `ntc` defaulted twice | fix | Y (pending) |
| 172 | aloci_negspec_planet:468,537-545; klip:387,403-411 | WRONG | output written as `res0_np_.fits` but code reads `res0.fits`/`res0_collapsed.fits` → stale file from an earlier fwdmod run | read `_np_` | Y (pending) |
| 173 | aloci_negspec_planet:490,501; negspec_disk:378,389 | WRONG | chains 2/3 call sdialoci with `/fwdmod` (copy-paste) → negative-planet data treated as a model | drop `/fwdmod` | Y (pending) |
| 174 | aloci_negspec_planet:149, klip:123 | CRASH | `case method` no default/else → unset `method` is "undefined", `method=2` "no match" | `else` | Y (pending) |
| 175 | klip_negspec_planet:279 | CRASH | calls `myaper` (undefined; should be `charis_myaper`) | fix | Y |
| 176 | klip_negspec_planet:340,356 | WRONG | `nonorthup` not declared; writes `outputcube.fits` to CWD every frame | fix | Y (pending) |
| 177 | aloci_disk_fwdmod:28,66,107-110 | CRASH | `reductype,sdi,adipsdi,sdipadi,prefname,help` used but not declared → only ADI path can run; help example uses non-keywords | declare | Y (pending) |
| 178 | aloci_disk_fwdmod:162 vs 205 | WRONG | cube read from `reducdirorg` but header from `reducdir` (differ when picked via GUI) | same dir | Y (pending) |
| 179 | aloci_disk_fwdmod:364,376 | WRONG | `snrlim` presence tested but value never passed to `charis_generate_roi` | pass it | Y (pending) |
| 180 | klip_disk_fwdmod:548 | WRONG | `sxpar(header,'zrsub')` — key is `'rsub'` → `zrsub` always 0 | fix key | Y (pending) |
| 181 | disk_fwdmod:325-329 | WRONG | `maskrad` defaults to 0.1 so the "if you need to mask" block always runs | default 0 | Y (pending) |
| 182 | disk_fwdmod:282-291 | GEOM | Gaussian PSF `npixel=201` regardless of `dimx` | derive | Y (pending) |
| 183 | disk_fwdmod:747-748,793-794 | WRONG | FITS keys `'csq_col (binned image)'`, `'npts (binned image)'` invalid; written to `h1col` instead of `h1cube` | short keys, right header | Y (pending) |
| 184 | disk_fwdmod:827,101 | WRONG | format `f6.3` prints `******` for KSI0>99.999; header has 11 columns, format writes 10 | widen | Y (pending) |
| 185 | disk_fwdmod:710 | WRONG | χ² weight `1/sigdata²` with zeros → Inf | guard | Y (pending) |
| 186 | disk_fwdmod:587-603 | DEAD | same rotated model added `nsc=21` times /21 (no shift) — pure waste | once | Y (pending) |
| 187 | disk_fwdmod vs negspec_disk vs planet codes | WRONG | three different model-rotation conventions (`angoffset` / final `ROTANG` / `+ROTANG`) | Python: one convention (§2.3) | N |
| 188 | negspec_disk:356 | CRASH | `charis_adialoci,...,/np` — no `np` keyword; also no `/fc` | `/fc` | Y |
| 189 | negspec_disk:313 vs 320 | WRONG | `bad` reused: zero pixels of data in first iteration, NaN pixels of `scratch` later | separate vars | Y (pending) |
| 190 | negspec_disk:312,413-463 | DEAD | writes `psf.fits` per frame/channel/isc; unreachable block uses undefined vars | tidy | Y (pending) |
| 191 | recalc_disk_fwdmod_chisq:4 | CRASH | `N_PARAMS() eq 0` always true (keyword-only procedure) → always prints help and exits | `keyword_set(help)` | Y |
| 192 | recalc_disk_fwdmod_chisq:85 | CRASH | calls `snratio_sub` (undefined) | `charis_snratio_sub` | Y |
| 193 | all partialsub | DEAD | results (`synth_throughput*.txt`, `attenmap*.fits`, `psf.fits`, …) written to CWD; `plot` needs X; `findfile` deprecated | Python: `reduc/proc`, `savefig` | N |

### 9. planetspec/

| # | File:line | Sev | What is wrong | Fix | IDL patch |
|---|---|---|---|---|---|
| 194 | charis_insert_planet_into_cube.pro:221 | CRASH | no bounds check; planet near edge clips/wraps via `ipsf` indexing | clip | Y (pending) |
| 195 | charis_insert_planet_into_cube.pro:182-183,208 | WRONG | `parangle` wraps at ±180° near transit → `dpa` blows up; smear sign opposite to callers (harmless only because symmetric) | unwrap | Y (pending) |
| 196 | charis_insert_planet_into_cube.pro:141,66 | CRASH | with `empspectrum`, `diagnostic` uses undefined `starspec`; pick mode overwrites positional `inputpsfcube` with a filename | guard | Y (pending) |
| 197 | charis_planet_photometric_calibration_calculation.pro:80-95 | CRASH | extension test `strsplit(modelspec,'.')[1]` crashes on no-dot names, misclassifies paths with `./` | use `file_basename`/suffix | Y (pending) |
| 198 | charis_planet_photometric_calibration_calculation.pro:141-142 | WRONG | `where(lt 0)` → −1 zeroes last channel | guard | Y (pending) |
| 199 | charis_planet_photometric_calibration_calculation.pro:36,84,94,127-130 | CRASH | `case filtname` no `else`; comment says Å but code needs µm; `stop` inside library function | `else`; fix comment; `message` | Y (pending) |
| 200 | charis_convert_planetspec.pro:6 | WRONG | tests `keyword_set(filtname)` but variable is `filtname0` → filter always reset to `'lowres'` | fix | Y (pending) |
| 201 | charis_convert_planetspec_to_lowres_calculation.pro:34 | WRONG | 3-row test `sz[1] eq 3` should be `sz[0]`; `readcol` asks 3 columns on 2-column files | fix | Y (pending) |
| 202 | …_calculation, …_calculation2, planetphot | WRONG | positional `lambda` ignored (overwritten by `get_charis_wvlh`); `PSF_GAUSSIAN(Npixel=0)` when model doesn't cover grid edges | honour arg; guard | Y (pending) |
| 203 | …_calculation:100 | WRONG | `goto skipoverme` makes mJy conversion dead → returns F_λ though comments say mJy | fix | Y (pending) |
| 204 | convert vs photcal | WRONG | inconsistent FITS layouts (`[2,N]` vs `[N,2]`) | Python: one reader | N |
| 205 | planetphot:261 | DEAD | `close,/all` after `return`; `stop` on unknown filter | fix | Y (pending) |

### 10. pdi/ (old tree + fork rewrite)

| # | File:line | Sev | What is wrong | Fix | IDL patch |
|---|---|---|---|---|---|
| 206 | charis_pdi.pro:105,143-152 (fork: `qbad`,`ubad`,`outside_fov`) | WRONG | `where()` without count; default FOV = full frame → `outside_fov=-1` → last pixel of every slice set NaN | guards | Y (pending) |
| 207 | charis_pdi.pro:114-118,244 | WRONG | bad cycle triggers `break` but processing continues; later reads `Qout[ci]` never written | abort/skip cycle | Y (pending) |
| 208 | charis_pdi.pro:121-124 | WRONG | missing HWP position → `where` −1 → `cycle_inds[-1]` wrong frame, no error | error | Y (pending) |
| 209 | charis_pdi.pro:64 | WRONG | assumes `reduc_left.log` row order matches file list; indexes `allpa` with HWP-table indices | match by number | Y (pending) |
| 210 | charis_pdi.pro:289,291 | CRASH | `sxaddpar` with `meanadd`/`first_order_ip_correct` undefined when flag unset (verify in fork) | `keyword_set()` | Y (pending) |
| 211 | charis_pdi.pro:284-293 | WRONG | FITS keys >8 chars (`FIRST_ORDER_IP_CORRECT`, `FOV_MASK_PARS`, `HWPINFO_NAME`) truncated/HIERARCH | short keys | Y (pending) |
| 212 | charis_pdi.pro:341-344 (old), fork `telluric_channels=[5,6,14,21]` | GEOM | 22-element `lam_offsets`/telluric list hard-coded | derive from wavelengths | Y (pending) |
| 213 | charis_pdi.pro:182,337 | GEOM | `dist(sz[0])`, `angarr(dimx)` assume square | use both dims | Y (pending) |
| 214 | charis_pdi.pro (fork) | DEAD | `coverage_map.fits` to CWD; `spawn rm` | `file_delete`; `reduc/proc` | Y (pending) |
| 215 | charis_pdi_flatfield.pro:71 | WRONG | `flatslice[where(flatslice eq 0)]=1` no count check | guard | Y |
| 216 | charis_pdi_flatfield.pro:65,86,53,91 | GEOM/WRONG | FOV ROI `[125,125,-26.9166]` hard-coded; suffix `'e_flat'` ignores `suffname`; unclear error with no dark flats; no guard on near-zero flat before division | fix | Y (pending) |
| 217 | charis_pdi_split_pols.pro:47-55 (old) / fork register | WRONG | `n_elements(test) gt 1` — single file treated as none | `count=` | Y (pending) |
| 218 | charis_pdi_split_pols.pro:95 | CRASH | divides by `dy0_pol` (zero possible) | guard | Y (pending) |
| 219 | charis_pdi_register_cube.pro:70-80 | WRONG | `spawn mv/cp` Unix-only; help lists `guessoffsets=` not a keyword; forwards `help=help` | `file_move`/`file_copy` | Y (pending) |
| 220 | charis_pdi_specphot_cal.pro:72-93 | WRONG | `where` without count; default FOV → last pixel NaN | guards | Y (pending) |
| 221 | charis_pdi_specphot_cal.pro:63,160-161 | DEAD | temp prefix `'n'` hard-coded; `spawn rm`; `modamp`/`notimeratio` silently ignored with `/unsat`; help lists non-keywords | fix | Y (pending) |
| 222 | charis_pdi_hwp_match.pro:124-125 | WRONG | `utchwp == utcexp` exactly → `hwpind` stale/undefined; `argmin=0` with later move → `hwpind=-1` → `angs[-1]` | handle ties/−1 | Y (pending) |
| 223 | charis_pdi_hwp_match.pro:89-100,85 | WRONG | `time_format` from last log row only; midnight UT rollover unhandled | parse per row; day wrap | Y (pending) |
| 224 | charis_pdi_hwp_match.pro:110,149,156 | WRONG | comments contradict code on start-vs-end exposure time | expose option | Y (pending) |
| 225 | charis_pdi_hwp_match.pro:80 | WRONG | `uni_angs` from all log entries → >4 angles gives `hwp_pos>3`, which `charis_pdi` assumes 0-3 | restrict to cycle angles | Y (pending) |
| 226 | charis_pdi_hwp_match.pro:234,240,174 | DEAD | `((360+a-b)-360) mod 360` wrap does nothing; `indgen` vs `lindgen` | tidy | Y (pending) |
| 227 | charis_pdi_rlims_gridsearch.pro | WRONG | pfname hard-coded to `V2254Ori_low.info` | parameter | Y (pending) |

### 11. tools/

| # | File:line | Sev | What is wrong | Fix | IDL patch |
|---|---|---|---|---|---|
| 228 | contrast_calc/charis_calc_final_contrast.pro:117,219 | CRASH | typo `numchannnels` → crashes with `/meanadd` | fix | Y |
| 229 | charis_calc_final_contrast.pro:123 | CRASH | `bbimage` undefined → `databb` keyword never works | fix | Y (pending) |
| 230 | charis_calc_final_contrast.pro:177,191,235 | WRONG | checks `maskcoords` but keyword is `maskcoord` → mask never applied | fix | Y |
| 231 | charis_calc_final_contrast.pro:181,197-203,236 | DEAD | `nofinite`,`nozero`,`nofilt` never used (always on) | wire up | Y (pending) |
| 232 | charis_calc_final_contrast.pro:108-110,138-155,211-213,308,169,73-74 | GEOM | lowres channel ranges; legend assumes 22; square frame; attenmaps read from primary HDU only | derive | Y (pending) |
| 233 | charis_calc_raw_contrast.pro:294 | WRONG | `mask le aperrad` compares whole 22-element array → only first 22 pixels masked | `aperrad[il]` | Y (pending) |
| 234 | charis_calc_raw_contrast.pro:311 | WRONG | spot-2 noise index (#67) | fix | Y |
| 235 | charis_calc_raw_contrast.pro:483-487 | CRASH | `profplotname` vs `prefplotname` → crashes unless `outputfile` set | fix | Y |
| 236 | charis_calc_raw_contrast.pro:281 | WRONG | `sat*noise` used before computed at :310-313 → error underestimated | reorder | Y (pending) |
| 237 | charis_calc_raw_contrast.pro:11,203,459,524,254 | DEAD | `help`,`verbose`,`nosub` not keywords; `suffname` undefined with `/pick` | declare | Y (pending) |
| 238 | charis_calc_raw_contrast.pro:301,357,401,434,394-396,525 | WRONG/GEOM | non-sub path uses `abs(imcal)` not masked slice; `imcol.fits` to CWD; output names contain `.fits_`; channels 4/11/18 hard-coded | fix | Y (pending) |
| 239 | detection/charis_check_for_detection.pro:68-79 | WRONG | no `LOCI_*`/`KLIP_*` keys → rmin=rmax=0 → never detects anything (classical ADI/RDI) | fallback | Y (pending) |
| 240 | charis_check_for_detection.pro:101-153,45,179 | GEOM/WRONG | 201 hard-coded; `strsplit(datacube,'.')` breaks on dotted paths; `good` not checked for −1 | fix | Y (pending) |
| 241 | astrometry/charis_est_prop_mot.pro:12,14 | CRASH | calls `hci_querysimbadf`/`hci_queryvizierf` (undefined; repo has `charis_*`) | rename | Y |
| 242 | charis_est_prop_mot.pro:123,215-216 | WRONG | ecliptic offsets (`euler,…,3`) plotted as E/N equatorial with no frame rotation; RA proper motion added to ecliptic longitude → physics wrong except near pole/equator | Python: astropy `apply_space_motion` | N |
| 243 | charis_est_prop_mot.pro:299-302 | CRASH | runs unconditionally using `predastrom`,`sigepos` → crashes without `posfile` or with <2 epochs | guard | Y (pending) |
| 244 | charis_est_prop_mot.pro:24,36,43,195,156,293-294 | WRONG/DEAD | `plxset gt 0` (#1); typo `maxoffset`; O(N²) unused array; hard-coded "Oct 2015" point; `aberration` unused | fix | Y |
| 245 | atmosfit/charis_atmosmodelfit.pro:51-54,70-86 | CRASH | `filtername0` undefined if no wavelength case matches; `case band` no `else`, case-sensitive | `else` | Y (pending) |
| 246 | charis_atmosmodelfit.pro:59,109-110 | GEOM | masked channels 5/6/14/15 and loops 0..21 hard-coded | derive | Y (pending) |
| 247 | charis_atmosmodelfit.pro:177,196,211,15/100,283-292 | DEAD/WRONG | `metals` overwritten; fixed LUN 2; `cut` not checked for −1; covar read twice; `stop` in skipped block; `sigfig`/`rplot` missing (dead path) | tidy | Y (pending) |
| 248 | atmosfit/arch/charis_atmosmodelfittest.pro:156 | CRASH | bare URL line → won't compile | comment out | Y |
| 249 | empcomp/charis_get_numspt.pro:19-51 | WRONG | letter checks in order K,M,L,T match any letter in string; `strsplit(sg,sspectype)` takes last token → `'L1 VL-G'` parses as 0 | regex parse | Y (pending) |
| 250 | empcomp/charis_calcsptgravind.pro:38,53,55,13,90 | WRONG | index definitions reuse numerator channel in denominator (contradicts comments); `c0` hand-changed to 12.19; `mode` set → nothing computed then crash | fix per comments | Y (pending) |
| 251 | empcomp/charis_empbdplanspec.pro:189 | WRONG | `fname.contains(name)` substring match → false matches, last wins | exact match | Y (pending) |
| 252 | charis_empbdplanspec.pro:161,124-133,32,97,242,307 | WRONG/GEOM | always reads Montreal `.tsv` so Bonnefoy fits get empty spt/gravity; unknown library leaves `librarydir` undefined; masked channels differ from atmosfit; 0..21 loop; `modeldir='/'`; dof `ngood-1` | fix | Y (pending) |
| 253 | empcomp/plotchisqdist_charis.pro:43-53 | WRONG | second sample reuses `sample`/`filein` (same file read twice); read into `h202` overwriting first sample's H2O-2 index | separate vars | Y (pending) |
| 254 | plotchisqdist_charis.pro:170,172 | CRASH | `badarray` undefined (definitions commented) → crashes with `/maskbad` | define | Y (pending) |
| 255 | plotchisqdist_charis.pro:17 | WRONG | dir `Montreal_spectral_library` vs `Montreal_Spectral_library` written by empbdplanspec (case matters on Linux) | consistent | Y (pending) |
| 256 | plotchisqdist_charis.pro:269,368,397,460 | DEAD | hard-coded data tweak `numerical_spt[good[1]]=10.25`; κ And b values; "top 10" loop over wrong size | remove tweaks | Y (pending) |
| 257 | empcomp/charis_est_lum_bd.pro:28 | WRONG | `spt` required even when `bck` given; results only printed | return value | Y (pending) |
| 258 | filters/charis_get_phot_bandwidth.pro:49 | DEAD | `stop` on unknown filter | `message` | Y (pending) |
| 259 | imagequality/charis_calc_ee.pro:90,3-13,16,18 | CRASH | calls `profradxy` (missing); help always prints; `rbk` unused; `il=0` unsettable | use `profrad_tc` | Y |
| 260 | imagequality/charis_planefit_ex.pro | CRASH | calls missing `get_xycind` (#79); 201 and [100,62] hard-coded; writes `blah*.fits` | example only; fix via #79 | Y (pending) |
| 261 | psf/charis_resamplepsf.pro:33,39,40,35 | DEAD | `typ` undefined; writes `tmp.fits`; no-op `subarr`; larger `scalefac` shrinks (name misleading) | tidy; document | Y |
| 262 | diskanalysis/other/charis_diskmodel_resample.pro:5 | CRASH | extra `)` → syntax error | fix | Y |
| 263 | charis_diskmodel_resample.pro:50-54,18 | WRONG | help path reaches `end` without `return`; half-size comparison crops 401-px input | fix | Y (pending) |

### 12. subroutines/ (top level + chisq_calcs)

| # | File:line | Sev | What is wrong | Fix | IDL patch |
|---|---|---|---|---|---|
| 264 | profrad_tc.pro:88 | WRONG | `interpol(rayont,p1dt,…)` value/abscissa arguments swapped (`/nan` p2d path) | swap | Y |
| 265 | profang_tc.pro:21-22,35 | WRONG | radii 10-35 hard-coded ignoring `r1`/`r2`; loop stops at n−2 so last azimuth bin is 0 | honour args; loop to n−1 | Y (pending) |
| 266 | mkspider.pro:54-55 | WRONG | duplicated `+=` → horizontal spider vanes added twice | remove dup | Y |
| 267 | fixpix_rs.pro:109 | WRONG | `xs/5*ys/5` evaluates as `((xs/5)*ys)/5` → `sigma` longer than filled length; trailing zeros bias `medsig`/`sigsig` | parenthesise | Y |
| 268 | fixpix_rs.pro:208-216,173/203 | WRONG | mask sign +1/−1 inconsistent between branches; `mparm` inconsistent | unify | Y (pending) |
| 269 | charis_northup.pro:26,50 | WRONG | `hrot2` rotates the header, then `rotat(…,hdr=)` rotates it again by `angoffset` → CD/CRPIX off by `angoffset`; centre `sz/2` vs `(n−1)/2` half-pixel for even sizes | rotate header once | Y (pending) |
| 270 | rotat_cube.pro:26-27 | WRONG | `yi=igood/dimx` not taken `mod dimy` for 3-D indices → slice 0 alone sets NaN footprint | `mod dimy` | Y |
| 271 | rotat.pro | WRONG | when `missing` not set, edges extrapolated | Python: NaN outside | N |
| 272 | shift_sub.pro | WRONG | integer shifts use `shift()` (wrap-around) while fractional shifts clamp at edges → inconsistent edge behaviour | Python: NaN pad both | N |
| 273 | charis_speckle_alignsub.pro:90 | CRASH | typo `cubin` → `/fft` path crashes | fix | Y |
| 274 | charis_speckle_alignsub.pro:44,32 | GEOM | centre fixed at [100,100]; modifies `locs` in place | use PSFCENT; copy | Y (pending) |
| 275 | charis_alignspeckle.pro:55,61 | DEAD | `stop` when a spot is NaN; `badvalue` unused | `message` | Y (pending) |
| 276 | charis_getrelshift.pro:9-10 | CRASH | `round(sxg)` runs even when `sxg` not passed → default lost, crash | guard | Y |
| 277 | charis_snratio_sub.pro:66,69,71,23 | DEAD/GEOM | writes `roimask.fits`/`immask.fits` to CWD; 201 hard-coded; mutates caller's `fwhm` | fix | Y (pending) |
| 278 | snratio.pro:44-46 | WRONG | `expfilt` overwritten unconditionally; no callers | document | N |
| 279 | charis_calc_chisqbin.pro:11 | WRONG | `congrid` does nearest-neighbour sampling, not binning as the comment says | Python: block-mean | N |
| 280 | charis_combinefits.pro:67 | WRONG | `1/i1` → Inf where no valid pixels; assumes square | guard | Y (pending) |
| 281 | charis_generate_roi.pro:80,86 | WRONG | slice `dx+1` wide for even `dx`; `if … then A & B` ambiguous; `bad=-1` last pixel; default size 201 | fix | Y (pending) |
| 282 | charis_myaper.pro:55 | DEAD | `findgen(1.,ny)` works by accident; no sky subtraction | document | N |
| 283 | sumaper_im.pro:59,39 | WRONG | no bounds check → rows wrap; `rmax` undefined if not given (comment says default 100) | guard; default | Y (pending) |
| 284 | svd_invert.pro | WRONG | divide by zero if all singular values 0; `vecWKernel` undefined when kernel empty | guard | Y (pending) |
| 285 | nbr2txt.pro:21, nbrlist.pro:43,57 | WRONG | `alog10(0)` breaks file number 0; `fix()` overflows above 32767 | `long()` | Y |
| 286 | param.pro:47 | WRONG | rewrites the `.info` file on every `/get` (and creates it if missing) | write only on set/del | Y (pending) |
| 287 | hrot2.pro:124-125,135-152 | DEAD | 2-D check always passes; interactive `read` prompts if <8 args | tidy | Y (pending) |
| 288 | angarr.pro | WRONG | centre is `n/2`, not `(n−1)/2` (inconsistent with `dist_circle`) | Python: `image_center` | N |
| 289 | rot_ratef.pro | WRONG | `acos` azimuth has quadrant ambiguity | `atan2` | Y (pending) |
| 290 | intersect.pro | WRONG | counts duplicates within one array as intersections | Python: `np.intersect1d` | N |
| 291 | fftscale.pro, inter_spider.pro | CRASH | depend on missing `smartshift`,`mtfpix`,`lineint`; only dead callers | document; not ported | N |
| 292 | charis_queryvizierf.pro:194 | CRASH | calls missing `strcompress2` (constraint path) | `strcompress` | Y |
| 293 | chisq_calcs/chisqcalcwbin.pro:28 | WRONG | `widthobs(diff/eobs)` is a subscript, not a product | `*` | Y |
| 294 | chisqcalcwbin.pro:35 | CRASH | normalises by `widthphot` in the branch where photometry is *not* given → undefined | fix | Y (pending) |
| 295 | chisqcalcwcovar.pro:36-46,51 | WRONG | excluded channels zeroed in `diff` but full covariance kept (not equivalent to dropping channels); matrix inverted 2e4 times in loop | Python: sub-block + analytic α | N |
| 296 | chisqcalc.pro | WRONG | `ndof` ignores photometry points and the fitted parameter | fix | Y |
| 297 | chisqcalcwcovarwbin.pro:65 | WRONG | normalisation includes excluded channels; `phot`,`widthphot` unused; no callers | document | N |

### 13. scripts/ and cross-cutting

| # | File:line | Sev | What is wrong | Fix | IDL patch |
|---|---|---|---|---|---|
| 298 | charis_autoreduce.pro:101 | WRONG | `readcol,'obs*txt'` relies on a wildcard inside `openr` | `obslog.txt` | Y (pending) |
| 299 | charis_autoreduce.pro:107-108 | DEAD | `step4:` label inside an `if` block; labels never jumped to except `step11` | tidy | Y (pending) |
| 300 | charis_autoreduce.pro step 10 | WRONG | >1 detection → `fwdmod` writes several `synth_throughput*.txt` → `extract_1d_spectrum` opens a GUI (#103) → "auto" reduction stops | pass file explicitly | Y (pending) |
| 301 | charis_autoreduce.pro:34,73-74 | DEAD | unclosed quote in a print; "search SIMBAD if no target" stub; `quicklook` unused | tidy | Y (pending) |
| 302 | repo-wide (~30 files) | GEOM | 201×201 frames, centre 100, pixscale 0.0162, 22-channel ranges | Python: derive everything; IDL: only crashers patched | N |
| 303 | repo-wide (31 files, ~91 calls) | DEAD | `dialog_pickfile` GUI selection; `read` prompts; `set_plot,'x'`/`window`; `stop` | Python: explicit args, `savefig`, exceptions | N |
| 304 | repo-wide | WRONG | `median()` without `/even` in 30 call sites; `aper` argument-order slips | Python: `idl_median(even=)` where output-affecting | N |
| 305 | `.DS_Store` ×25, `*.bk`, `arch/`, `backup/`, `charis_rdiklip_old.pro`, `charis_autoreducebk.pro` | DEAD | cruft committed | delete; add to `.gitignore` | Y (pending) |

### 14. Found by the constants inventory (fork tree)

| # | File:line | Sev | What is wrong | Fix | IDL patch |
|---|---|---|---|---|---|
| 306 | imprep:141, register_cube:217, psfsub ×6, get_charis_mueller_matrix:24 | WRONG | three latitudes (19.825504 / 19.825 / 19.823806) and two longitudes (−155.47602 / −155.4802) in use | one `[site]` value; IDL reads `charis_get_constant('lat'/'lng')` | Y |
| 307 | pdi/mueller_matrix/component_models/m3_model.pro:17-21 | WRONG | post-2026-recoat M3 fit used for all dates, no epoch switch | coefficients from the epoch table (`m3_*`, `mjd=`) | Y |
| 308 | pdi/mueller_matrix/charis_pdi_mueller_matrix.pro:34 | WRONG | parity-diagonal omitted "because flip already applied", but the flip only happens for MJD>60933 → pre-NBS PDI gets wrong matrix | tie to `parity` | Y (pending) |
| 309 | pdi/mueller_matrix/component_models/hwp_model.pro:8-13 | GEOM | precomputed 22-entry No/Ne tables used when `n_elements(lambda) eq 22`; other grids use Sellmeier → two code paths | always Sellmeier | Y (pending) |
| 310 | pdi/charis_pdi_register_cube.pro:130-131 | WRONG? | `PARITY` flips polcent **y**, whereas an x-flip should mirror x (and swap L/R) | confirm with author; config key | N (document) |
| 311 | pdi/charis_pdi.pro:185 | GEOM | `telluric_channels=[5,6,14,21]` includes the last channel; differs from atmosfit `[5,6,14,15]` and empcomp `[5,6,7,14,15,16]` | wavelength ranges in config | Y (pending) |
| 312 | tools/contrast_calc/charis_calc_raw_contrast.pro:249 | WRONG | `[2.6]*aperradf` — typo for `[2,6]` | fix | Y |
| 313 | charis_specphot_cal.pro:887-888 | WRONG | header `Sky_In=2`/`Sky_Out=6` hard-coded, not tied to the actual `sat_skyrad` | write actual values | Y (pending) |
| 314 | charis_register_cube.pro:145 vs help | DEAD | `fittrim` default 4, help says 2; :1135 trims `2:n-3` ignoring `fittrim` | honour `fittrim` | Y (pending) |
| 315 | charis_register_cube.pro:1118 | GEOM | `PSFcensrefcorr[*,10]` channel 10 hard-coded | refslice | Y (pending) |
| 316 | prelims/charis_nbs_corr.pro:12-14 | WRONG | `rotzp` alternatives 34.7/32.3/33.5 left as comments; 113° magic | config `rotzp_deg`, documented | Y |
| 317 | charis_specphot_cal.pro:255, _unsat:205, raw_contrast:149 | WRONG | missing-MJD fallback `mjd_in=58000` after a y/n prompt | require MJD or explicit `--mjd` | Y (pending) |
| 318 | pdi/mueller_matrix/get_charis_mueller_matrix.pro:12 | WRONG | re-opens `./data/raw/CRSA<origname>_cube.fits` ext 3 for `D_IMRANG`/`PARANG` | imprep copies keys into prep header | Y (pending) |
| 319 | setup/charis_path.pro:14,22 (fork) | WRONG | `/home/thomasmc/charis-dpp/idl/`, `/Volumes/HAL/models/` | env var | Y |
| 320 | psfsub/*: `na` read from .info but never written; `spang*`,`spmask`,`cenmask` likewise | DEAD | schema gaps | config defaults | Y (pending) |
| 321 | psfsub/charis_adialoci.pro:322, charis_adiklip.pro:257,549,554, charis_rdiklip_old.pro:272,616-627, subroutines/{snratio:157, charis_snratio_sub:79,184, profrad_tc:30,45, setcolors:298, craigm/cmrebin:30-31, craigm/cmsv_rvtype:475,485} | CRASH | the fork's auto-formatter rewrote `mean(`, `max(`, `finite(`, `size(`, `total(` as `[...]` indexing; under `compile_opt idl2` these index an undefined variable | restored to function calls | Y |
| 322 | setup/charis_get_constant.pro (all ~100 callers) | WRONG | without `mjd=` the epoch table returns the **latest** epoch, so 2023 data get Chen+2023 pixscale/angoffset, post-NBS parity and post-recoat M3 unless every caller passes `mjd=`; only imprep and get_charis_sat_contrast did | session default via `charis_set_epoch` (common block), called by imprep (header MJD) and every top-level PDI/registration/specphot routine (`.info` OBSDATE); `charis_northup` uses the frame MJD | Y |
| 323 | setup/charis_get_constant.pro, charis_dpp/data/instrument.toml | WRONG? | satellite-spot contrast epochs stop at 2017-08-30; README (T. Currie priv. comm.) lists further changes at 2021-12-31, 2022-05-20 and separate 500 Hz-modulation values from 2022-10-20 | add epochs once the values are confirmed in magnitudes; keep `satcon*` keys | N (document) |
| 324 | subroutines/craigm/ (whole directory, 114 files) | CRASH | the fork's auto-formatter added `compile_opt idl2` to Craig Markwardt's library; files that index common-block variables with parentheses after an `@fxfilter`/`@dxcommon` include (fxgfiltered, fxgread, fxgseek, fxgwrite, fxpclose, dx*) then compile `filterflag(unit)` as a function reference, after which every later file assigning `filterflag(unit) = …` (fxgopen, fxpopenr) fails with a syntax error; the failure depends on compile order, so it only shows when the tree is compiled as a whole | the directory is third-party and carried no intentional patch, so it is restored verbatim to upstream `a07976a` (this also undoes the #321 fixes in cmrebin/cmsv_rvtype the same way); `charis_compile_check` now skips `@include` fragments and requires every declared routine of a file to compile | Y |
| 325 | reg/charis_getcenter.pro + subroutines/charis_getcenter.pro; specphotcal/gbfunc.pro + charis_calc_spec_covar.pro:413 | MINOR | two routines are defined twice in the tree (both copies already in upstream `a07976a`); which copy runs depends on `!path` order and compile order. The two `charis_getcenter` copies share the algorithm but `reg/` writes a stray `im.fits` and `subroutines/` carries dead `rotat` code; the two `gbfunc` copies differ in parameterisation (`P[4]` vs a `rho` variable) and only the in-file copy is ever current when `charis_calc_spec_covar` runs | reported as `DUP` by `charis_compile_check`; delete one copy of each once the registration validation says which `charis_getcenter` ran | N (document) |
