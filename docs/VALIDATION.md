# Validation: old IDL → new IDL → Python on real 2023-07-30 CHARIS PDI data

Status: **set-up written, no runs yet** (2026-11-05). Results tables are filled in
as each comparison completes.

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
   deliberate differences (`/even` medians, even-size centre convention,
   north-up WCS header, derotation kernel).

Tool: `charis compare-trees A B --mode identical|port [--allow-list f] [--out d]`
(`charis_dpp/testing/compare_trees.py`), which pairs `prep/`, `reg/`, `proc/`,
`reduc*.log` and `hwp_info.txt` by relative path.

## Layout on mueller

```
~/charis-val/
  old-idl/      git worktree at tag old-idl-2025-04   (IDL 8.9, ~/idl/idl89)
  new-idl/      git worktree at branch integrate-port
  datasets/
    frames-2023-07-30.md          frame lists worked up from the observing log + headers (copy in docs/validation/)
    HD236894_bb/                  polarized standard, 16 s broadband PDI
      data/raw/CRSA*_cube.fits    symlinks to /home/shared/exoserver/Subaru-Europa/2023-07-30/HD236894-broadband-pol/
      pdi_flats/                  extracted PDI flat cubes
      HD236894_bb.info            hand-written; the same file for every run
      HD236894_bb.toml            charis migrate-info HD236894_bb.info
      old/ new/ py/               one working directory per run (each gets its own reduc/)
    Europa_bb/                    same shape; raw = extracted cubes 123628–123847
  compare/                        compare-trees outputs (compare.csv / compare.md per pair)
```

Set-up commands (run once connectivity to mueller is available):

```bash
ssh mueller
cd ~/charis-dpp && git fetch origin integrate-port old-idl-2025-04 old-idl-2024-06
mkdir -p ~/charis-val/datasets ~/charis-val/compare
git worktree add ~/charis-val/old-idl old-idl-2025-04
git worktree add ~/charis-val/new-idl integrate-port
# old branch has the hard-coded path in idl/setup/charis_path.pro: point it at the worktree (do not commit)
sed -i "s|charispath=.*|charispath='$HOME/charis-val/old-idl/'|" ~/charis-val/old-idl/idl/setup/charis_path.pro
# Python side
conda activate charisenv && pip install -e ~/charis-val/new-idl
# compile gate for the new tree (IDL 8.9)
IDL_PATH="+$HOME/charis-val/new-idl/idl:+$HOME/IDLAstro:<IDL_DEFAULT>" ~/idl/idl89/bin/idl -e charis_compile_check
```

Per run, IDL is started with `IDL_PATH` pointing at exactly one of the two
worktrees (plus IDLAstro), and `CHARIS_DPP_PATH` set to the new worktree for
the new-IDL run.

## Step list (identical for every run)

```
charis_newobs, date='20230730', object='HD236894', filter='low'    ; Europa: arbitrary name, sptype -999
charis_imprep                                                      ; then set fnum_sat / fnum_sky by hand
charis_pdi_flatfield, 'x.info', brightflatdir='../pdi_flats/'
charis_pdi_subtract_sky, 'x.info', /scalesky
charis_pdi_register_cube, 'x.info', astrogrid=11.2                 ; only if the astrogrid was on
charis_pdi_specphot_cal, 'x.info', modamp=25, ap_factor=3, starlib=1   ; standard only
charis_pdi_hwp_match, 'x.info', /use_headers
charis_pdi, 'x.info', /meanadd, /save_derot_seq                    ; and a second run with /no_mm_ip_corr, /no_mm_crosstalk_corr
```

Python: the `charis imprep --pdi … charis pdi` sequence from `WORKFLOW.md` on the
migrated `.toml`.

Controlled variables (identical in every run, recorded here when fixed):

| variable | value | why |
|---|---|---|
| `fnum_sat`, `fnum_sky` | HD 236894: sat 123560–123583 (6 clean cycles), no sky; Europa: sat 123628–123847 (55 cycles), sky 123848–123855 | from `docs/validation/frames-2023-07-30.md`; the DPP mis-tags PDI frames as sky, so set by hand |
| flats, darks | bright PDI lamp flat 123412–123416 (5 × 5.9 s, no lamp-off frames); darks 123417–123420 (200.6 s); twilight sky flats 123926–933 as the alternative | the only PDI flats taken that night |
| HWP angle source | raw-frame `RET-ANG1` (cubes carry none); the angle at the start of the ramp | the waveplate moves during the last read(s) of each ramp, so `hwp_time_ref=start` semantics; a per-read check is an open question |
| registration method, `polcent_offset` | non-spot: `X_GRDST=OFF` on every frame of the night (`docs/validation/frames-2023-07-30.md`); IDL `charis_pdi_register_cube` with the PSF-centroid method for the star, `method=xcorr` (Python) / cross-correlation for Europa; `polcent_offset` left at the instrument default | no satellite spots anywhere, so `astrogrid=` is not used |
| `hwp_time_ref` | `mid` (Python default; IDL per code) | HWP log vs exposure timing |
| `pdi.mueller.flip_matrix` | `false` | fork behaviour on both sides; BUGS.md #308 revisited afterwards |
| `angoffset` | 2.2 (epoch `satspot-recal-2017`) | 2023 epoch; changing it means rerunning `charis_pdi` |
| Mueller | two runs: corrected and `/no_mm_*` | isolates the Mueller term |
| `unres_pol_sub`, `bg_sub` | off | |
| out-of-FOV pixels | NaN before `charis_pdi` | coverage NaNs must agree (`finite_agree`) |

## Deliberate differences expected in `identical` mode (allow-list candidates)

Filled in from the first old-vs-new run. Every entry must cite a BUGS.md row or
a commit on `integrate-port`.

| product | reason |
|---|---|
| | |

## Results

### HD 236894 (polarized standard)

| pair | mode | files | pass | fail | notes |
|---|---|---|---|---|---|
| old IDL vs new IDL | identical | | | | |
| new IDL vs Python | port | | | | |

Recovered P and AoLP vs catalogue, per run: TBD.

### Europa

| pair | mode | files | pass | fail | notes |
|---|---|---|---|---|---|
| old IDL vs new IDL | identical | | | | |
| new IDL vs Python | port | | | | |
| new IDL vs 2023 products | (gross structure only) | | | | |

Checks specific to Europa: the c0019–c0036 cycle-pairing anomaly around frame
123705 in the 2023 products; the sign of disk-median Q; integrated Q and U per
cycle against time for the corrected and uncorrected runs (a trend in the
corrected run means the Mueller model or the non-sidereal bookkeeping is off);
bootstrap uncertainties on q, u, P and AoLP over the 55 cycles.

## Verdicts

One line per difference found, with its attribution (BUGS.md row, commit, or
"unexplained → new BUGS.md row").
