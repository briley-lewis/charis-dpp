# Handoff: state of the Python port (branch `claude/gracious-turing-m3b1zl`)

Written for a local session that has IDL and the real data. Read this first,
then `docs/WORKFLOW.md` and `BUGS.md`.

## What exists

* `charis_dpp/` — a complete Python port of the IDL pipeline (`idl/`), on the
  branch `claude/gracious-turing-m3b1zl` of `briley-lewis/charis-dpp`
  (9 commits, tip `f7b5ec4`). It does **not** need IDL. Install with
  `pip install -e ".[dev]"`; CLI is `charis <step>` (51 subcommands,
  `charis autoreduce` for the chain). Tutorial: `docs/tutorial.ipynb`.
* Source of the port: the fork `thomasmcintoshh/charis-dpp-mm-updates`,
  which was copied into `idl/` and then patched.
* Configuration: dated instrument epochs in `charis_dpp/data/instrument.toml`
  (sat-spot contrast change MJD 57995, Chen+2023 constants from 2024-10-18,
  NBS parity flip MJD 60933, IR M3 recoat 2026-04-16), one `.toml` per dataset
  (`charis newobs` writes it; `charis migrate-info` converts an old `.info`;
  the legacy `.info` is still written for IDL), and an optional
  `[instrument]` override table per dataset.
* New features asked for: `charis register --method xcorr` (no satellite spots;
  `--xcorr-center symmetry|centroid|given`), hardened HWP matching
  (`charis pdi-hwp-match`: log / `hwp_angles*.csv` / `RET-ANG1` headers,
  midnight roll-over, warnings for moves during an exposure, per-frame report,
  `--time-ref start|mid|end`, default `mid`).
* Tests: `pytest -m "not slow"` (unit, seconds) and `pytest -m slow`
  (synthetic end-to-end chains: ADI, KLIP, RDI, SDI, disks, PDI, autoreduce).
  98 pass. `tests/unit/test_idl_parity.py` checks the IDL constant table
  against `instrument.toml`.
* `BUGS.md`: 320 IDL defects found by reading the code. All are handled in
  Python. 65 are patched in `idl/` (grep `BUGS.md #n`); 220 are marked
  `Y (pending)`; the rest documented only.
* IDL changes: `setup/charis_get_constant.pro` now holds the epoch table
  (`mjd=` keyword; without it the latest values are returned, as before);
  `setup/charis_path.pro` reads `$CHARIS_DPP_PATH` / `$CHARIS_DPP_BTSETTL`;
  `cals/get_charis_sat_contrast.pro` uses the table; new
  `subroutines/{get_xycind,calc_centroid_spots4,medfitsme}.pro`.

## What was NOT possible here and needs the local machine

1. **No IDL was available.** None of the IDL patches have been compiled or
   run. First thing to do locally: `.compile` the patched files (list them with
   `grep -rl "BUGS.md" idl --include=*.pro`) and run one dataset through the
   IDL chain; fix any syntax slip.
2. **No real data.** The Python pipeline was validated only against synthetic
   data with known truth (`charis make-synthetic`). The regression harness
   against IDL products is scaffolded but has never run:
   put a dataset under `$CHARIS_DPP_TESTDATA/<name>/` with `raw/CRSA*_cube.fits`,
   the `.info`/`.toml`, `idl/reduc.log`, `idl/reduc/reg/*reg_cal.fits`,
   `idl/reduc/proc/*`, then `pytest tests/regression`
   (tolerances in `tests/regression/tolerances.toml`; `charis compare-fits A B`
   compares two files by hand). Expect differences from the deliberate fixes
   (north-up WCS, even-size centre convention, `/even` medians) — see
   `docs/WORKFLOW.md`.
3. **Observing logs** (two Google Sheets) could not be read from the cloud
   session. Use them locally to build each night's `.toml` (frame lists, HWP
   log, modulation, astrogrid spacing, flats).
4. `charis newobs` needs network for SIMBAD/2MASS; offline, fill `[target]`
   by hand (magnitudes, spectral type, parallax) before `specphot-cal`.

## Open questions worth settling with real data

* Mueller matrix parity diagonal (BUGS.md #308): config flag
  `pdi.mueller.flip_matrix`, default off (the fork's behaviour).
* PDI `polcent_offset` parity handling (BUGS.md #310): the fork flipped the
  y offset post-NBS; the port mirrors x and swaps beams. Check on post-2025 PDI data.
* `hwp_time_ref` default `mid` (the fork's behaviour); the IDL comments
  disagreed about start vs end.
* Registration `method=xcorr` reference centring: `symmetry` is the default
  but on synthetic data `centroid` was more reliable; check on real spot-less data.
* The 220 pending IDL patches: apply as needed; each row in `BUGS.md` says what to change.
