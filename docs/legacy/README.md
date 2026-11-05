# Legacy material from the `2025-updates` branch

- `local-registration-edits-2025.patch` — uncommitted edits (spring 2025) to
  `idl/reg/charis_register_cube.pro` and `idl/pdi/charis_pdi_register_cube.pro`
  that added `fitgauss`/`fittrim`/`rin_corr`/`rout_corr`, session logging and a
  combined `sqrt(R^2+L^2)` frame. The `integrate-port` branch carries the same
  functionality from the McIntosh fork (combined frame kept as L + R, the
  single-sum total intensity). Kept for the record only.

Tags: `old-idl-2025-04` (tip of `2025-updates`, last code run before the port)
and `old-idl-2024-06` (closest to the state used for the June 2024 Europa
reduction: upstream + Kellen Lawson's Mueller-matrix PDI tree).
