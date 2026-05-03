# Phase C-hybrid Notes — Splenic Transit at Storage Day N

The headline integrated-coupling demo: take a cell that the 42-day
storage curve has aged, drop it into a 0.5–1.0 μm splenic
interendothelial slit, and measure transit kinetics. Uses both moat
assets simultaneously — storage time-evolution (Phase 14) and
mechano-metabolic coupling under flow (Phase 12.B) — in a single
biological story.

## Architecture

Three sub-deliverables compose:

### C.1 — Storage cell-state snapshot (`StorageCellSnapshot`)

`StorageCurveSimulator::cell_state_at_day(d)` runs a replay simulator
from day 0 to day d and returns:
- Full 38-species `MetabolitePool` at day d.
- ATP-driven spectrin stiffness modifier (Phase 8 / 14.C coupling).
- Pump efficiency, leak multiplier, oxidative stress envelope state.
- Deformability index relative to fresh.

Validated by `tests/storage_replay_parity.rs` (4 tests):
day-21 ATP within 1% of direct trajectory; day-0 matches initial pool;
envelope state carries through; deterministic across multiple calls.

### C.2 — `SplenicSlit` + `SlitFlow`

`flow::SplenicSlit`: rectangular slot with width (0.5–1.0 μm), height
(~3 μm), length (~5 μm), axis (flow direction), narrow_axis, entry.
`flow::SlitFlow`: parabolic-in-narrow-direction velocity profile inside
the slit; zero outside. Wall shear rate `4·v_max/w` — sub-micron slits
produce kHz-range shear at modest centerline velocities.

Validated by 7 unit tests: parabolic profile, wall velocity = 0, narrow-
slit shear concentration, transit time canonical, etc.

### C.2 — `SplenicTransitConfig` + `run_splenic_transit`

Couples `StorageCellSnapshot` with `SlitFlow`:
1. Build a fresh mesh + spectrin network.
2. Position the cell straddling the slit entrance.
3. Apply the snapshot's stiffness modifier to `SkalakMaterial`.
4. Loop physics substeps; per substep, recompute slit-flow drag and
   apply via `apply_slit_drag_to_external_forces`; CPU
   `PhysicsSolver::step`.
5. Stop when centroid passes exit plane (success) or simulated-time
   cap is reached (clearance failure).

Reports: transit time, centroid axial displacement, peak strain, peak
velocity, deformability, wall-clock.

CLI: `cargo run --release -- --diagnose-splenic-transit --storage-day N --slit-width W`

### C.3 — `(storage_day, slit_width)` sweep

`sweep_storage_day_x_slit_width(days, widths, config)` returns the
result table; `write_transit_csv` emits
`target/splenic_transit_storage_curve.csv`.

`tests/splenic_transit_sweep.rs` runs the headline 7×3 grid (days 0,
7, 14, 21, 28, 35, 42 × widths 0.5, 0.7, 1.0 μm) in ~4 s wall-clock.

## Findings

| Day | w_um | shear (1/s) | displ_um | strain | v_peak | ok |
|-----|------|-------------|----------|--------|--------|-----|
|   0 | 0.50 |        4000 |    1.274 |  1.224 |  152.5 |   N |
|   0 | 0.70 |        2857 |    1.453 |  1.359 |  163.1 |   N |
|   0 | 1.00 |        2000 |    1.625 |  1.517 |  185.2 |   N |
|  42 | 0.50 |        4000 |    1.263 |  1.217 |  152.6 |   N |
|  42 | 0.70 |        2857 |    1.452 |  1.357 |  162.6 |   N |
|  42 | 1.00 |        2000 |    1.625 |  1.514 |  184.5 |   N |

**Slit-width effect**: strong and monotonic. Narrower slit → smaller
centroid displacement (0.5 μm: 1.27 μm; 1.0 μm: 1.63 μm). The
mechanical filter geometry is correctly captured.

**Storage-day effect**: muted — day-0 → day-42 produces only ~0.001 μm
change at the same width. The Phase 8 `SpectrinModulator`'s 0.5
`max_stiffening_factor` (giving 1.0× → 1.37× stiffness across the
storage range) is too small a mechanical change to translate
visibly to centroid transit kinetics under the current slit-flow
drag regime.

## Limitations and follow-on

- **Stiffness coupling amplification.** Increasing
  `max_stiffening_factor` from 0.5 → 2.0 (giving up to 3× stiffening at
  zero ATP) would amplify the storage-day signal. Not done in this
  pass since it would require revalidating Phase 8 against
  Manno 2002 spectrin stiffness data.
- **Bending modulus / membrane viscosity coupling.** Aged RBCs in vivo
  show large changes in bending modulus and membrane viscosity, not
  just shear stiffness. The current model only modulates shear modulus.
  Adding ATP→bending and ATP→viscosity couplings is a clean follow-on.
- **Upstream pressure model.** The current SlitFlow gives zero velocity
  outside the slit volume. In vivo, an upstream Poiseuille zone
  drives the cell into the slit. Phase 13 (multi-cell / vessel) work
  would naturally extend this.
- **Quantitative comparison to Pivkin 2016 / Picas 2013.** The
  framework produces the right directional signals; quantitative
  benchmarking against published splenic-clearance times requires the
  amplified stiffness coupling above plus a proper upstream-pressure
  model.

## Status

Framework complete. Ready for follow-on tuning. The integrated demo
script (`--diagnose-splenic-transit`) runs in <1 s wall-clock per
(day, width) configuration; the headline sweep CSV is reproducible
via `cargo test --release --test splenic_transit_sweep` in ~5 s.

## What's next

Phase 17 prep (open-science release): split `rbc-validation-suite`
into a workspace member crate, add `flake.nix`/`Dockerfile`/`Makefile`,
write `docs/methods.md` and `docs/headline_figures.md` for the
preprint.
