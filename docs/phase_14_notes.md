# Phase 14 Notes — Storage Lesion at Physiological Timescale

Date: 2026-05-02

## Headline (Phase 14.A)

The simulator can now run a 42-day storage curve, sampling metabolite
state at each storage day. Three-tier multi-rate scheme:

- **Inner (1 ms)**: `FullyIntegratedSolver` 38-species ODE.
- **Middle (~seconds)**: biochemistry equilibration window per outer
  step.
- **Outer (1 day)**: `StorageLesionModel` advances; pump efficiency,
  leak conductance, and oxidative stress recompute per Hess 2010 /
  Luten 2008 / Dumaswala 1999 envelopes.

This is the foundation the headline-paper deliverable rests on. Phase
14.B closes the gap between simulator output and Hess 2010's quantitative
ion targets.

## What landed (Phase 14.A)

### `src/storage/` — new module

- `StorageSimConfig { end_day, days_per_step, seconds_of_bio_per_step,
  bio_dt_sec, force_atp_dpg_targets }` — knobs.
- `StorageSample` — per-day metabolite snapshot (ATP, 2,3-DPG, GSH,
  GSSG, H₂O₂, NADPH, Na⁺, K⁺, lactate, pump efficiency, leak multiplier,
  oxidative stress).
- `StorageCurveSimulator` — orchestrator. Holds `FullyIntegratedSolver`,
  `MetabolitePool`, `StorageLesionModel`. `step_to_day()` updates the
  envelope, applies it to the solver (oxidative stress, pump Vmax, leak
  conductances), optionally forces ATP / 2,3-DPG to the empirical
  targets, then runs biochemistry for `seconds_of_bio_per_step` to
  equilibrate the rest of the pool.
- `run()` — loops day-by-day from 0 to `end_day`.
- `sample_at_day()` — nearest-neighbour lookup on the recorded curve.
- `write_csv()` — exports the full trajectory (one row per day, 13
  columns).

### Validation targets (Hess 2010 / Luten 2008 / Zimrin 2009)

| Day | ATP (mM) | 2,3-DPG (mM) | Na+ (mM) | K+ (mM) |
|-----|----------|--------------|----------|---------|
| 0   | 2.0      | 5.0          | 10       | 140     |
| 14  | 1.5      | 0.5          | 25       | 120     |
| 42  | 0.5      | 0.0          | 60       | 90      |

### Tests added

`tests/storage_curve.rs`:

| Test | Status | Result |
|------|--------|--------|
| `day_0_targets` | ✅ | ATP 2.00, DPG 5.00, Na 10.0, K 140.0 |
| `atp_decays_per_hess_2010` | ✅ | d21 ATP = 1.024 (~1.0), d42 ATP = 0.502 (~0.5) |
| `dpg_depletes_per_zimrin_2009` | ✅ | d14 < 1 mM, d42 < 0.05 mM |
| `ion_gradients_trend_correctly` | ✅ | Na monotonic up, K monotonic down |
| `full_curve_has_one_sample_per_day` | ✅ | 43 samples (days 0–42) |
| `write_csv_succeeds` | ✅ | 1 header + 43 data rows |
| `ion_gradients_long_equilibration` | ⚠️ `#[ignore]` | d42 Na=25, K=130 with 30 s/day bio |

3 unit tests in `storage::simulator::tests` (lib).

ATP and 2,3-DPG match Hess 2010 quantitatively because the simulator
forces them to the envelope. The biochemistry-driven metabolites (ions,
GSH, NADPH) follow the right trends but don't yet quantitatively match
Hess 2010 — see "Phase 14.B" below.

## Test counts after 14.A

- 200 lib tests (was 197; +3 from `storage::simulator::tests`).
- 148 integration tests (was 142; +6 from `storage_curve`).

## What Phase 14.A demonstrates

1. **The simulator can run the full 42-day storage trajectory** in
   under a second of wall-clock with default config (1 s of bio per
   day = ~1 s wall-clock total).
2. **Per-day sampling works** and the CSV export gives the same data
   format as the published storage metabolomics tables (one row per
   day, columns per metabolite).
3. **The envelope-forcing scheme lets the rest of the metabolite pool
   find a self-consistent steady state** at each storage day. NADPH /
   GSH / lactate / glycolysis intermediates all sit in physiological
   ranges and respond to the modified pump / leak / oxidative-stress
   parameters.
4. **The directional response is correct**: Na⁺ rises monotonically,
   K⁺ falls monotonically, ATP and 2,3-DPG follow Hess 2010, GSH
   responds to the rising oxidative stress.

## What Phase 14.A does **not** yet demonstrate

The headline-paper deliverable would be a curve indistinguishable from
Hess 2010's reported values across all key metabolites. Phase 14.A
falls short in two specific places:

### Ion-gradient quantitative match

With `seconds_of_bio_per_step = 30` (already 0.04 × 30 × 42 ≈ 50 s of
wall-clock per curve), Na⁺ at day 42 reaches 25 mM. Hess 2010 reports
~60 mM. The gap is real: pump+leak kinetics have a roughly 1-hour
equilibration timescale in this simulator (set by enzyme Vmax and leak
conductance values that themselves came from in-vivo data), and 30 s
of simulated time is short compared to that.

Three concrete options for Phase 14.B:

- **Longer equilibration window** (`seconds_of_bio_per_step ≈ 300`):
  ~7 minutes of wall-clock, full ion equilibration. Acceptable for one-
  off curve generation; too slow for unit tests.
- **Analytic quasi-steady-state for ions**: solve `3·pump_rate(Na) =
  leak_rate(Na)` at each storage day via bisection, set Na/K directly,
  let biochem only equilibrate the fast subsystems (NADPH / GSH).
  Mathematically clean; requires changes inside the simulator.
- **Retune pump / leak rate constants** so the equilibration timescale
  matches the chosen integration window. Risks breaking the day-0
  homeostasis numbers that Phase 11.2.E validated against.

The right answer is probably (2) — analytic QSS — because it preserves
the day-0 calibration and removes the wall-clock dependency entirely.

### GSH / GSSG response

The simulator preserves GSH at ~2.5 mM through day 42 because the
glutathione cycle's Phase 10-tuned Vmax values keep up with the rising
oxidative stress. Storage literature reports GSH falls to ~1 mM by day
42 (D'Alessandro et al. 2020). Closing this gap is part of Phase 14.B.

## Why this is the headline foundation

Once Phase 14.B closes the quantitative gap, the simulator produces a
**publishable storage curve**: one CSV per simulation run, one figure
per simulation run, in seconds of wall-clock. With the open-source
release of Phase 17, any reviewer can run `cargo test storage_curve` to
verify the result. That's the durable scientific moat the strategic
plan identified — long-timescale disease modeling at multi-week
physiological scale that no competitor simulator (HemoCell, uDeviceX,
Pivkin–Karniadakis) currently provides.

## What's next

- ✅ **Phase 14.B (this commit)**: analytic quasi-steady-state for the
  ion subsystem. `solve_ion_qss(pump_eff, leak_mult, atp_mM)` does
  bisection on the Na+ root of `3·pump_rate(Na) = leak_Na(Na)` and
  derives K+ from the same pump rate. Gated by
  `force_ion_qss = true` in `StorageSimConfig` (default `true`). With
  QSS enabled, day-0 lands at the pump+leak equilibrium (Na ≈ 10 mM,
  K ≈ 145 mM) and day-42 lands at Na ≈ 96 mM / K ≈ 52 mM in seconds of
  wall-clock — a measurable storage trajectory without minutes of bio
  equilibration per day.
- **Phase 14.B' — envelope re-fitting**: Phase 14.B reveals a
  parameter-identifiability finding. With the simulator's enzyme
  kinetics fixed (validated in Phase 11.2.E), the linear envelope
  parameters (`pump_efficiency_decay_per_day = 0.02`,
  `leak_increase_per_day = 0.015`) drive the QSS to Na ≈ 96 mM at day
  42 — pathologically beyond Hess 2010's reported ~60 mM. Algebraically,
  no positive-decay-rate / positive-leak-rate pair can simultaneously
  match Hess 2010 day-14 (Na ≈ 25, K ≈ 120) AND day-42 (Na ≈ 60, K ≈ 90)
  with linear envelopes. Either nonlinear envelopes or a different
  enzyme calibration would close the gap. Documenting this as a
  parameter-identifiability finding (mirroring the Phase 10 NADPH/PPP
  refit) is the right move; full fitting is Phase 14.B'.
- **Phase 14.C**: deformability-decline coupling. The plan called for
  reproducing the 42-day deformability decline curve. With ATP, 2,3-DPG,
  GSH, and oxidative stress all tracked over time, the next step is to
  feed those into a deformability index (e.g., via the existing
  `SpectrinModulator` ATP→stiffness coupling) and produce a single
  scalar deformability vs. day curve.
- **Phase 14.D**: sweep the storage parameter space — supercooled
  vs. standard storage, additive solutions (AS-3, SAGM, PAGGSM),
  comparator runs.

## Phase 14.B addendum — what landed

### `solve_ion_qss(pump_efficiency, leak_multiplier, atp_mM) -> IonQss`

Bisection-based root finder for the Na+ ion quasi-steady-state under
modified pump+leak parameters. K+ derives from the same pump rate. Three
unit tests:

- `ion_qss_day_0_matches_physiological` — at envelope = (1, 1, 2 mM
  ATP), QSS lands at Na ≈ 10 mM, K ≈ 140 mM (within 5 mM).
- `ion_qss_day_42_in_pathological_range` — at envelope = (0.16, 1.63,
  0.5 mM ATP), QSS lands at Na > 40 mM, K < 130 mM. Quantitative match
  to Hess 2010 (Na ≈ 60, K ≈ 90) is gated on Phase 14.B' envelope re-fit.
- `ion_qss_monotone_in_pump_efficiency` — Na strictly increases as pump
  efficiency drops (sanity check on the bisection).

### Storage simulator integration

`StorageCurveSimulator::step_to_day()` now (when `force_ion_qss = true`)
calls `solve_ion_qss` and writes the result to the metabolite pool
before running the bio equilibration window. The bio equilibration is
then a small correction (other species: NADPH / GSH / lactate)
rather than a multi-hour drive on the ion system.

`run()` no longer takes a special "initial physiological" day-0 sample
— day 0 goes through the same `step_to_day(0.0)` path so all samples
come from the same procedure. Day-0 K is now ~145 mM (the QSS value)
instead of the bare-pool 140 mM.

### Default behaviour change

`StorageSimConfig::default()` enables `force_ion_qss = true`. Existing
callers that want the legacy slow-equilibration behaviour can set
`force_ion_qss: false` explicitly (the
`ion_gradients_long_equilibration` `#[ignore]` test exercises this
path).

### Test counts after 14.B

- 203 lib tests (was 200; +3 from `ion_qss_*`).
- 148 integration tests (unchanged from 14.A; existing tests still
  pass with the QSS-on-by-default config).
