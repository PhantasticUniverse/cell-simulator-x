# Phase 11.2.D Notes — Hemoglobin Adair Binding + pH Buffer on GPU

Date: 2026-05-01

## Headline

The post-RK4 hemoglobin (Adair 4-site cooperative binding) update and the
pH buffer (algebraic from lactate) are folded into the inner-step loop of
the same compute kernel as 11.2.C. After 1 s of simulated time (1000 RK4
ms-steps × 38 species + 1000 Hb post-Euler steps), every species agrees
with the f64 CPU reference to **0.0023% relative error worst-fit** (Na⁺ —
unchanged from 11.2.C) and Hb saturation matches to **2e-6 absolute**
(0.937480 CPU vs 0.937478 GPU starting from 0.75).

The full pipeline through 11.2.D matches the CPU `FullyIntegratedSolver`
end-to-end with the only remaining caveat that the three inline
homeostasis hacks (basal NADPH, basal GSH, ATP regen) are still skipped on
both CPU reference and GPU kernel — those are 11.2.E.

## What landed

### Files modified

- `shaders/compute/biochem_full.wgsl` — extended in place (kept as the
  single shared kernel for both 11.2.C and 11.2.D paths):
  - **Uniforms**: 22 new fields covering pO₂/pCO₂/T, total Hb, base P50,
    Bohr/DPG/T/CO2 allosteric coefficients, kinetic rate constants,
    pH-buffer parameters, base Adair K₁..K₄, plus an `enable_hb` toggle.
  - **Bind group**: a third storage buffer at binding 2 holds per-cell
    Hb saturation. Always present so the bind-group layout is stable
    across 11.2.C and 11.2.D call paths.
  - **New helpers**: `compute_ph(lactate)` (Van Slyke buffer), then
    `effective_p50(ph, dpg)` (Bohr + DPG + van't Hoff + CO₂),
    `saturation_adair(po2, k1..k4)` (4-site cooperativity), and
    `hb_post_step(lactate, dpg, sat)` (Euler update of saturation
    toward equilibrium).
  - **Inner loop**: each RK4 ms-step is followed by `hb_post_step` when
    `enable_hb = 1`. The CPU side does the same: one `RK4Integrator::step`
    + one `HemoglobinSolver::step` per millisecond.

- `src/compute/biochem.rs`:
  - `FullBiochemUniforms` and `FullBiochemBatchConfig` extended with the
    22 new fields. Defaults match `HemoglobinSolver::default()`,
    `AllostericParameters::default()`, and `PhBufferModel::default()`.
  - `run_full_biochem_batch` is now a thin wrapper over
    `run_full_biochem_internal(... None ...)` (Hb path off by default —
    the 11.2.C parity test keeps working unchanged).
  - **New** `run_full_biochem_batch_with_hb(ctx, pools, hb_saturations,
    config)` enables the Hb path and updates `hb_saturations` in place.

### Files added

- `tests/biochem_gpu_hb_parity.rs` — runs the same 1-second sim on CPU
  and GPU with `INITIAL_HB_SAT = 0.75` and `PO2_MMHG = 100`, asserts
  per-species ≤ 1% rel-err **and** `|Δsaturation| < 0.01`.

### Files unchanged but verified

- `tests/biochem_gpu_full_parity.rs` (Phase 11.2.C test) — still passes
  identically (Na⁺ at 0.0023% worst-fit). The new `enable_hb` defaults to
  `false`, the kernel skips the Hb code path, and the dummy `hb_state`
  buffer is allocated but never read back.

## Numerical details

- **Adair constants**: precomputed once on the CPU side via
  `AdairConstants::default()` (which calibrates K₁..K₄ to give P50 = 26.8
  mmHg and Hill ≈ 2.7) and passed in as four uniform f32s. The kernel
  never runs the bisection — it just multiplies the four constants by a
  per-step `p50_ratio = base_p50 / effective_p50` to get the conditional
  binding curve.
- **Henry's law approximation**: O₂ → mM via `pO₂/760`, matching the CPU
  side verbatim (Roughton 1972 simplification, not actual solubility).
- **`STANDARD_PH = 7.4`** is hardcoded in WGSL as a `const f32` to mirror
  the CPU code's compile-time constant.
- The post-step uses Euler — not RK4 — and runs at the same dt as the
  inner ODE step. The CPU reference does the same; the test thus checks
  algorithmic parity rather than higher-order accuracy.

## Architectural decision: shared kernel, conditional code path

Originally I considered making `biochem_full.wgsl` a 38-species-only kernel
and adding a second `biochem_full_hb.wgsl` for the Hb-extended version.
That would have meant ~500 lines of WGSL duplication (helpers + rate
functions + RK4) since WGSL has no preprocessor.

Instead, the kernel was extended in place with an `enable_hb` toggle and
the `hb_state` buffer is always bound. The 11.2.C path passes a dummy
`hb_state` and disables the post-step; the 11.2.D path passes real
saturations and enables it. One kernel, one host wrapper (`internal`), two
public entry points. The 11.2.C parity test required only a one-line
change (replace explicit struct literal with `..Default::default()`) and
still passes bit-identically against its prior CPU reference.

## Tests added

- `tests/biochem_gpu_hb_parity.rs::full_biochem_with_hb_cpu_gpu_parity_one_second`:
  1-cell, 1-second run with PO₂ = 100 mmHg and initial Hb saturation 0.75.
  CPU and GPU both relax to ≈ 0.937 within 2e-6 absolute. Skipped if no
  GPU adapter.

## Test counts after 11.2.D

- 186 lib tests (unchanged).
- 133 integration tests (was 132; +1 from `biochem_gpu_hb_parity`).
- Validation suite: 7/8 (unchanged).

## What's next (Phase 11.2.E)

Add the three inline homeostasis hacks from
[full_integration.rs:251-294](../src/biochemistry/full_integration.rs):

1. Basal NADPH consumption (coefficient 0.0003 with NADPH/NADP⁺ ratio
   modulation).
2. Basal GSH oxidation (coefficient 0.001 with GSH/GSSG ratio modulation).
3. ATP-deficit regen (coefficient 0.2, deficit toward 1.8 mM ATP).

These are pure-functional expressions of state — no control flow needed.
After 11.2.E the GPU path produces bit-equivalent (within f32 noise)
trajectories to the full `FullyIntegratedSolver::step` exactly. That
milestone closes Phase 11.2 (biochemistry) and clears the runway for
Phase 11.3 (Skalak + WLC + DPD).
