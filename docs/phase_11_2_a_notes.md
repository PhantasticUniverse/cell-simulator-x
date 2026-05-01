# Phase 11.2.A Notes — Glycolysis on GPU

Date: 2026-05-01

## Headline

The 11-enzyme glycolytic backbone runs on a wgpu compute shader. After 1 s
of simulated time (1000 RK4 steps × 17 species × 11 enzymes per derivative
call), every species agrees with the f64 CPU reference to **0.0000% relative
error to four decimal places**. Worst-fit species: Pyruvate at 0.0002%
relative error — well below the 1% Phase 11 acceptance threshold.

```
species      cpu (mM)       gpu (mM)   rel-err
Glucose      4.984250       4.984253    0.0000
G6P          0.081440       0.081440    0.0000
F6P          0.006359       0.006359    0.0000
F1,6BP       0.002920       0.002920    0.0000
DHAP         0.070583       0.070583    0.0000
GAP          0.003023       0.003023    0.0000
1,3-BPG      0.000169       0.000169    0.0000
3-PG         0.180578       0.180578    0.0000
2-PG         0.004911       0.004911    0.0000
PEP          0.004858       0.004858    0.0000
Pyruvate     0.038512       0.038512    0.0000
Lactate      1.510700       1.510700    0.0000
ATP          2.028105       2.028102    0.0000
ADP          0.221895       0.221895    0.0000
NAD+         0.099472       0.099472    0.0000
NADH         0.000528       0.000528    0.0000
Pi           0.949964       0.949964    0.0000
```

This validates the toolchain end-to-end on the most algorithmically
substantive deliverable in Phase 11.

## What landed

### Files added

- `shaders/compute/biochem_glycolysis.wgsl` (~280 lines):
  - 11 enzyme rate functions ported from `src/biochemistry/glycolysis.rs`
    (HK, GPI, PFK, Aldolase, TPI, GAPDH, PGK, PGM, Enolase, PK, LDH).
  - 5 kinetic helpers (Michaelis-Menten, Hill, reversible MM,
    ordered bi-bi, reversible ordered bi-bi).
  - `derivatives()` function that mirrors `MetabolismSolver`'s closure:
    enzyme rates → stoichiometric scatter to dydt → world-level
    GLUT1 / MCT1 / external ATP demand.
  - `rk4_step()` compute entry: one workgroup-thread per cell, runs
    `n_steps` RK4 ms-steps in registers without buffer round-trips.
- `src/compute/biochem.rs` (~250 lines): `GlycolysisBatchConfig` and
  `run_glycolysis_batch()` host wrapper. Handles f64↔f32 conversion at
  upload/download.
- `tests/biochem_gpu_parity.rs` (~150 lines): 1-second CPU vs GPU
  simulation with per-species 1% / 1e-3 mM tolerance.

### Files modified

- `src/compute/mod.rs` — re-export `run_glycolysis_batch`,
  `GlycolysisBatchConfig`.

## Sub-phase scope (11.2.A)

This is the glycolytic *backbone* (gluconeogenesis-side metabolism not
included in any phase yet). Specifically:

- **Included**: 11 glycolysis enzymes, GLUT1 facilitated diffusion,
  MCT1 lactate export, external ATP consumption (membrane pumps).
- **Excluded** (deferred to 11.2.B): 2,3-BPG shunt (BPGM, BPGP).
- **Excluded** (deferred to 11.2.C): PPP (G6PDH, 6PGDH, RPE, RPI, TK,
  TA, PGL), glutathione cycle (GR, GPx, gamma-GCS, GS), Piezo1, ion
  homeostasis (Na/K-ATPase, leaks).
- **Excluded** (deferred to 11.2.D): Hemoglobin Adair binding + Bohr/DPG
  modulation (post-RK4 Euler), pH buffer (algebraic).
- **Excluded** (deferred to 11.2.E): three inline homeostasis hacks in
  `full_integration.rs::step` (basal NADPH consumption, basal GSH
  oxidation, ATP regen).

The 17 species ported in 11.2.A are the first 17 indices of the
38-species `FullyIntegratedIndices`. The full kernel for 11.2.B–E is a
mechanical extension of the same template.

## Numerical strategy

Per the Phase 11 plan: f32 throughout the GPU kernel, f64 on CPU,
conversion at host boundaries. Glycolytic species span 0.001–5 mM, well
within f32 precision (24-bit mantissa, ~7 sig figs). Results vindicate
the choice — bit-near-perfect agreement at 4 decimals.

The Kahan summation escape hatch was *not* needed for glycolysis. It
remains reserved for the lower-concentration species in 11.2.C (Ca²⁺ at
1e-4 mM, H₂O₂ at ~1 µM, NADP⁺ at ~0.002 mM) where f32 precision could
matter.

## Architecture

`run_glycolysis_batch` is a one-call dispatch: take N pools, upload, run,
download. It recreates the pipeline + bind group layout per call which
is wasteful for a hot loop but appropriate for Phase 11.2.A where the
kernel API is still settling. A long-lived `GlycolysisBackend` struct
that caches the pipeline will land in 11.2.B once the kernel surface
stops changing.

## Tests added

- `tests/biochem_gpu_parity.rs::glycolysis_cpu_gpu_parity_one_second`:
  builds a CPU "glycolysis-only" RK4 loop (mirrors the GPU kernel
  exactly so the comparison is not confounded by the 2,3-BPG shunt that
  `MetabolismSolver` includes) and asserts per-species 1% / 1e-3 mM
  match. Skipped if no GPU adapter is available.

## Test counts after 11.2.A

- 186 lib tests (unchanged from 11.0).
- 131 integration tests (was 130; +1 from `biochem_gpu_parity`).
- Validation suite: 7/8 (unchanged).

## What's next (Phase 11.2.B)

Add the 2,3-BPG shunt (BPGM, BPGP — 2 enzymes from
`src/biochemistry/rapoport_luebering.rs`) and the species at index 17.
This is small (~50 lines of WGSL) and demonstrates that the kernel
extension pattern works. After 11.2.B the simulator's GPU path matches
`MetabolismSolver` exactly, which is the natural unit-of-comparison.

After that, 11.2.C (PPP + glutathione + Piezo1 + ions, ~300 lines of
WGSL) and 11.2.D (hemoglobin + pH, ~80 lines + post-RK4 Euler) close
out the 38-species port.

## Why Phase 11.1 was skipped

The Phase 11 plan called for an SoA refactor (`PhysicsState` →
parallel `Vec<f32>`, `MetabolitePoolBatch` species-major) before any
GPU work. Phase 11.2.A skipped that refactor entirely: the host wrapper
performs AoS→flat conversion *at upload time* rather than restructuring
the CPU-side data permanently. For multi-cell biochem this costs a
single O(N×17) gather per dispatch — negligible compared to the
dispatch latency. The structural refactor remains valuable if the GPU
path becomes the inner hot loop, but it's not load-bearing for a
correct first port.
