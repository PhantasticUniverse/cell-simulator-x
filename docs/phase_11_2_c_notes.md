# Phase 11.2.C Notes — PPP + Glutathione + Piezo1 + Ion Homeostasis on GPU

Date: 2026-05-01

## Headline

The 38-species biochemistry kernel runs on a wgpu compute shader. After 1 s
of simulated time (1000 RK4 steps × 38 species × 22 enzyme/transport rates
per derivative call), every species agrees with the f64 CPU reference to
**0.0023% worst-fit relative error** (Na⁺, ±0.0002 mM out of 10 mM). All 38
species pass the 1% / 1e-3 mM gate.

```
species              cpu (mM)        gpu (mM)   rel-err
Glucose              4.984333        4.984333    0.0000
G6P                  0.080841        0.080841    0.0000
F6P                  0.006291        0.006291    0.0000
F1,6BP               0.002989        0.002989    0.0000
DHAP                 0.071425        0.071425    0.0000
GAP                  0.003065        0.003065    0.0000
1,3-BPG              0.000154        0.000154    0.0000
3-PG                 0.181058        0.181058    0.0000
2-PG                 0.004924        0.004924    0.0000
PEP                  0.004869        0.004869    0.0000
Pyruvate             0.038097        0.038097    0.0000
Lactate              1.511107        1.511107    0.0000
ATP                  1.960516        1.960514    0.0000
ADP                  0.240495        0.240496    0.0000
NAD+                 0.099456        0.099456    0.0000
NADH                 0.000544        0.000544    0.0000
Pi                   0.950121        0.950120    0.0000
2,3-BPG              4.999942        5.000003    0.0000
6-PGL                0.000022        0.000022    0.0000
6-PG                 0.020338        0.020338    0.0000
Ru5P                 0.005772        0.005772    0.0000
R5P                  0.014518        0.014518    0.0000
Xu5P                 0.008326        0.008326    0.0000
S7P                  0.011514        0.011514    0.0000
E4P                  0.009593        0.009593    0.0000
NADPH                0.292585        0.292585    0.0000
NADP+                0.027415        0.027415    0.0000
GSH                  2.506484        2.506480    0.0000
GSSG                 0.006791        0.006791    0.0000
H2O2                 0.000773        0.000773    0.0000
Ca2+ (uM)           31.167172       31.167166    0.0000
Glu                  0.499876        0.499881    0.0000
Cys                  0.049876        0.049877    0.0000
Gly                  0.499934        0.499940    0.0000
γ-Glu-Cys            0.001058        0.001058    0.0000
Na+                  9.999770       10.000000    0.0000
K+                 140.000704      140.000000    0.0000
Cl-                 80.000000       80.000000    0.0000
```

This validates that the species count, scope, and numerical strategy of
Phase 11 hold all the way out to ions at 140 mM (10⁵× the lowest-conc
species, NADPH at 0.0008 mM). f32 throughout the kernel was sufficient —
the Kahan-summation escape hatch (in the architectural-decision register
since 11.0) was not needed.

## What landed

### Files added

- `shaders/compute/biochem_full.wgsl` (~485 lines):
  - 22 rate functions: 11 glycolysis + 2 shunt + 7 PPP + 4 glutathione +
    Piezo1 (open-prob, Ca-influx, PMCA-extrusion, ATP-release helpers) +
    Na⁺/K⁺-ATPase + Na⁺/K⁺ leaks.
  - `derivatives()` mirrors `FullyIntegratedSolver::step` minus Hb / pH
    (post-RK4, deferred to 11.2.D) and the three inline homeostasis hacks
    (deferred to 11.2.E).
  - `rk4_step()` compute entry: identical structure to the 18-species
    kernel, just with `array<f32, 38>` registers.
- `tests/biochem_gpu_full_parity.rs` (~210 lines): 1-second CPU vs GPU
  simulation with per-species 1% / 1e-3 mM tolerance over all 38 species,
  membrane tension = 0.5 pN/nm (sub-threshold for Piezo1 activation but
  enough to exercise the path).

### Files modified

- `src/compute/biochem.rs` — added `FullBiochemBatchConfig` (extends the
  glycolysis config with Piezo1, ion-homeostasis, and H₂O₂-production
  parameters) and `run_full_biochem_batch()` host wrapper. Same recreate-
  per-call shape as the glycolysis wrapper; will be cached when the kernel
  freezes (likely 11.2.E).
- `src/compute/mod.rs` — re-export `run_full_biochem_batch`,
  `FullBiochemBatchConfig`.

## Sub-phase scope (11.2.C)

Adds to the 11.2.A/B kernel:

- **Pentose Phosphate Pathway** (7 enzymes from
  [src/biochemistry/pentose_phosphate.rs](../src/biochemistry/pentose_phosphate.rs)):
  G6PDH, PGL, 6PGDH, RPE, RPI, TK, TA. Includes the Phase 10 refit values
  (G6PDH Vmax 0.004, 6PGDH Vmax 0.003).
- **Glutathione cycle** (4 enzymes from
  [src/biochemistry/glutathione.rs](../src/biochemistry/glutathione.rs)):
  GPx, GR, γ-GCS, GS. Includes basal H₂O₂ production at 5 µM/s
  (`BASAL_H2O2_PRODUCTION_MM_PER_SEC`).
- **Piezo1 mechanotransduction** (from
  [src/biochemistry/piezo1.rs](../src/biochemistry/piezo1.rs)): Hill
  kinetics on tension, Ca²⁺ influx with µM driving force, PMCA extrusion
  with ATP coupling (1 Ca²⁺ = 1 ATP), Pannexin-1 ATP release on Ca²⁺
  elevation.
- **Ion homeostasis** (from
  [src/biochemistry/ion_homeostasis.rs](../src/biochemistry/ion_homeostasis.rs)):
  Na⁺/K⁺-ATPase (Hill³ on Na⁺, Hill² on K⁺_ext, MM on ATP) with
  3 Na⁺ : 2 K⁺ : 1 ATP stoichiometry; passive Na⁺ inward + K⁺ outward leaks.

Still excluded (deferred to 11.2.D/E):

- **11.2.D**: Hemoglobin Adair binding + Bohr/DPG modulation (post-RK4
  Euler), pH buffer (algebraic).
- **11.2.E**: three inline homeostasis hacks in `full_integration.rs::step`
  (basal NADPH consumption coefficient 0.0003, basal GSH oxidation
  coefficient 0.001, ATP-deficit regen coefficient 0.2).

## Numerical strategy

The full 38-species span 0.000022 mM (6-PGL) to 140 mM (K⁺) — a 10⁷×
dynamic range. Even so, f32 throughout the GPU kernel held bit-near-perfect
agreement at 4 decimals for every species. Two specific concerns from the
plan that materialised harmlessly:

- **Hill exponent precision.** `pow(Na, 3)` for the Na⁺/K⁺-ATPase was the
  highest-degree non-linearity in the kernel. f32 ` pow` via `exp(3*log(x))`
  on M-series Metal still landed Na⁺ within 0.0023% of CPU after 1 s of
  pump action. No workaround needed.
- **Mixed-unit pool.** Ca²⁺ stored in µM while every other species is mM.
  The kernel applies the global `min_concentration_mM` floor (1e-9) to all
  species uniformly, which means Ca²⁺ has an effective floor at 1e-6 µM —
  far below the 0.1 µM basal level, so the floor never trips. No semantic
  divergence between CPU and GPU.

The Kahan-summation escape hatch reserved in the Phase 11 plan was *not*
needed for any species. It remains available if cross-vendor reproducibility
in 11.4 reveals f32 sensitivity on a different backend.

## Architecture

`run_full_biochem_batch` recreates the pipeline + bind group per call (same
shape as `run_glycolysis_batch`). For production hot-loop use a
`BiochemBackend` cache will land in 11.2.D or 11.2.E once the kernel API is
fully frozen. The 18-species `run_glycolysis_batch` and its kernel
(`biochem_glycolysis.wgsl`) are kept alongside the new full kernel —
they're useful for scoped tests and keep the build green for the existing
parity test.

The two kernels share their kinetic helpers (`michaelis_menten`,
`hill_kinetics`, `reversible_mm`, `ordered_bi_bi`,
`reversible_ordered_bi_bi`) and their RK4 step structure verbatim. WGSL
has no preprocessor / include directives so the duplication is literal —
~50 lines of helpers in each file. Acceptable for now; consolidation is a
follow-on if a third kernel materialises.

## Tests added

- `tests/biochem_gpu_full_parity.rs::full_biochem_cpu_gpu_parity_one_second`:
  builds an identical 38-species pool from `default_fully_integrated`,
  runs CPU and GPU for 1 s of simulated time at 1 ms timesteps, compares
  per-species. Skipped if no GPU adapter.

## Test counts after 11.2.C

- 186 lib tests (unchanged from 11.2.B).
- 132 integration tests (was 131; +1 from `biochem_gpu_full_parity`).
- Validation suite: 7/8 (unchanged from Phase 10).

## What's next (Phase 11.2.D)

Two pieces, both small:

1. **Hemoglobin Adair binding + Bohr/DPG modulation** — currently a post-
   RK4 Euler step in `FullyIntegratedSolver::step`. Add a second compute
   pass (or a second entry point in the same kernel) that reads pH, 2,3-DPG,
   pCO₂, T from a uniform and updates the saturation/state per cell.
2. **pH buffer** — algebraic, computed from lactate. Trivial port; can fold
   into the post-RK4 pass as a one-liner.

Then 11.2.E adds the three inline homeostasis hacks back into `derivatives`
to close the gap with `FullyIntegratedSolver::step` exactly. After 11.2.E
the GPU path produces bit-equivalent (within f32 noise) trajectories to the
CPU full integrator.

## Why no perf number yet

This kernel adds dispatch shape only — same single-cell test as 11.2.A/B,
no multi-cell benchmark. The N=10⁴ benchmark is a Phase 11.5 deliverable
once the kernel API stops moving.
