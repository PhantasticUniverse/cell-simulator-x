# Phase 11.2.E Notes — Inline Homeostasis Corrections + Full-Solver Parity

Date: 2026-05-01

## Headline

The three inline homeostasis correction terms from
[full_integration.rs:251-294](../src/biochemistry/full_integration.rs)
are now ported to the GPU kernel, closing the last gap between
`FullyIntegratedSolver::step` (CPU) and the GPU compute kernel. After 1 s
of simulated time with the full production parameter set, every species
agrees with CPU to **0.0023% worst-fit relative error** (Na⁺) and Hb
saturation to **2e-6 absolute** (0.937480 vs 0.937478 from initial 0.75).

This closes Phase 11.2 (biochemistry on GPU). The 38-species ODE solver
including hemoglobin and pH buffer is now a single GPU compute pass that
matches the CPU production code path bit-near-perfectly.

## What landed

### Files modified

- `shaders/compute/biochem_full.wgsl`:
  - New uniform `enable_homeostasis_corrections: u32` (Phase 11.2.E
    section in the U struct).
  - In `derivatives`, gated by the flag, three additions:
    1. **Basal NADPH consumption** — methemoglobin reductase, thioredoxin,
       cytochrome b5 reductase, etc. Coefficient 0.0003 with NADPH/NADP⁺
       ratio modulation capped at 30.
    2. **Basal GSH oxidation** — non-enzymatic, protein-S-glutathionylation.
       Coefficient 0.001 with GSH/GSSG ratio modulation capped at 5×.
    3. **ATP-deficit regen** — pulls ATP toward 1.8 mM via mass-action on
       ADP. Coefficient 0.2.
  - Pad realigned: 42 fields × 4 = 168 bytes + 8 bytes pad = 176 bytes
    (multiple of 16, wgpu requirement).

- `src/compute/biochem.rs`:
  - `FullBiochemUniforms` extended with `enable_homeostasis_corrections`.
  - `FullBiochemBatchConfig` extended with the same field, defaulting to
    `false`. The 11.2.C and 11.2.D parity tests stay bit-stable because
    they don't enable corrections.

### Files added

- `tests/biochem_gpu_full_solver_parity.rs` — the strictest possible
  parity test. CPU side runs the production `FullyIntegratedSolver::step`
  for 1000 1-ms steps. GPU side runs the same number of steps with all
  flags on (Piezo1 + ions + Hb + pH + corrections). Asserts per-species
  ≤ 1% rel-err **and** `|Δsaturation| < 0.01`.

## Numerical fidelity

The CPU and GPU paths now share the entire stack:

| Layer | CPU | GPU | Match |
|-------|-----|-----|-------|
| Glycolysis (11 enzymes) | `GlycolysisSolver` | `rate_*` in WGSL | Verbatim |
| 2,3-BPG shunt (BPGM, BPGP) | `RapoportLueberingSolver` | `rate_bpg{m,p}` | Verbatim |
| PPP (7 enzymes) | `PentosePhosphatePathway` | `rate_g6pdh`–`rate_ta` | Verbatim |
| Glutathione (4 enzymes + H₂O₂) | `GlutathioneCycle` | `rate_g{px,r,amma_gcs,s}` | Verbatim |
| Piezo1 + PMCA + Pannexin-1 | `Piezo1System` | `piezo1_*`, `pmca_*`, `atp_release_*` | Verbatim |
| Na/K-ATPase + leaks | `IonHomeostasisSystem` | `nakatpase_rate` + leak terms | Verbatim |
| Hb Adair + Bohr/DPG/T/CO₂ | `HemoglobinSolver::step` | `hb_post_step` | Verbatim |
| pH buffer | `PhBufferModel::calculate_ph` | `compute_ph` | Verbatim |
| Basal NADPH brake | `step` line 251 | derivatives() line 654 | Verbatim |
| Basal GSH brake | `step` line 266 | derivatives() line 660 | Verbatim |
| ATP regen | `step` line 276 | derivatives() line 666 | Verbatim |
| RK4 integrator | `RK4Integrator::step` | `rk4_step` | Verbatim |
| Per-step max-change clamp + min-conc floor | `IntegratorConfig` | `apply_clamp_and_floor` | Verbatim |

f32 throughout the GPU vs f64 throughout the CPU. After 1 s of integration
the worst-fit species (Na⁺) is 0.0023% off — well under the 1% gate. f32
held perfectly across the full 10⁷× dynamic range from 22 µM (6-PGL) to
140 mM (K⁺). Kahan summation was reserved as an escape hatch in the Phase
11 plan but never needed.

## Test counts after 11.2.E

- 186 lib tests (unchanged).
- 134 integration tests (was 133; +1 from `biochem_gpu_full_solver_parity`).
- Validation suite: 7/8 (unchanged).

## Phase 11.2 retrospective

Phase 11.2 unfolded as five sub-phases of growing scope, each one a tested
delta on the previous:

| Sub | Scope | Species | Worst-fit rel-err | Lines added |
|-----|-------|---------|--------------------|-------------|
| A | Glycolysis (11 enzymes) | 17 | 0.0002% (Pyruvate) | ~480 (kernel + host + test) |
| B | + 2,3-BPG shunt (BPGM/BPGP) | 18 | 0.0009% (2,3-BPG) | +60 |
| C | + PPP, glutathione, Piezo1, ions | 38 | 0.0023% (Na⁺) | +700 |
| D | + Hb Adair + pH buffer (post-RK4 Euler) | 38 + sat | 2e-6 abs (Hb sat) | +250 |
| E | + 3 inline homeostasis corrections | 38 + sat | 0.0023% (Na⁺) | +120 |

The single biggest discovery — Phase 10's aldolase-Keq unit-error fix —
was carried into the GPU kernel with literal numeric values
(`mass_action / (fbp * 7.8e-2)` in WGSL), so the parity tests would have
caught any regression on either side.

Cumulatively, Phase 11.2 added one WGSL kernel (~680 lines), two host
wrapper functions, and four integration tests. Each test takes < 0.3 s
to run on Apple M4 Max.

## What's next (Phase 11.3)

The hard port: Skalak fused with TensionComputer, WLC with vertex-centric
CSR adjacency (replacing the HashMap aggregation in
[wlc.rs](../src/physics/wlc.rs)), DPD membrane forces with stateless PCG
hash for RNG, Velocity-Verlet position update. All dispatched inside a
single `ComputePass` per substep.

The plan called this 3 weeks of work because the data-layout refactor
(SoA, vertex-centric CSR) blocks the actual port. Some of the SoA
groundwork landed in Phase 10.5 already; the rest needs to happen as
prep for 11.3.

The kill-switch for Phase 11.3 is the speedup gate: ≥ 10× wall-clock at
N=100 vs CPU baseline. If achieved < 3×, descope to "GPU biochem only,
physics stays on CPU" and ship Phase 11.2 as the Phase 11 deliverable.
The biochemistry kernel as it stands in 11.2.E is already the more
algorithmically substantive deliverable, so this descope is acceptable
as a fallback.
