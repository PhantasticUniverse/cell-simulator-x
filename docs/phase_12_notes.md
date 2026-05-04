# Phase 12 Notes — External Fluid (Single-Cell Capillary)

Date: 2026-05-02

## Headline (Phase 12.A)

A new `flow` module implements analytic Poiseuille flow in a cylindrical
channel and a Stokes-style drag coupling that pushes per-vertex drag
forces into `PhysicsState::external_forces_uN`. The existing CPU
`PhysicsSolver::step` picks them up natively, so a cell placed in a
channel with non-zero centerline velocity drifts along the flow axis.

This is the first scientific capability the simulator gains beyond the
"baseline RBC + biochemistry" framework — the start of the journey
toward parachute shape under flow (Phase 12.C) and tank-treading
frequency validation against Fischer 2007.

## What landed (Phase 12.A)

### `src/flow/mod.rs` and `src/flow/poiseuille.rs`

- `CylindricalChannel { radius_um, axis, center }` — geometry.
- `Poiseuille { channel, max_velocity_um_per_sec }` — analytic field
  with `velocity_at(p) = v_max * (1 - (r/R)²) * axis` inside the channel,
  zero outside.
- `Poiseuille::mean_velocity_um_per_sec()` and
  `Poiseuille::wall_shear_rate_per_sec()` — derived quantities for
  validation.
- `drag_force_uN(flow, pos, vel, drag_coeff) = drag_coeff * (v_fluid - v_vertex)`
  — Stokes-form drag.
- `apply_drag_to_external_forces(state, mesh, flow, drag_coeff)` —
  populates `state.external_forces_uN` so the existing CPU solver picks
  it up unchanged.

### Tests

- 8 unit tests in `flow::poiseuille::tests` covering profile shape,
  mean / wall-shear derivations, drag direction, and external-force
  population.
- 2 integration tests in `tests/flow_deformation.rs`:
  - `cell_translates_along_flow_axis` — cell at channel center drifts
    +0.018 µm along Z over 1000 substeps with 1 mm/s centerline flow.
    Lateral drift stays < 0.05 µm.
  - `no_flow_no_translation` — zero centerline velocity → centroid
    drift < 1e-3 µm.

### Counts

- 197 lib tests (was 189; +8 from `flow::poiseuille::tests`).
- 142 integration tests (was 140; +2 from `flow_deformation`).

## What's deferred

### Phase 12.B — `PhysicsBackend` integration

The CPU path uses `PhysicsState::external_forces_uN`, which the existing
`PhysicsSolver::step` already integrates. The GPU `PhysicsBackend` does
**not** yet have an external-force input buffer. Adding it requires:

1. A new storage buffer `external_forces` bound at the kernel level.
2. An update to `skalak_init_from_baseline` (or a new kernel) that adds
   `external_forces[v]` to `forces[v]` alongside the WLC baseline.
3. A host-side helper that updates `external_forces` between substeps
   (since drag depends on vertex velocity, which changes each step).

The kernel change is small. Defer to 12.B because it requires another
parity test against the CPU path with flow.

### Phase 12.C — Quantitative validation

Two canonical benchmarks:

- **Parachute shape under Poiseuille flow.** A cell of radius ≈ 4 µm in
  a channel of radius ≈ 5 µm at moderate flow adopts a "parachute"
  profile (concave on the upstream side). Compare equilibrium shape to
  Skalak 1973 / Fung 1993 figures and Pivkin–Karniadakis 2007.
- **Tank-treading frequency vs shear rate** (Fischer 2007). At fixed
  high shear rate, the membrane "tank-treads" — rotates around the
  cell's center of mass — at a frequency proportional to shear rate.
  Track the angular position of an identifiable vertex over time;
  measure rotation rate; compare to Fischer's 0.1–0.4 × γ̇ band.

Both require longer simulations and post-processing. Phase 12.C will
likely produce a `--diagnose-tank-treading` CLI mode and a notebook /
script that runs a sweep of shear rates.

## What's still unknown / open questions

- **Drag coefficient calibration.** The default `drag_coeff` in tests is
  hand-tuned (5 µN·s/µm) so drag dominates membrane tension at
  reasonable flow speeds. A first-principles value from Stokes drag on a
  small effective sphere comes out ~1.1 µN·s/µm — close, but
  per-vertex distribution depends on the local membrane normal and
  curvature, which we ignore. Phase 12.C will pin this down.
- **No-slip and inlet/outlet boundary conditions.** The current channel
  is infinite and the cell never reaches a wall. For tighter geometries
  (e.g., the splenic 0.5 µm interendothelial slit envisioned for Phase
  13's microcirculation), we'll need a wall repulsion force.
- **No two-way coupling.** The fluid pushes the cell, but the cell does
  not perturb the fluid. For a small RBC in a much larger vessel that's
  fine; for capillaries (vessel ≈ cell size) we'd need IBM or LBM. Defer
  to Phase 13 / future work.

## Rationale for choosing Poiseuille first

The plan listed "immersed-boundary or simplified Stokes" — Poiseuille is
simpler than either, but unlocks the same scientific question: how does
the membrane respond to a known velocity gradient? Tank-treading and
parachute shape both have analytic Poiseuille drivers in the literature
(Fischer 2007 used a rheoscope with quasi-Couette flow; parachute
papers use Poiseuille). Starting analytic keeps Phase 12.A small while
delivering real scientific capability — and decouples the
fluid-solver complexity from the membrane-response physics that's the
actual headline.

## What's next

If staying sequential: Phase 12.B (`PhysicsBackend` integration) →
Phase 12.C (validation) → Phase 13 (multi-cell microcirculation) →
**Phase 14 (storage lesion at 42-day timescale).**

If pivoting toward the headline: Phase 14 directly. Storage lesion is
single-cell, biochemistry-dominated, and doesn't need flow. Phase 12.B
+ 12.C are valuable for the broader paper but not gating for the
headline metabolomics curves (D'Alessandro et al., Rolfsson et al.,
storage time-courses).

---

## Phase 12.B: PhysicsBackend external_forces buffer + World drag helper

Phase 14 has shipped (storage lesion + 14.D extensions). Phase 12.B now
closes the GPU side of the flow story so single-cell-in-flow scenarios
can be GPU-accelerated.

### 12.B.1 — PhysicsBackend external_forces buffer

`shaders/compute/physics_step.wgsl` gains a 9th storage buffer (binding
9, read-only, n_vertices*3*f32):
```
@group(0) @binding(9) var<storage, read> external_forces: array<f32>;
```
Inside `skalak_init_from_baseline`, each vertex's force aggregation now
reads `external_forces[v]` alongside `wlc_baseline[v]` and the sum of
its CSR-aggregated Skalak per-element contributions. Default-zero
buffer makes the new code path a no-op until the host writes drag.

`PhysicsBackend` exposes `set_external_forces(&[Vec3])` and
`clear_external_forces()`. The buffer is COPY_DST so
`queue.write_buffer` works without staging.

`ComputeContext::new_headless` requests
`max_storage_buffers_per_shader_stage = 16` (wgpu default 8 was the
Vulkan mobile baseline; Apple Silicon, NVIDIA, AMD all support ≥16).

#### Validation

`tests/physics_external_forces_parity.rs` (2 tests):
- `external_forces_uniform_constant_parity`: GPU under uniform +z drag
  (1 nN/vertex) matches CPU reference within 1e-5 μm after 100
  substeps.
- `external_forces_zero_matches_baseline_parity`: zero-buffer
  regression-safety check.
All 6 prior physics GPU parity tests still pass.

### 12.B.2 — World::apply_poiseuille_drag

`World::apply_poiseuille_drag(&Poiseuille, drag_coeff)` runs
`flow::apply_drag_to_external_forces` per-cell in parallel via rayon.
The CPU `PhysicsSolver::step` already reads
`physics_state.external_forces_uN`, so this completes the CPU path
end-to-end.

`tests/flow_cpu_gpu_parity.rs::poiseuille_drag_cpu_gpu_parity` proves
the GPU path matches CPU within tolerance: per substep, drag is
computed host-side from the GPU's read-back position+velocity and
uploaded via `set_external_forces`. After 100 substeps:
- per-vertex |Δpos| < 1e-5 μm
- cell-centroid drift < 1e-3 μm (the Phase 12.B plan target)

### Deferred to 12.B.3 (post-12.C)

- Per-substep drag computation entirely on GPU (currently the host
  reads back positions/velocities each substep). Acceptable for now
  because PhysicsBackend already pays for one read-back per substep.
  A pure-GPU drag kernel becomes attractive once Phase 13 multi-cell
  scaling pushes that read-back over budget.
- Per-cell PhysicsBackend integration in `World` (currently World only
  uses the CPU PhysicsSolver). The `flow_cpu_gpu_parity` test exercises
  the GPU path at the PhysicsBackend layer; wiring it through
  `World::step` is a Phase 13 architectural change since N>1 cells need
  shared-buffer coordination.
