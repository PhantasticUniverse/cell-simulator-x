# Phase 11.3 Notes — Physics on GPU (Skalak + WLC + DPD + Verlet)

Date: 2026-05-01

## Headline

The four physics force/integrator components — `SkalakSolver::compute_forces`,
`WLCSolver::compute_forces`, `DPDSolver::compute_membrane_forces`, and
`VelocityVerlet::half_step_velocity` + `step_position` — are all ported
onto wgpu compute kernels with bit-near-perfect CPU/GPU parity.

| Sub | Kernel | CPU/GPU agreement | Algorithmic novelty |
|-----|--------|-------------------|----------------------|
| A | `skalak.wgsl` (per-element + per-vertex CSR aggregate) | 0.0127% rel-err worst-fit | Replaces serial scatter with deterministic CSR |
| B | `wlc.wgsl` (per-edge + per-node CSR aggregate) | 1e-11 μN abs worst-fit | Replaces HashMap aggregation with CSR |
| C | `dpd.wgsl` (per-vertex with stateless PCG hash RNG) | 5e-8 μN abs worst-fit | Replaces `StdRng` with bit-deterministic PCG |
| D | `integrator.wgsl` (per-vertex Verlet) | 0 absolute (bit-identical) | Direct port |

All four kernels run on Apple M4 Max via Metal, on top of wgpu 0.20.

## What landed

### Files added

- `shaders/compute/skalak.wgsl` (~145 lines): two-kernel pipeline (per-
  element force compute + per-vertex CSR aggregation).
- `shaders/compute/wlc.wgsl` (~110 lines): same pattern for spectrin edges/
  nodes; replaces the CPU-side HashMap aggregation that was the documented
  GPU blocker.
- `shaders/compute/dpd.wgsl` (~80 lines): per-vertex DPD with stateless
  PCG hash + Box-Muller for deterministic noise.
- `shaders/compute/integrator.wgsl` (~70 lines): two compute entries for
  half-step velocity and position update.
- `src/compute/physics.rs` (~830 lines): host wrappers + Pod uniforms +
  CSR builders for all four kernels.
- `tests/skalak_gpu_parity.rs`, `tests/wlc_gpu_parity.rs`,
  `tests/dpd_gpu_parity.rs`, `tests/verlet_gpu_parity.rs`.

### Files modified

- `src/physics/dpd.rs` — added `pcg_hash_u32`, `pcg_uniform`,
  `pcg_gaussian` free helpers and a new `compute_membrane_forces_pcg`
  method that replaces `StdRng` with the same stateless PCG used in WGSL.
- `src/compute/mod.rs` — re-exports for the new public items.

## Architectural decisions

### Two-pass scatter via CSR (not atomic float add)

WGSL has no native f32 atomic add. The two options were:

1. **Atomic CAS loops on `atomic<u32>` with bit-cast** — works, but breaks
   bit-identity due to non-deterministic ordering and adds overhead.
2. **Pre-bin per-vertex CSR adjacency** — each thread writes once to its
   own per-element slot, then a per-vertex aggregation pass sums.

We chose option 2. CSR is built once per mesh (and again per spectrin
network) at host time. The aggregation order is fixed (CSR offset → end),
so f32 rounding is deterministic — given the same input, we get the same
output every run.

This pattern repeats verbatim across Skalak (per-element-per-vertex) and
WLC (per-edge-per-node).

### Stateless PCG hash for DPD RNG

`StdRng::from_entropy()` made the existing CPU code non-reproducible run-
to-run. The plan called for stateless PCG hash keyed on `(vertex_id,
step_count, axis, seed)`. We implemented it on both CPU and GPU with
bit-identical math, so:

- Existing legacy code paths (using `StdRng`) keep working.
- New `compute_membrane_forces_pcg` is the deterministic path.
- The GPU kernel uses the same hash + Box-Muller, producing 5e-8 μN
  absolute agreement after a single dispatch (i.e. f32 epsilon, attributable
  only to driver-level differences in `sin`/`cos`/`log`/`sqrt`).

### One-shot dispatch (per call) — persistent backends deferred

Every Phase 11.3 host wrapper currently allocates buffers, builds the
pipeline + bind group, dispatches, downloads, and frees per call. This is
appropriate while the kernel surface is still settling but **wasteful for
the substep hot loop** (1000 substeps/sec at typical settings).

Persistent `PhysicsBackend` caching is deferred to Phase 11.3.E together
with the multi-substep mega-pass. The kernels themselves are stable and
won't change shape; the host wrapper layer is what gets the cache.

## Phase 11.3.E — Integrated pipeline (closed)

The integration layer landed. New module: `PhysicsBackend` (in
`src/compute/physics.rs`) owns:

- All persistent buffers (positions, velocities, forces, mesh elements,
  CSR adjacency, WLC baseline, scratch).
- All compiled compute pipelines.
- A single `step()` method that dispatches six kernels in one
  `ComputePass`, mirroring `PhysicsSolver::step` exactly.

### How the six dispatches map to the CPU pipeline

| # | Kernel | CPU equivalent | Threads |
|---|--------|----------------|---------|
| 1 | `verlet_half_step_velocity` | `integrator.half_step_velocity(forces)` | n_vertices |
| 2 | `verlet_step_position` | `integrator.step_position(velocities)` | n_vertices |
| 3 | `skalak_per_element` | `skalak.compute_forces` (inner loop) | n_elements |
| 4 | `skalak_init_from_baseline` | scatter + WLC contribution | n_vertices |
| 5 | `dpd_add` | `dpd.compute_membrane_forces_pcg` | n_vertices |
| 6 | `damping_and_half_step` | `forces -= vel*damping` + final half | n_vertices |

### WLC is a static baseline, not a per-substep dispatch

Looking at `WLCSolver::compute_forces`, it reads `network.nodes[i].position`
which never updates as the mesh moves (the spectrin network struct is
borrowed `&immutable` throughout the simulation). WLC therefore contributes
a **constant force per substep** — equal to its initial value. Phase 11.3.E
exploits this: at backend creation we run the existing `run_wlc_forces`
once, project the per-spectrin-node forces into per-mesh-vertex
contributions via `mesh_vertex_idx`, and store the result in a
`wlc_baseline` storage buffer. Every substep folds it back in by having
`skalak_init_from_baseline` initialize `forces[v] = wlc_baseline[v]` before
adding the Skalak aggregate. No per-substep WLC kernel is required.

(If a future change moves spectrin nodes with the mesh, this baseline
buffer becomes per-substep instead — but as the CPU code stands today,
the static baseline matches CPU semantics exactly.)

### Parity result

`tests/physics_backend_parity.rs::physics_backend_cpu_gpu_parity` runs
the full CPU `PhysicsSolver::step`-equivalent loop (with PCG-DPD on both
sides for determinism) and the GPU `PhysicsBackend.step()` for 10 substeps
each. Worst |Δpos| = **1.9e-15 μm** (essentially f64 noise on near-zero
y-coordinates; non-trivial axes match bit-identically), worst |Δvel| =
**8.3e-10 μm/s** — well below the 1e-3 threshold.

This is the integrated physics pipeline running end-to-end on the GPU
with bit-near-perfect agreement against CPU. From here, Phase 11.3.F (or
the eventual N-cell scaling phase) can drive multiple substeps via a
host-side `for _ in 0..N { backend.step() }` loop with no readback.

## What's still on the floor

- **Persistent multi-substep dispatch** — currently the host calls
  `backend.step()` N times. Each call submits its own command encoder,
  which has ~µs of dispatch latency. For 1000 substeps/sec at typical
  settings that's ~1 ms of overhead — negligible at single-cell N=1, but
  worth measuring at multi-cell scale.
- **Multi-cell support** — `PhysicsBackend` is single-cell (one mesh,
  one spectrin network). Multi-cell scaling requires either many backends
  or a single backend with an outer "cell" axis on every buffer.
- **Render-compute buffer sharing** (Phase 11.4) — currently the only way
  to visualize is `read_state()` followed by uploading positions to the
  render path. Direct buffer sharing eliminates that round-trip.
- **Mesh normals recomputation** — the CPU side calls `recompute_normals`
  after `step_position` for diagnostic / rendering use. The GPU pipeline
  doesn't yet update normals; positions are correct, but if normals are
  needed the host must re-derive them after `read_state`.

## Numerical precision retrospective

- **Skalak** (largest absolute force magnitudes, ~10 μN): worst-fit 0.0127%
  relative error, 1.9e-4 μN absolute. Order-of-summation rounding from CSR
  aggregation, no driver-level math differences observed.
- **WLC** (smallest force magnitudes, ~1e-5 μN): worst-fit 1e-11 μN
  absolute. Marko-Siggia is monotone and well-conditioned; aggregation
  order doesn't matter at this scale.
- **DPD** (driven by Box-Muller through `sin`/`cos`/`log`/`sqrt`): worst-
  fit 5e-8 μN absolute. Non-trivial transcendentals; observed difference
  is at the f32 epsilon floor.
- **Verlet** (pure arithmetic): bit-identical (0.0 absolute).

f32 throughout was sufficient. No species or kernel needed Kahan
summation, log-space transforms, or any other escape hatch from the
Phase 11 plan.

## Test counts after 11.3.D

- 186 lib tests (unchanged).
- 138 integration tests (was 134; +4 from skalak/wlc/dpd/verlet GPU parity).
- Validation suite: 7/8 (unchanged).

## What's next

Phase 11.3.E is the integration layer: persistent buffers + one-pass
multi-substep dispatch + membrane damping + spectrin↔mesh coupling. After
that, Phase 11.4 (render-compute buffer sharing) and 11.5 (timestamp
profiling, deprecation of CoupledSolver) close out Phase 11.

Beyond Phase 11, the plan's roadmap continues with Phase 12 (external
fluid, single-cell capillary), Phase 13 (multi-cell microcirculation),
and **Phase 14 — storage lesion at physiological timescale**, which the
user has identified as the headline scientific differentiator.
