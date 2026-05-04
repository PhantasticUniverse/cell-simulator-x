# Phase 11.0 Notes — GPU Compute Scaffolding

Date: 2026-05-01

## Headline

Phase 11.0 is the GPU compute pipeline foundation. A `ComputeContext`
manages a headless wgpu device + queue, dispatches a sentinel `vec_add`
kernel, and asserts CPU = GPU bitwise on 1024 f32 elements. Apple M4 Max
detected, dispatch in 116 ms (cold start), result correct.

```
$ cargo run --release -- --diagnose-gpu
=== Cell Simulator X — GPU Compute Diagnostic (Phase 11.0) ===

Adapter:    Apple M4 Max (IntegratedGpu on Metal)
Setup:      21.43ms
Dispatch:   116.39ms (1024 f32 elements)
Result:     CPU = GPU bitwise on all 1024 elements ✓

=== Phase 11.0 sentinel: PASS ===
```

## Scope decisions

The Phase 11 plan opened with a dependency-upgrade sub-task: bump
wgpu 0.20 → 25.x, drag egui 0.28 → 0.34, winit 0.29 → 0.30. The day-one
pre-flight surfaced five distinct egui API breakages affecting the HUD
(`Context::run` → `run_ui`, `id_source` → `id_salt`, `Area::new` requires
`Id`, `Painter::rect` needs `StrokeKind`, `Rounding` → `CornerRadius`),
estimated at half a day to a day of mechanical migration plus revalidation.

The user pre-approved (during plan review) accepting HUD freeze if the
upgrade was painful. With wgpu 0.20 already shipping full compute support
(WGSL, atomics, workgroup memory, all backends), the practical reason to
upgrade is cleanup, not capability. We chose to **stay on wgpu 0.20** and
defer the dependency bump until a future phase forces it (likely
Phase 12 if newer wgpu features are needed for lattice-Boltzmann fluid).

This trade keeps Phase 11.0 focused on what's actually load-bearing —
proving the compute pipeline works end-to-end — and lets later sub-phases
inherit a working scaffold instead of a half-migrated one.

## Architecture

`ComputeContext` is intentionally separate from `RenderState`. Reasons:

- `RenderState` requires a window (wgpu 0.20's `create_surface` needs a
  surface-providing handle). Compute does not.
- CLI diagnostics (`--diagnose-gpu`, `--validate`, `--diagnose-multi-cell`)
  must run headless.
- Sharing buffers between render and compute requires sharing a `Device`,
  which we'll do in **sub-phase 11.4** when the simulator visualizes
  multi-cell GPU state. For 11.0 the contexts are independent.

Public surface:

```rust
pub struct ComputeContext {
    pub adapter: wgpu::Adapter,
    pub device: Arc<wgpu::Device>,
    pub queue:  Arc<wgpu::Queue>,
}

impl ComputeContext {
    pub async fn new_headless() -> Result<Self>;
    pub fn new_headless_blocking() -> Result<Self>;
    pub fn adapter_summary(&self) -> String;
}
```

`Arc<Device>` and `Arc<Queue>` so that downstream kernel helpers can be
written as free functions taking `&ComputeContext` without forcing a
particular borrow lifetime.

## Files added

- `src/compute/mod.rs` — `ComputeContext` (~110 lines).
- `src/compute/diagnostics.rs` — `vec_add_gpu` host wrapper, `run_diagnose_gpu`
  CLI entrypoint, in-module CPU=GPU bitwise test (~210 lines).
- `shaders/compute/vec_add.wgsl` — minimal `@compute @workgroup_size(64)`
  kernel.

## Files touched

- `src/lib.rs` — `pub mod compute;`.
- `src/main.rs` — `--diagnose-gpu` CLI flag.

## Test counts

- 186 lib unit tests (was 185; +1 from `compute::diagnostics::tests::vec_add_cpu_gpu_bitwise`).
- 130 integration tests (unchanged across mechanics, metabolism, oxygen,
  redox, ion, integration, disease, coupling, world, validation).
- Validation suite: 7/8 pass (unchanged from Phase 10/10.5).

## Exit criteria status

From the Phase 11 plan, sub-phase 11.0 exit criteria:

- ✅ Existing render path still works (no changes to render code).
- ✅ Validation suite still 7/8 (no metabolism/mechanics changes).
- ✅ `cargo run --release -- --diagnose-gpu` reports adapter info,
   dispatches `vec_add` on 1024 elements, asserts CPU=GPU bitwise,
   prints timestamped duration.

## What's next (sub-phase 11.1)

Phase 11.1 lands the SoA refactor on CPU (deferred Phase 10.5 items):

- `PhysicsState` from `Vec<Vec3>` to parallel `Vec<f32>` arrays.
- `MetabolitePoolBatch` species-major variant alongside the existing
  `MetabolitePool`.
- Double-buffered (read/write) per-cell state in `World`.

Exit: validation suite still 7/8, world tests still 100%, no GPU code
touched, `--diagnose-multi-cell 10` perf within 5% of pre-refactor
baseline.

## Decisions deferred

- **Dependency upgrade** (wgpu 0.20 → 29.x, egui/winit alongside).
  Triggered by future phases needing wgpu features not in 0.20, or by
  egui ecosystem decisions outside Phase 11. Pre-flight researched
  versions: target triple is `wgpu 29.0`, `egui 0.34`, `egui-wgpu 0.34`,
  `egui-winit 0.34`, `winit 0.30`.
- **Render-compute device sharing**. Phase 11.4 deliverable; depends on
  Phase 11.2 (biochem on GPU) and Phase 11.3 (mechanics on GPU) landing
  first so we know what state needs to be shared.
- **Timestamp queries / wgpu-profiler**. Phase 11.5 deliverable.
