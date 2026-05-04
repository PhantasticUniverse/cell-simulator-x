# Phase 11.4 Notes — Render-Compute Device Sharing

Date: 2026-05-02

## Headline

`RenderState` and `ComputeContext` can now share a single `Arc<Device>` +
`Arc<Queue>`. The `ComputeContext::from_shared(...)` constructor and the
`RenderState::shared_device()` / `shared_queue()` accessors together let
the render pipeline and the compute pipeline allocate buffers against the
same physical device — a precondition for direct render-compute buffer
sharing.

This is the **API-level** piece of Phase 11.4. The next piece (the
vertex-shader binding to a compute-managed position buffer) is mechanical
but touches the existing `cell.wgsl` render shader; we land that as a
follow-on so this commit stays focused on the device-sharing primitives.

## What landed

### `RenderState` (src/render/pipeline.rs)

- `device` and `queue` were `wgpu::Device` / `wgpu::Queue` (owned). Now
  they are `Arc<wgpu::Device>` / `Arc<wgpu::Queue>`. Internal call sites
  (`&self.device`, `&self.queue`) are unchanged because `&Arc<T>`
  deref-coerces to `&T`.
- Two new public accessors:
  - `RenderState::shared_device() -> Arc<wgpu::Device>`
  - `RenderState::shared_queue() -> Arc<wgpu::Queue>`

### `ComputeContext` (src/compute/mod.rs)

- `adapter` was `wgpu::Adapter`; now `Option<wgpu::Adapter>`. Headless
  contexts (`new_headless`) still own one; shared contexts don't.
- New constructor: `ComputeContext::from_shared(device, queue) -> Self`.
- `adapter_summary()` handles both cases.

### Test

`tests/compute_context_shared.rs::from_shared_dispatches_correctly` —
clones `Arc<Device>`/`Arc<Queue>` from a headless context, builds a
shared context, dispatches `vec_add_gpu`, asserts correctness.

## Architectural notes

- The render and compute paths bind `Arc<Device>` and `Arc<Queue>`
  directly into command encoders/passes; `Arc::clone` is the only cost,
  and dispatching against the same device guarantees that buffers
  created in one context are usable in the other.
- wgpu **cannot** share buffers across separate `Device` instances even
  on the same `Adapter` — that's why this refactor is required.
- The HUD freeze accepted in the Phase 11 plan stays: egui 0.28 +
  egui-wgpu 0.28 still bind to whichever `Device` the render pipeline
  passes them, which is now an `Arc<Device>` (deref-coerces to `&Device`
  for the egui APIs).

## What's still deferred

- **Vertex shader binding to compute-managed position buffer.**
  Currently the render pipeline binds a vertex buffer that's filled from
  CPU-side `Mesh::vertices` via `update_mesh()`. To eliminate the CPU
  round-trip during the substep hot path, the vertex shader would read
  positions from a `storage` buffer owned by `PhysicsBackend`. That
  needs:
  - A new vertex shader entry point (or a shared shader) that uses
    `@vertex fn vs_main(@builtin(vertex_index) vid: u32)` and reads
    positions from a storage binding.
  - `RenderState::bind_compute_buffer(buf: &wgpu::Buffer)` to wire the
    backend's positions buffer into the render bind group.
  - Synchronization: render reads should occur after compute writes
    finish (wgpu handles this within a single command encoder; if they
    cross encoders we need a queue.submit sync point).
  This is mechanical and adds no algorithmic risk, but it requires
  refactoring the existing vertex layout. It will land in a follow-on
  patch (Phase 11.4.B) once we wire `PhysicsBackend` into the `World`
  step path.

- **`World` integration.** `World::step()` currently calls
  `Cell::step()`, which uses CPU `PhysicsSolver`. Migrating it to
  `PhysicsBackend` is a small surgery but interacts with how cells own
  their state. We defer that to the multi-cell scaling phase.

## Test counts after 11.4

- 186 lib tests (unchanged).
- 140 integration tests (was 139; +1 from `compute_context_shared`).
