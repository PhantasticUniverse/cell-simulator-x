# Phase 10.5 Notes

Date: 2026-05-01
Commit: `04d1912` ("Phase 10.5 milestone 1: introduce multi-cell World API")

## What landed

A new top-level `world` module exposing `World`, `Cell`, and `CellHandle`.
`World` owns N independent cells and steps them in parallel via rayon's
`par_iter_mut`. Each `Cell` bundles mesh, spectrin network, physics
state, physics solver, biochemistry solver, metabolite pool, and tension
history into one self-contained object.

`World::step` is the canonical multi-cell API for Phase 11+. The
existing single-cell `CoupledSolver` continues to work unchanged so all
existing callers, tests, and CLI commands keep functioning.

## Exit criteria

The Phase 10.5 plan stipulated three exit criteria. Status:

| Criterion | Status | Evidence |
|---|---|---|
| `World::step` runs N=1, N=10, N=100 cells on CPU | ✅ | `tests/world_tests.rs` covers N=1, 10, 100. `--diagnose-multi-cell N` runs at any N. |
| All Phase 10 validation tests pass at N=1 | ✅ | Validation suite still 7/8; world tests verify ATP parity across cells. |
| DPD neighbor-list correctness verified | N/A | The current DPD on membrane vertices is O(N) per-vertex (`compute_membrane_forces`), not the O(N²) `compute_forces` path. The O(N²) path is unused by membrane mechanics, so the spatial-hash refactor is unblocked but not yet needed. Deferred to Phase 12 (external fluid). |

## Performance

`--diagnose-multi-cell 10 -d 1.0` on Apple Silicon, release build:

```
N=1  baseline:  60.64s wall-clock for 1.0s simulated
N=10:          126.61s wall-clock for 1.0s simulated × 10 cells
Scaling:       2.09× baseline (vs 10× without parallelism)
Per-cell:      0.21× baseline cost per cell
```

This is ~5× speedup from rayon parallelism on a 10-thread CPU,
demonstrating the embarrassingly-parallel structure works as intended.
Per-cell ATP across the 10 cells: 2.018 mM mean, with negligible spread
(only DPD random forces vary).

## Scope decisions made during execution

The Phase 10.5 plan from
`~/.claude/plans/it-s-clear-what-this-functional-galaxy.md` listed five
architectural changes:

1. **`World` + `CellHandle` API** — done.
2. **SoA per-cell physics state** — *deferred to Phase 11*.
3. **Species-major SoA biochemistry state** — *deferred to Phase 11*.
4. **Double-buffered state** — *deferred to Phase 11*.
5. **Eliminate per-substep CPU readback in `TensionComputer`** —
   *deferred to Phase 11*.

The deferrals are intentional. Items 2–5 only deliver CPU-side value
when there is a GPU consumer. Adding them now would be premature
abstraction with no measurable improvement on the CPU baseline. They
are natural deliverables to land *with* the wgpu compute kernels in
Phase 11, where the storage-buffer layout and dispatch chaining are
the actual reasons for the refactor.

The `World` API was the irreducible Phase 10.5 deliverable: it
restructures the *control flow* (how cells are addressed and stepped)
in a way that is impossible to backfill once GPU kernels exist. The
*data layout* changes (SoA, double-buffered) and the *coupling
elimination* (buffer-based tension) can be done in lockstep with the
kernels.

## What's surfaced for Phase 11

- Replace per-cell `Vec<Vec3>` (in `PhysicsState`) with parallel
  `Vec<f32>` arrays for x/y/z components. wgpu storage buffers want
  this layout.
- Replace `MetabolitePool::concentrations_mM: Vec<f64>` with
  species-major SoA for multi-cell GPU dispatch (one buffer per
  metabolite of length `n_cells`, not one buffer per cell of length
  `n_metabolites`).
- Add double-buffered state (read/write) so single-dispatch shaders
  don't alias.
- Migrate `TensionComputer` to write into a per-cell tension buffer
  that biochemistry reads via storage-buffer binding rather than via
  host code.
- Add neighbor-list / spatial-hash to the unused `DPDSolver::compute_forces`
  path before any inter-cell DPD interactions land in Phase 13.

## Files

- `src/world/mod.rs` — `World`, `CellHandle`, parallelism.
- `src/world/cell.rs` — `Cell` (per-cell state + step).
- `src/main.rs` — `--diagnose-multi-cell N` flag.
- `tests/world_tests.rs` — 7 integration tests (N=1, N=10, N=100, ATP
  parity, handle stability, reset, tension override).
- `src/world/mod.rs` (in-file) — 5 unit tests (empty world, add_cell,
  step at N=1 and N=10, handle independence).

## Test counts

- 185 lib unit tests (was 180; +5 from `world::tests`).
- 7 integration tests in `tests/world_tests.rs` (new).
- Validation suite: 7/8 passing (no change from Phase 10).
