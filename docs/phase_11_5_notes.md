# Phase 11.5 Notes — Backend Switch + CoupledSolver Deprecation

Date: 2026-05-02

## Headline

Three pieces from the Phase 11.5 plan landed:

1. **`WORLD_BACKEND=cpu|gpu` env var infrastructure** — `world::backend`
   module exposes a `Backend` enum and `env_default()` reader. Tests and
   future call sites can flip between CPU and GPU paths without
   recompiling.
2. **`CoupledSolver` / `CoupledConfig` `#[deprecated]`** — formal warning
   emitted at every consumer use site, scheduled for removal in Phase 12.
3. **wgpu-profiler integration** — **deferred** (see "What's deferred"
   below). Adding a new dependency + instrumenting every dispatch on
   wgpu 0.20 is mechanical but substantial enough to warrant its own
   commit; this Phase 11.5 commit ships the cheaper deliverables first.

## What landed

### `world::backend` module (`src/world/backend.rs`)

```rust
pub enum Backend { Cpu, Gpu }
pub const ENV_VAR: &str = "WORLD_BACKEND";
pub fn env_default() -> Backend;
```

Default is `Backend::Cpu`. Unknown values silently fall back to default
(no panic on typos). Three new lib tests exercise parsing.

Re-exported at `world::Backend` and `world::backend_from_env`.

### CoupledSolver deprecation

`#[deprecated(since = "1.1.0", note = "...")]` attribute on
`CoupledSolver` and `CoupledConfig`. The note points to `World` +
`compute::PhysicsBackend` as the migration target. Removal is scheduled
for Phase 12.

To keep build output clean, `#![allow(deprecated)]` is set at module
scope on:

- `src/coupling/coupled_solver.rs` — the file that defines them.
- `src/coupling/mod.rs` — the re-export.
- `src/lib.rs` (per-export) — the lib-root re-export.
- `src/world/cell.rs` and `src/world/mod.rs` — `World`/`Cell` still take
  `CoupledConfig` (parameter parity until Phase 12).
- `src/main.rs` — `--diagnose-coupled` and `--diagnose-multi-cell` CLI
  modes still use `CoupledSolver`.
- `tests/coupling_tests.rs` and `tests/world_tests.rs` — guarded test
  paths for the legacy code.

Net effect: `cargo build --tests` emits **0** deprecation warnings, but
any new external user of `CoupledSolver` or `CoupledConfig` will see the
warning at their use site (which is the point of the deprecation).

## What's deferred

### wgpu-profiler integration

The plan called for wgpu-profiler ≥ 0.24 with `Features::TIMESTAMP_QUERY`,
exposed in the HUD as per-dispatch timing
(`physics: 12 ms · biochem: 4 ms · copy: 0.3 ms · total: 17.4 ms/step`).
This is genuinely useful for the multi-cell scaling phase (Phase 13)
where measuring kernel-vs-dispatch overhead matters.

Why deferred:

- wgpu 0.20 pairs with wgpu-profiler 0.16 (not 0.24). The version graph
  stays consistent only if we keep things simple.
- Instrumenting every dispatch in `PhysicsBackend.step()` and
  `run_full_biochem_batch_with_hb()` is one-by-one work touching ~10
  call sites.
- Surfacing timing in the HUD requires wiring through the egui panels,
  which interacts with the HUD freeze.

The right time to land wgpu-profiler is when the speedup gate becomes
measurable — i.e., once Phase 13 (multi-cell) is in flight. Until then,
overall wall-clock measurement (`std::time::Instant`) is sufficient.

### Vertex-shader binding for direct render-compute buffer sharing

This was the second piece of Phase 11.4 (the API piece landed; the
shader piece deferred). Same comment as Phase 11.4 notes — it's
mechanical but better landed once `World` + `PhysicsBackend` integrate.

### `World::step` migration to `PhysicsBackend`

Currently `World::step` calls `Cell::step` which uses CPU
`PhysicsSolver`. Wiring `Backend::Gpu` to use `PhysicsBackend` is the
forcing function for that migration; it gates Phase 13 multi-cell
scaling. Touched lightly here (env-var reads it but no dispatch path
checks it yet).

## Test counts after 11.5

- 189 lib tests (was 186; +3 from `world::backend::tests`).
- 140 integration tests (unchanged).
- Validation suite: 7/8 (unchanged).

## Phase 11 retrospective

With 11.0 → 11.5 landed, every CPU subsystem in this simulator now has a
working GPU equivalent that's parity-tested against its CPU baseline:

| Component | CPU | GPU | Parity |
|-----------|-----|-----|--------|
| Glycolysis (11 enzymes + 2,3-BPG shunt) | `MetabolismSolver` | `run_glycolysis_batch` | 0.0009% rel-err |
| 38-species biochem | `FullyIntegratedSolver` | `run_full_biochem_batch_with_hb` | 0.0023% rel-err |
| Skalak forces | `SkalakSolver` | `run_skalak_forces` | 0.0127% rel-err |
| WLC forces | `WLCSolver` | `run_wlc_forces` | 1e-11 μN abs |
| DPD forces | `DPDSolver` (PCG path) | `run_dpd_membrane_forces` | 5e-8 μN abs |
| Velocity-Verlet | `VelocityVerlet` | `run_verlet_step` | bit-identical |
| Integrated step | `PhysicsSolver::step` | `PhysicsBackend::step` | 1.9e-15 μm pos drift / 10 substeps |

Phase 11 is now **closed** for the algorithmic deliverables. What
remains is integration glue (wgpu-profiler, vertex-shader binding,
`World` migration) — none of which is gating for the headline
scientific result.

## What's next

Phase 12 (external fluid in single-cell capillary) is the next entry on
the strategic roadmap. After that, Phase 13 (multi-cell microcirculation)
and **Phase 14 — storage lesion at 42-day timescale**, the headline
differentiator.

The biochemistry GPU kernel (Phase 11.2) and the integrated
`PhysicsBackend` (Phase 11.3.E) are the foundation Phase 14 will build
on — μs mechanics inner loop, s metabolism middle loop, hours-days
disease outer loop using analytic envelopes for the slow processes.
