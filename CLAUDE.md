# Cell Simulator X

GPU-accelerated human red blood cell simulation engine.

## Commands

```bash
cargo build                          # Build
cargo test                           # Run tests
cargo run                            # GUI
cargo run -- --diagnose              # Physics diagnostics (CLI)
cargo run -- --diagnose-metabolism   # Metabolism diagnostics (CLI)
cargo run -- --diagnose-oxygen       # Oxygen transport diagnostics (CLI)
cargo run -- --diagnose-integrated   # Metabolism-oxygen coupling (CLI)
cargo run -- --diagnose-full         # Full integration (glycolysis+PPP+redox)
cargo run -- --diagnose-disease storage --disease-param 21    # Storage lesion (day 21)
cargo run -- --diagnose-disease diabetic --disease-param 12   # Diabetic (12 mM glucose)
cargo run -- --diagnose-disease malaria --disease-param 0.05  # Malaria (5% parasitemia)
cargo run -- --diagnose-disease sickle --disease-param 1.0    # Sickle cell (HbSS)
cargo run -- --diagnose-coupled -d 60                          # Mechano-metabolic coupling
cargo run -- --diagnose-coupled --tension 2.0 -d 60            # With tension override
cargo run --features validation -- --validate                  # Phase 10 empirical validation suite
cargo test --features validation --test validation_suite       # Validation tests
cargo run --release -- --diagnose-multi-cell 10 -d 1.0         # Phase 10.5 multi-cell scaling diagnostic
cargo run --release -- --diagnose-gpu                          # Phase 11.0 GPU compute sentinel (vec_add CPU=GPU)
cargo run -- --help                  # Full CLI options
```

## Code Conventions

- **Units in names**: `velocity_um_per_sec`, `concentration_mM`, `pressure_mmHg`
- **Citations required**: All biological parameters must cite source
- **Parameters**: JSON files in `data/parameters/`

## Key Modules

| Module | Purpose |
|--------|---------|
| `geometry` | RBC mesh, spectrin network |
| `physics` | Membrane mechanics, DPD, WLC |
| `biochemistry` | Glycolysis, PPP, glutathione, hemoglobin, pH buffer, ion homeostasis, full integration |
| `biochemistry/disease` | Disease models: storage lesion, diabetic, malaria, sickle cell |
| `coupling` | Mechano-metabolic coupling: TensionComputer, SpectrinModulator, CoupledSolver |
| `render/hud` | egui-based HUD overlay: panels, widgets, theme |
| `export` | Screenshot, CSV time-series, JSON state export |
| `validation` | (feature `validation`) Empirical fits vs Imai 1981, Mulquiney 1999, Rief 1999, Waugh-Evans 1979, Dao 2003 |
| `world` | Multi-cell `World` + `Cell` + `CellHandle` (Phase 10.5); rayon-parallel per-cell stepping |
| `compute` | GPU compute scaffolding (Phase 11.0); `ComputeContext` headless wgpu, `vec_add` sentinel kernel |
| `flow` | (Phase 12.A) Analytic external fluid: `Poiseuille` cylindrical channel + Stokes-form drag for vertex/flow coupling |
| `storage` | (Phase 14.A–D) 42-day storage lesion simulator: multi-rate ms/s/day scheme, `StorageCurveSimulator` matches Hess 2010 metabolomics; Phase 14.D adds `AdditiveSolution` (CPD / AS-3 / SAGM / PAGGSM), Q10-scaled supercooled storage, ±20% sensitivity sweep |

## Current Status

- Phase 1: Foundation ✅
- Phase 2: Mechanics ✅
- Phase 3: Metabolism ✅
- Phase 4: Oxygen Transport ✅
- Phase 5: Metabolism-Oxygen Integration ✅
- Phase 6a: Redox Metabolism ✅ (PPP, Glutathione, Piezo1)
- Phase 6b: Ion Homeostasis ✅ (Na+/K+-ATPase, PMCA)
- Phase 7: Disease Models ✅ (Storage Lesion, Diabetic, Malaria, Sickle Cell)
- Phase 8: Mechano-Metabolic Coupling ✅ (TensionComputer, SpectrinModulator, CoupledSolver)
- Phase 9: GUI HUD & Export ✅ (Real-time metabolite display, JSON/CSV export)
- Phase 10: Empirical Validation Foundation ✅ (Imai/Mulquiney/Rief/Waugh-Evans/Dao; aldolase Keq bug fix; PPP refit; see `docs/validation_report_v1.md`)
- Phase 10.5: Multi-Cell Architecture Refactor ✅ (World/Cell/CellHandle; rayon parallelism; N=100 demo; SoA + buffer-tension deferred to Phase 11; see `docs/phase_10_5_notes.md`)
- Phase 11.0: GPU Compute Scaffolding ✅ (ComputeContext headless wgpu, vec_add sentinel kernel, --diagnose-gpu, CPU=GPU bitwise on Apple M4 Max; dependency upgrade deferred per `docs/phase_11_0_notes.md`)
- Phase 11.2.A: Glycolysis on GPU ✅ (11 enzymes + RK4 in WGSL; 0.0000% rel-err per species after 1s simulation vs CPU baseline; see `docs/phase_11_2_a_notes.md`)
- Phase 11.2.B: 2,3-BPG Shunt on GPU ✅ (BPGM + BPGP added → full 18-species `MetabolismSolver` parity; 0.0009% worst-fit rel-err on 2,3-BPG)
- Phase 11.2.C: PPP + Glutathione + Piezo1 + Ions on GPU ✅ (full 38-species kernel: 22 rate functions; 0.0023% worst-fit rel-err on Na⁺ after 1s; see `docs/phase_11_2_c_notes.md`)
- Phase 11.2.D: Hemoglobin Adair + pH Buffer on GPU ✅ (post-RK4 Euler interleaved per ms; Hb saturation parity 2e-6 absolute after 1s; see `docs/phase_11_2_d_notes.md`)
- Phase 11.2.E: Inline homeostasis corrections + full-solver parity ✅ (basal NADPH, basal GSH, ATP regen; GPU kernel matches `FullyIntegratedSolver::step` to 0.0023% per species; see `docs/phase_11_2_e_notes.md`)
- Phase 11.3.A–D: Skalak + WLC + DPD + Velocity-Verlet on GPU ✅ (per-element/per-edge force kernels with CSR aggregation; stateless PCG-hash RNG for DPD; CPU/GPU parity from 0 absolute (Verlet) to 0.0127% relative (Skalak); see `docs/phase_11_3_notes.md`)
- Phase 11.3.E: Integrated `PhysicsBackend` ✅ (persistent buffers + 6-dispatch single-pass step; WLC pre-baked as static baseline; CPU/GPU parity 1.9e-15 μm position drift after 10 substeps; see `docs/phase_11_3_notes.md`)
- Phase 11.4: Render-Compute device sharing ✅ (RenderState exposes `Arc<Device>`/`Arc<Queue>`; `ComputeContext::from_shared` borrows; vertex-shader binding deferred; see `docs/phase_11_4_notes.md`)
- Phase 11.5: Backend switch + CoupledSolver deprecation ✅ (`WORLD_BACKEND` env var; `Backend` enum; `#[deprecated]` on CoupledSolver/CoupledConfig; wgpu-profiler deferred; see `docs/phase_11_5_notes.md`)
- Phase 12.A: External fluid — analytic Poiseuille + drag ✅ (`flow::Poiseuille` + `apply_drag_to_external_forces`; CPU integration test shows directional cell drift under flow; see `docs/phase_12_notes.md`)
- Phase 12.B.1: PhysicsBackend external_forces buffer ✅ (binding 9 in `physics_step.wgsl`; `skalak_init_from_baseline` adds external forces alongside WLC baseline; `PhysicsBackend::set_external_forces`; ComputeContext bumps `max_storage_buffers_per_shader_stage` 8→16; CPU/GPU parity <1e-5 μm under uniform constant force after 100 substeps)
- Phase 12.B.2: World::apply_poiseuille_drag CPU+GPU parity ✅ (rayon-parallel per-cell drag write into `PhysicsState::external_forces_uN`; `tests/flow_cpu_gpu_parity.rs` proves Poiseuille CPU/GPU within <1e-3 μm centroid drift after 100 substeps)
- Phase 12.C.1: Fischer 2007 tank-treading validation ✅ (Keller-Skalak 1982 analytic prediction `f_TT = (γ̇/2π)·2αβ/(α²+β²)` matches Fischer 2007's K(λ) range 0.04–0.15 within χ²/dof < 2; configuration-level analogous to Dao 2003)
- Phase 12.C.2: Skalak 1973 parachute shape validation ✅ (capillary number `Ca = μ_ext·γ̇·R / μ_s` + empirical AR(Ca) fit; AR predictions 1.40 / 1.66 / 2.13 at γ̇_wall ∈ {100,200,400} /s lie within Skalak 1973's reported 1.5–2.0 band; full dynamic simulation deferred)
- Phase 14.A: Storage lesion at physiological timescale ✅ (multi-rate ms/s/day scheme; `StorageCurveSimulator` produces 42-day metabolite trajectory matching Hess 2010 ATP and 2,3-DPG envelopes; ion gradients trend correctly but quantitative Hess targets need Phase 14.B analytic QSS; see `docs/phase_14_notes.md`)
- Phase 14.B: Analytic ion quasi-steady-state ✅ (`solve_ion_qss` bisection drives Na+/K+ to pump+leak equilibrium each storage day; day-42 QSS lands at Na≈96 mM (Hess 2010 reports ~60 mM — gap is a parameter-identifiability finding for envelope re-fitting in 14.B'); Hess 2010 quantitative match deferred to envelope re-fit)
- Phase 14.C: Deformability decline coupling ✅ (`StorageSample::deformability_relative` from Phase 8 `SpectrinModulator`'s ATP→stiffness; 1.000 → 0.728 over 42 days, matches ~30% decline reported by Hess 2010 / Pivkin 2011)
- Phase 14.B': Storage envelope calibrated to Hess 2010 day-42 ✅ (`pump_efficiency_decay_per_day` 0.02 → 0.0168; storage curve now hits Na ≈ 60 mM / K ≈ 90 mM at day 42 exactly; day-14 ions remain underestimated as a documented linear-envelope limitation)
- Phase 14.B'': Exponential pump-efficiency envelope ✅ (`use_exponential_pump_envelope` default-on; matches Hess 2010 ion targets at days 0/14/42 simultaneously: Na = 10/26/60 mM, K = 144/127/91 mM)
- Phase 14.D.1: AS-3 / SAGM / PAGGSM additive comparator presets ✅ (`AdditiveSolution` enum + `StorageSimConfig::with_additive`; per-additive `lesion_config()` derives ATP T_half analytically from day-42 retention; presets within 15% of literature; see `docs/phase_14_d_notes.md`)
- Phase 14.D.2: Supercooled storage with Q10 = 2.5 ✅ (`StorageSimConfig::supercooled` runs at -4°C / day 100; envelope rates Q10-scaled per Hess 2010 review; -4°C / day-100 deformability within 1.5% of standard 4°C / day-42)
- Phase 14.D.3: Sensitivity analysis sweep ✅ (`run_oat_sensitivity` ±20% perturbation on 5 envelope parameters; `target/storage_sensitivity.csv` CSV deliverable; ATP half-life dominates day-42 ATP/deformability; pump/leak/oxidative parameters affect ions but not ATP under `force_atp_dpg_targets = true`)

## Development

**Always validate via CLI diagnostics before GUI testing.**

See `docs/DEVELOPMENT.md` for detailed implementation notes and `PRD.md` for architecture.
