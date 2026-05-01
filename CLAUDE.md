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

## Development

**Always validate via CLI diagnostics before GUI testing.**

See `docs/DEVELOPMENT.md` for detailed implementation notes and `PRD.md` for architecture.
