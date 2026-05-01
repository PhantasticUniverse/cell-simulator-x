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
- Phase 10.5: Multi-Cell Architecture Refactor 🚧 (in progress — World/Cell/CellHandle; rayon parallelism; N=100 demo)

## Development

**Always validate via CLI diagnostics before GUI testing.**

See `docs/DEVELOPMENT.md` for detailed implementation notes and `PRD.md` for architecture.
