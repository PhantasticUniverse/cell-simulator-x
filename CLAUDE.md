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

## Current Status

- Phase 1: Foundation ✅
- Phase 2: Mechanics ✅
- Phase 3: Metabolism ✅
- Phase 4: Oxygen Transport ✅
- Phase 5: Metabolism-Oxygen Integration ✅
- Phase 6a: Redox Metabolism ✅ (PPP, Glutathione, Piezo1)
- Phase 6b: Ion Homeostasis ✅ (Na+/K+-ATPase, PMCA)

## Development

**Always validate via CLI diagnostics before GUI testing.**

See `docs/DEVELOPMENT.md` for detailed implementation notes and `PRD.md` for architecture.
