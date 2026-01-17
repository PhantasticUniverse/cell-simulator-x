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
| `biochemistry` | Glycolysis, enzyme kinetics, hemoglobin |

## Current Status

- Phase 1: Foundation ✅
- Phase 2: Mechanics ✅
- Phase 3: Metabolism ✅
- Phase 4: Oxygen Transport ✅
- Phase 5: Integration (next)

## Development

**Always validate via CLI diagnostics before GUI testing.**

See `docs/DEVELOPMENT.md` for detailed implementation notes and `PRD.md` for architecture.
