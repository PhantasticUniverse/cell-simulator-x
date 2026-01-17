# Cell Simulator X

GPU-accelerated human red blood cell simulation engine integrating mechanics, metabolism, and oxygen transport.

## Tech Stack

- **Core**: Rust
- **GPU Compute**: Metal (primary), CUDA (secondary), WebGPU (web)
- **UI**: TypeScript
- **Build**: Cargo, wasm-pack

## Commands

```bash
cargo build              # Build project
cargo test               # Run all tests
cargo bench              # Run benchmarks
cargo check              # Fast type checking
cargo clippy             # Lint
cargo fmt                # Format code
cargo run                # Run simulator
```

## Code Conventions

- **Units in names**: `velocity_um_per_sec`, `concentration_mM`, `pressure_mmHg`
- **Citations required**: All biological parameters must cite source (journal, year)
- **Parameter format**: JSON files in `data/parameters/` with citation metadata
- **Module structure**: See `src/` layout in PRD.md Section 11

## Key Modules

| Module | Purpose |
|--------|---------|
| `geometry` | RBC mesh, spectrin network, protein placement |
| `physics` | DPD solver, membrane mechanics, WLC model |
| `biochemistry` | Metabolic pathways, enzyme kinetics, hemoglobin |
| `integration` | Mechano-metabolic coupling |
| `compute` | Metal/CUDA/WebGPU backends |

## Architecture

See `PRD.md` for:
- Detailed architecture diagrams (Section 4)
- Parameter values and data sources (Section 5)
- Validation targets (Section 10)
- Phase roadmap (Section 7)

## Current Phase

**Phase 1: Foundation** - Project setup and geometry
- See PRD.md Section 7 for detailed phase checklist

## Critical Validation Targets

| Metric | Target |
|--------|--------|
| P50 oxygen | 26.8 ± 1 mmHg |
| ATP concentration | 1.5-2.5 mM |
| Membrane shear modulus | 5.5 μN/m |
| Simulation speed | >10 fps |
