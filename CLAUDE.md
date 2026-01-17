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

**Phase 1: Foundation** - ✅ COMPLETE

### Phase 1 Features
- Fung-Tong parametric RBC geometry (biconcave disc)
- Surface mesh generation (~10K vertices)
- Spectrin network graph (~33K tetramers, hexagonal lattice)
- WebGPU/Metal rendering pipeline with Phong shading
- Orbit camera with mouse controls
- Parameter system with JSON configs and citations
- State structures for biochemistry, physics, environment

---

**Phase 2: Mechanics** - ✅ COMPLETE

### Phase 2 Features
- **WLC Spectrin Elasticity**: Marko-Siggia force-extension model
  - Persistence length: 20 nm
  - Contour length: 200 nm
  - Rest length: 75 nm
- **Skalak Membrane Mechanics**: Strain energy model
  - Shear modulus: 5.5 μN/m (Evans & Waugh 1977)
  - Area modulus: 450 mN/m
  - Bending modulus: 0.18 pN·μm (Evans 1983)
- **DPD Solver**: Dissipative Particle Dynamics
  - Conservative forces (soft repulsion)
  - Dissipative forces (friction)
  - Random forces (thermal fluctuations)
  - Fluctuation-dissipation theorem: σ² = 2γk_BT
- **Velocity-Verlet Integrator**: Symplectic time integration
  - Adaptive timestep control
  - Velocity/displacement clamping for stability
- **Dynamic Mesh Rendering**: GPU buffer updates for real-time deformation

### New Physics Module Structure
```
src/physics/
├── mod.rs        # PhysicsSolver combining all forces
├── wlc.rs        # WLC spectrin elasticity
├── membrane.rs   # Skalak strain energy model
├── dpd.rs        # DPD fluid dynamics
└── integrator.rs # Velocity-Verlet integration
```

### Controls
- **Mouse drag**: Orbit camera
- **S key**: Toggle spectrin network overlay
- **R key**: Reset camera
- **P key**: Toggle physics simulation (starts paused)
- **F key**: Apply force to center vertex (micropipette simulation)
- **+/- keys**: Adjust physics substeps per frame
- **Escape**: Exit

### Next Phase
**Phase 3: Core Metabolism** - ODE solver and glycolysis implementation
- See PRD.md Section 7 for detailed phase checklist

## Critical Validation Targets

| Metric | Target |
|--------|--------|
| P50 oxygen | 26.8 ± 1 mmHg |
| ATP concentration | 1.5-2.5 mM |
| Membrane shear modulus | 5.5 μN/m |
| Simulation speed | >10 fps |
