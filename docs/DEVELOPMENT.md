# Development Guide

## Project Status

### Phase 1: Foundation ✅
- Fung-Tong parametric RBC geometry (biconcave disc)
- Surface mesh generation (~10K vertices)
- Spectrin network graph (~33K tetramers, hexagonal lattice)
- WebGPU/Metal rendering pipeline with Phong shading
- Orbit camera with mouse controls
- Parameter system with JSON configs and citations

### Phase 2: Mechanics ✅
- **WLC Spectrin Elasticity**: Marko-Siggia force-extension (Lp=20nm, Lc=200nm)
- **Skalak Membrane**: Shear 5.5 μN/m, Area 450 mN/m, Bending 0.18 pN·μm
- **DPD Solver**: Conservative, dissipative, random forces (σ²=2γk_BT)
- **Velocity-Verlet**: Symplectic integration with adaptive timestep

### Phase 3: Core Metabolism ✅
- **RK4 ODE Integrator**: 4th-order Runge-Kutta, ~155K steps/sec
- **Enzyme Framework**: Michaelis-Menten, Hill, reversible bi-bi mechanisms
- **Glycolysis**: 11 enzymes (HK→LDH) with proper kinetics
- **Rapoport-Luebering Shunt**: BPGM/BPGP for 2,3-DPG regulation
- **Validated**: ATP 1.5-2.5 mM, 2,3-DPG 4.5-5.5 mM, P50 ~27 mmHg

### Phase 4: Oxygen Transport ✅
- **Adair 4-Site Model**: Sequential O2 binding with cumulative constants
- **Allosteric Effects**: Bohr (pH), 2,3-DPG, temperature (van't Hoff), CO2
- **Dynamic Kinetics**: O2 uptake/release with on/off rates
- **Validated**: P50 26.8±1 mmHg, Hill n 2.7±0.1, Bohr -0.48±0.05

### Phase 5: (Next)
- Integration of oxygen transport with metabolism
- Full cell simulation with coupled systems

## Module Structure

```
src/
├── geometry/       # RBC mesh, spectrin network
├── physics/        # DPD, membrane mechanics, WLC, integrator
├── biochemistry/   # Metabolism, enzyme kinetics
├── render/         # WebGPU/Metal rendering
├── config/         # Parameters, JSON loading
└── state/          # Cell state management
```

## CLI Diagnostics

### Physics
```bash
cargo run -- --diagnose              # Default: 1000 steps, 5 μN
cargo run -- --diagnose -n 5000 -f 50.0
```
- Micropipette simulation (downward force on center vertex)
- Reports displacement, velocity, energy

### Metabolism
```bash
cargo run -- --diagnose-metabolism           # Default: 60s
cargo run -- --diagnose-metabolism -d 30.0   # Custom duration
cargo run -- --diagnose-metabolism --glucose-step 5.0  # Perturbation
```
- Glycolysis + shunt simulation
- Reports ATP, 2,3-DPG, glucose, lactate

### Oxygen Transport
```bash
cargo run -- --diagnose-oxygen               # Default: standard conditions
cargo run -- --diagnose-oxygen --ph 7.2      # Test Bohr effect
cargo run -- --diagnose-oxygen --dpg 7.0     # Test 2,3-DPG effect
cargo run -- --diagnose-oxygen --temp 40     # Test temperature effect
```
- Adair equation with allosteric effects
- Generates OEC (oxygen equilibrium curve)
- Validates P50, Hill coefficient, Bohr coefficient

## GUI Controls

| Key | Action |
|-----|--------|
| Mouse drag | Orbit camera |
| S | Toggle spectrin overlay |
| R | Reset camera |
| P | Toggle physics |
| F | Apply force |
| +/- | Adjust substeps |
| Esc | Exit |

## Validation Targets

| Metric | Target | Source |
|--------|--------|--------|
| P50 oxygen | 26.8 ± 1 mmHg | Imai 1982 |
| Hill coefficient | 2.7 ± 0.1 | Imai 1982 |
| Bohr coefficient | -0.48 ± 0.05 | Imai 1982 |
| 2,3-DPG sensitivity | ~2.4 mmHg/mM | Benesch 1969 |
| ATP | 1.5-2.5 mM | Beutler 1984 |
| 2,3-DPG | 4.5-5.5 mM | Benesch 1969 |
| Shear modulus | 5.5 μN/m | Evans 1977 |

## References

See parameter JSON files in `data/parameters/` for full citations.
