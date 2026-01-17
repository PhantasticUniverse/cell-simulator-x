# Membrane Mechanics Domain Knowledge

> **Implementation Status**: ✅ Implemented in Phase 2 (`src/physics/`)

## DPD (Dissipative Particle Dynamics) Solver

**Implementation**: `src/physics/dpd.rs`

### Force Components
```
F_total = F_conservative + F_dissipative + F_random
```

**Conservative Force:**
```
F_C = a_ij * (1 - r_ij/r_c) * r̂_ij   for r < r_c
```

**Dissipative Force (drag):**
```
F_D = -γ * w_D(r) * (r̂_ij · v_ij) * r̂_ij
```

**Random Force (thermal):**
```
F_R = σ * w_R(r) * θ_ij * r̂_ij
```

**Fluctuation-Dissipation:**
```
σ² = 2 * γ * k_B * T
w_D(r) = [w_R(r)]²
```

### Implementation Parameters (DPDParameters)
| Parameter | Value | Description |
|-----------|-------|-------------|
| γ (gamma) | 4.5 | Dissipation coefficient |
| k_BT | 4.11e-9 μN·μm | Thermal energy at 37°C |
| r_c | 1.0 μm | Cutoff radius |
| a | 25.0 | Conservative force coefficient |

### Time Integration
- Velocity-Verlet with DPD correction (`src/physics/integrator.rs`)
- Default dt: 1e-5 s (10 μs)
- Temperature control via γ, σ balance

## Spectrin WLC (Worm-Like Chain) Model

**Implementation**: `src/physics/wlc.rs`

### Marko-Siggia Force-Extension
```
F = (k_B * T / L_p) * [1/(4*(1-x/L_c)²) - 1/4 + x/L_c]
```

### Implementation Parameters (WLCParameters)
| Parameter | Value | Description |
|-----------|-------|-------------|
| L_p | 20 nm (0.020 μm) | Persistence length |
| L_c | 200 nm (0.200 μm) | Contour length |
| L_rest | 75 nm (0.075 μm) | Rest length in situ |
| max_extension | 0.95 | Cap to prevent singularity |

Where:
- L_p = persistence length = 20 nm (updated from Rief et al. 1999)
- L_c = contour length = 200 nm (spectrin tetramer)
- x = end-to-end distance

### Network Topology
- ~33,000 spectrin tetramers
- Junctional complexes (actin + Band4.1)
- Hexagonal lattice (mostly)
- Defects: pentagons, heptagons

### Spring Constants
At rest length (70-80 nm):
- Effective k ≈ 10-20 pN/μm

## Skalak Strain Energy Function

**Implementation**: `src/physics/membrane.rs`

### Membrane Energy
```
W = (G_s/4) * (I₁² + 2*I₁ - 2*I₂) + (G_a/4) * I₂²
```

### Implementation Parameters (SkalakMaterial)
| Parameter | Value | Source |
|-----------|-------|--------|
| G_s (shear modulus) | 5.5 μN/m | Evans & Waugh 1977 |
| G_a (area modulus) | 450,000 μN/m (450 mN/m) | Evans et al. 1976 |
| κ (bending modulus) | 0.18 pN·μm | Evans 1983 |
| C₀ (spontaneous curvature) | 0 | Unstressed RBC |

Where:
- G_s = shear modulus = 5.5 μN/m (Dao et al. 2003)
- G_a = area dilation modulus ≈ 450 mN/m (near incompressibility)
- I₁ = α² + β² - 2 (strain invariant 1)
- I₂ = α²β² - 1 (strain invariant 2)
- α, β = principal stretch ratios

### Bending Energy (Helfrich)
Implementation uses discrete Laplacian approximation for bending forces.
```
W_bend = (κ/2) * (2H - C₀)² + κ_G * K
```

Where:
- κ = bending modulus = 0.18 pN·μm ≈ 44 k_BT (Evans 1983)
- H = mean curvature (computed via umbrella operator)
- C₀ = spontaneous curvature
- K = Gaussian curvature
- κ_G = Gaussian modulus ≈ -κ

## RBC Geometry

### Fung-Tong Parametric Surface
Biconcave disc in cylindrical coordinates:
```
z(r) = ±(1 - (r/R)²)^(1/2) * (C₀ + C₂*(r/R)² + C₄*(r/R)⁴)
```

Parameters:
- R = 3.91 μm (cell radius)
- C₀ = 0.81 μm
- C₂ = 7.83 μm
- C₄ = -4.39 μm

### Typical Values
| Property | Value | Source |
|----------|-------|--------|
| Volume | 90 fL | Canham 1970 |
| Surface area | 140 μm² | Canham 1970 |
| Excess area ratio | 1.4 | (SA > sphere of same V) |
| Thickness (rim) | 2.4 μm | |
| Thickness (center) | 1.0 μm | |

## Micropipette Aspiration Validation

### Force-Extension Response
At small deformations:
```
ΔL_p / R_p = (ΔP * R_p) / (2π * G_s)
```

Where:
- ΔL_p = aspiration length
- R_p = pipette radius
- ΔP = aspiration pressure

### Key Experimental Data (Dao et al. 2003)
- Linear response up to ~50% strain
- Shear modulus: 5.5 ± 1.1 μN/m
- Area modulus: 450-500 μN/m

## Viscosity

### Cytoplasmic Viscosity
- η_cytoplasm = 6 cP (Cokelet & Meiselman)
- ~5× water due to Hb concentration

### Membrane Viscosity
- η_membrane = 0.5-1.0 μN·s/m
- Important for dynamic deformation

## Key Citations

- Fedosov et al. (2010) Biophys J - DPD RBC model
- Dao et al. (2003) PNAS - Mechanical properties
- Evans (1983) Biophys J - Bending modulus
- Rief et al. (1999) Science - Spectrin mechanics
- Canham (1970) J Theor Biol - RBC geometry
