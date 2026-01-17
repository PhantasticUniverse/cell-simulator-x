# Membrane Mechanics Domain Knowledge

## DPD (Dissipative Particle Dynamics) Solver

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

### Time Integration
- Velocity-Verlet with DPD correction
- Typical dt: 0.001-0.01 τ (DPD time units)
- Temperature control via γ, σ balance

## Spectrin WLC (Worm-Like Chain) Model

### Marko-Siggia Force-Extension
```
F = (k_B * T / L_p) * [1/(4*(1-x/L_c)²) - 1/4 + x/L_c]
```

Where:
- L_p = persistence length = 7.5 nm (Rief et al. 1999)
- L_c = contour length = 194 nm (spectrin tetramer)
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

### Membrane Energy
```
W = (G_s/4) * (I₁² + 2*I₁ - 2*I₂) + (G_a/4) * I₂²
```

Where:
- G_s = shear modulus = 5.5 μN/m (Dao et al. 2003)
- G_a = area dilation modulus ≈ 500 μN/m
- I₁ = α² + β² - 2 (strain invariant 1)
- I₂ = α²β² - 1 (strain invariant 2)
- α, β = principal stretch ratios

### Bending Energy (Helfrich)
```
W_bend = (κ/2) * (2H - C₀)² + κ_G * K
```

Where:
- κ = bending modulus = 2×10⁻¹⁹ J (Evans 1983)
- H = mean curvature
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
