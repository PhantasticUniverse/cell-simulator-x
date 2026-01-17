# Scientific Validation Domain Knowledge

## Validation Methodology

### Levels of Validation

1. **Unit Validation**: Individual equations match analytical solutions
2. **Component Validation**: Modules reproduce published model results
3. **System Validation**: Full simulation matches experimental data
4. **Sensitivity Analysis**: Parameter uncertainty quantification

### Statistical Metrics

**Coefficient of Determination:**
```
R² = 1 - SS_res / SS_tot
Target: R² > 0.9 for metabolic flux
```

**Root Mean Square Error:**
```
RMSE = sqrt(mean((predicted - observed)²))
Report in same units as measurement
```

**Percent Error:**
```
%Error = |predicted - observed| / observed * 100
Target: <10% for mechanical properties
```

## Mechanical Validation

### Micropipette Aspiration (Dao et al. 2003)

**Protocol:**
1. Apply aspiration pressure ΔP = 0-10 Pa
2. Measure aspiration length L_p
3. Calculate shear modulus from slope

**Acceptance Criteria:**
- Shear modulus: 5.5 ± 1.5 μN/m
- Linear response up to 50% strain

**Data Format:**
```json
{
  "experiment": "micropipette_aspiration",
  "source": "Dao et al. 2003",
  "data": [
    {"pressure_Pa": 1.0, "length_um": 0.5, "error": 0.1},
    {"pressure_Pa": 2.0, "length_um": 1.1, "error": 0.1}
  ]
}
```

### Optical Tweezers (Mills et al. 2004)

**Protocol:**
1. Attach beads to cell surface
2. Apply stretching force 0-200 pN
3. Measure axial/transverse diameters

**Acceptance Criteria:**
- Axial elongation curve within 10% of experimental

### Osmotic Fragility

**Protocol:**
1. Expose RBCs to hypotonic NaCl (0.1-0.9%)
2. Measure hemolysis fraction
3. Report as % hemolysis vs osmolarity

**Acceptance Criteria:**
- 50% hemolysis at ~0.45% NaCl
- Sigmoid curve shape

## Metabolic Validation

### Steady-State Concentrations

**Reference: Beutler 1984**

| Metabolite | Target (mM) | Tolerance |
|------------|-------------|-----------|
| ATP | 1.5-2.0 | ±20% |
| 2,3-DPG | 4.0-5.0 | ±15% |
| GSH | 2.0-2.5 | ±20% |
| Lactate | 1.0-2.0 | ±25% |

### Flux Validation

**Glucose Consumption:**
- Target: 1.2-1.5 μmol/hr/mL RBC
- Method: Compare to NMR/MS data

**Lactate Production:**
- Target: ~2× glucose consumption
- Validates glycolytic stoichiometry

**ATP Turnover:**
- Target: ~3 μmol/hr/mL RBC
- Balance: production = consumption

### Dynamic Response

**Glucose Step Test (Joshi-Palsson):**
1. Step glucose 5 → 10 mM
2. Monitor metabolite transients
3. Compare time constants

## Oxygen Transport Validation

### Oxygen Equilibrium Curve (OEC)

**Standard Conditions:**
- pH 7.4, T = 37°C, PCO2 = 40 mmHg
- P50 = 26.8 ± 1 mmHg (Imai 1982)
- Hill coefficient n = 2.7

**Data Points:**
```json
{
  "experiment": "oec_standard",
  "conditions": {"pH": 7.4, "temp_C": 37, "pCO2_mmHg": 40},
  "data": [
    {"pO2_mmHg": 10, "saturation": 0.13},
    {"pO2_mmHg": 20, "saturation": 0.35},
    {"pO2_mmHg": 26.8, "saturation": 0.50},
    {"pO2_mmHg": 40, "saturation": 0.75},
    {"pO2_mmHg": 60, "saturation": 0.89},
    {"pO2_mmHg": 100, "saturation": 0.97}
  ]
}
```

### Bohr Effect

**Test Protocol:**
1. Measure OEC at pH 7.2, 7.4, 7.6
2. Calculate ΔlogP50/ΔpH

**Acceptance:**
- Bohr coefficient = -0.48 ± 0.05

### 2,3-DPG Effect

**Test Protocol:**
1. Measure OEC at [DPG] = 0, 2, 4, 6 mM
2. Report P50 shift per mM DPG

**Acceptance:**
- P50 increases ~3 mmHg per mM DPG

## Integration Validation

### Deformation → ATP Release (Wan et al. 2008)

**Protocol:**
1. Apply shear stress 0-10 Pa
2. Measure extracellular ATP
3. Compare dose-response curve

### Volume Regulation

**Protocol:**
1. Osmotic challenge (300 → 200 mOsm)
2. Monitor volume recovery
3. Time constant ~minutes

## Experimental Data Format

### Standard JSON Schema
```json
{
  "experiment_id": "string",
  "source": {
    "citation": "Author et al. Year",
    "doi": "10.xxxx/xxxxx"
  },
  "conditions": {
    "temperature_C": 37,
    "pH": 7.4,
    "other_params": {}
  },
  "data": [
    {
      "independent_var": 0.0,
      "dependent_var": 0.0,
      "error": 0.0,
      "n_samples": 10
    }
  ],
  "metadata": {
    "method": "description",
    "units": {"independent": "unit", "dependent": "unit"}
  }
}
```

## Sensitivity Analysis

### One-at-a-Time (OAT)
- Vary each parameter ±10%
- Report output sensitivity coefficient

### Morris Method
- Global screening
- Identify influential parameters

### Sobol Indices
- Variance-based global sensitivity
- First-order and total effects

## Key Citations

- Dao et al. (2003) PNAS - Mechanical validation
- Joshi & Palsson (1989) J Theor Biol - Metabolic validation
- Imai (1982) - Oxygen binding validation
- Wan et al. (2008) - Mechano-ATP release
- Beutler (1984) - Metabolite reference values
