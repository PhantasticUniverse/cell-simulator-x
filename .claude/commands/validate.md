# /validate - Run Validation Tests

Run validation tests against experimental data to verify simulation accuracy.

## Usage

```
/validate [module] [--verbose]
```

## Arguments

- `module` (optional): Specific module to validate
  - `mechanics` - Membrane mechanics, micropipette aspiration
  - `metabolism` - Steady-state metabolites, flux rates
  - `oxygen` - Oxygen equilibrium curve, Bohr effect
  - `integration` - Coupled mechano-metabolic responses
  - (none) - Run all validation suites

- `--verbose`: Show detailed comparison data

## Validation Suites

### Mechanics Validation
- Micropipette aspiration force-extension (Dao et al. 2003)
- Osmotic fragility curve
- Shear modulus: target 5.5 ± 1.5 μN/m

### Metabolism Validation
- Steady-state concentrations vs Beutler 1984
- Glucose consumption rate: 1.2-1.5 μmol/hr/mL
- ATP turnover balance

### Oxygen Validation
- OEC at standard conditions (P50 = 26.8 mmHg)
- Bohr coefficient = -0.48
- 2,3-DPG effect quantification

### Integration Validation
- ATP release under deformation
- Volume regulation kinetics

## Output

Reports:
1. Pass/Fail for each validation criterion
2. R² or % error vs experimental data
3. Plots saved to `validation/results/`

## Example

```
/validate mechanics --verbose
```

Output:
```
=== Mechanics Validation ===
Micropipette aspiration:
  Shear modulus: 5.3 μN/m (target: 5.5 ± 1.5) ✓
  R² vs Dao 2003: 0.94 ✓

Osmotic fragility:
  50% hemolysis: 0.44% NaCl (target: 0.45%) ✓

PASSED: 3/3 tests
```

## Implementation

Runs: `cargo test --features validation -- validation::{module}`

Data files: `data/validation/*.json`
Results: `validation/results/`
