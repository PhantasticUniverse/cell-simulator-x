# RBC Biochemistry Domain Knowledge

## Enzyme Kinetics Framework

### Michaelis-Menten Kinetics
Standard form for most RBC enzymes:
```
v = Vmax * [S] / (Km + [S])
```

### Hill Kinetics (Cooperative Binding)
For enzymes with sigmoidal kinetics (PFK, PK):
```
v = Vmax * [S]^n / (K0.5^n + [S]^n)
```

### Key Glycolytic Enzymes

| Enzyme | Km (mM) | Vmax | Regulation |
|--------|---------|------|------------|
| Hexokinase (HK) | 0.04 (glucose) | 2.8 | Inhibited by G6P |
| PFK | 0.1 (F6P) | 92 | Activated by ADP, Pi; Inhibited by ATP, citrate |
| Pyruvate Kinase | 0.15 (PEP) | 286 | Activated by F1,6BP |
| LDH | 0.14 (pyruvate) | 115 | Near equilibrium |

## Hemoglobin Oxygen Binding

### Adair Equation (4-site binding)
```
Y = (K1*pO2 + 2*K1*K2*pO2^2 + 3*K1*K2*K3*pO2^3 + 4*K1*K2*K3*K4*pO2^4) /
    (4 * (1 + K1*pO2 + K1*K2*pO2^2 + K1*K2*K3*pO2^3 + K1*K2*K3*K4*pO2^4))
```

### Allosteric Modifiers

**Bohr Effect (pH)**
- ΔlogP50/ΔpH = -0.48 (Imai 1982)
- T-state stabilized by H+ binding to His146β

**2,3-DPG Effect**
- Binds T-state Hb in β-chain central cavity
- Increases P50 by ~10 mmHg per mM
- Normal: 4-5 mM

**Temperature**
- ΔH = -14.5 kcal/mol O2 (Imai 1982)
- Higher T → right shift

## Glycolysis Pathway

```
Glucose → G6P → F6P → F1,6BP → DHAP + GA3P
                                    ↓
                              1,3-BPG → 3PG → 2PG → PEP → Pyruvate
                                  ↓
                              2,3-DPG (Rapoport-Luebering shunt)
```

### ATP Yield
- 2 ATP consumed (HK, PFK)
- 4 ATP produced (PGK, PK)
- Net: 2 ATP per glucose

### 2,3-DPG Regulation
- BPGM: Converts 1,3-BPG → 2,3-DPG
- BPGP: Hydrolyzes 2,3-DPG → 3-PG
- Balance determines Hb oxygen affinity

## Pentose Phosphate Pathway (PPP)

### Oxidative Branch
```
G6P → 6PG → Ribulose-5P + CO2
      ↓         ↓
   NADPH     NADPH
```

### Key Enzymes
| Enzyme | Km | Vmax | Regulation |
|--------|-----|------|------------|
| G6PDH | 0.007 mM (G6P) | 0.08 mM/s | Inhibited by NADPH (Ki 0.005 mM) |
| 6PGDH | 0.023 mM (6PG) | 0.04 mM/s | Product inhibition |

### Physiological Role
- Provides NADPH for glutathione reduction
- Normal PPP flux: 3-11% of glycolysis (model: ~60% due to high demand)
- NADPH/NADP+ ratio target: 10-20

## Glutathione Cycle

### Redox Reactions
```
H2O2 + 2GSH → GSSG + 2H2O  (GPx)
GSSG + NADPH + H+ → 2GSH + NADP+  (GR)
```

### Key Enzymes
| Enzyme | Function | Key Parameters |
|--------|----------|----------------|
| GPx | H2O2 detoxification | Vmax = 0.02 mM/s, Km_H2O2 = 0.002 mM, Km_GSH = 1.0 mM |
| GR | GSH regeneration | Vmax = 0.15 mM/s, Km_GSSG = 0.015 mM, Km_NADPH = 0.015 mM |

### Steady-State Targets
- GSH/GSSG ratio: >50 (model achieves ~2500 due to efficient GR)
- Total glutathione: 2-3 mM
- H2O2: <5 µM (basal production 5 µM/s in model)

## Ion Transport

### Na+/K+-ATPase
- 3 Na+ out, 2 K+ in per ATP
- Maintains membrane potential
- Rate: ~3 μmol/hr/mL RBC

### Piezo1 (Mechanosensitive)
- Hill-type activation by membrane tension
- Half-activation: ~1.5 pN/nm
- Hill coefficient: ~3
- Ca2+ conductance: ~10 pS
- Triggers:
  - Ca2+ influx (Ca_ext ~1.8 mM → cytosol)
  - Downstream: Gardos channel, ATP release
  - Volume regulation

### Band 3 (AE1)
- Cl-/HCO3- exchanger
- Facilitates CO2 transport
- 1.2 million copies/cell

## Metabolite Concentrations (Reference)

| Metabolite | Normal (mM) | Source |
|------------|-------------|--------|
| ATP | 1.5-2.0 | Beutler 1984 |
| ADP | 0.2-0.4 | Beutler 1984 |
| AMP | 0.01-0.02 | Beutler 1984 |
| 2,3-DPG | 4.0-5.0 | Benesch 1969 |
| GSH | 2.0-2.5 | Beutler 1984 |
| GSSG | 0.02-0.05 | Beutler 1984 |
| NAD+ | 0.05-0.07 | Beutler 1984 |
| NADH | 0.003-0.005 | Beutler 1984 |
| NADPH | 0.05-0.1 | Kirkman 2007 |
| NADP+ | 0.003-0.01 | Kirkman 2007 |
| H2O2 | <0.005 (5 µM) | Chance 1979 |
| Pi | 1.0-1.5 | Beutler 1984 |

## Key Citations

- Joshi & Palsson (1989) J Theor Biol 141:515 - Core RBC metabolic model
- Mulquiney et al. (1999) Biochem J 342:581 - pH-dependent kinetics
- Imai (1982) "Allosteric Effects in Haemoglobin" - Oxygen binding
- Beutler (1984) "Red Cell Metabolism" - Metabolite reference
- Kirkman & Gaetani (2007) Trends Biochem Sci - NADPH homeostasis
- Wu et al. (2004) J Biol Chem - Glutathione redox
- Chance (1979) Meth Enzymol - H2O2 measurement
