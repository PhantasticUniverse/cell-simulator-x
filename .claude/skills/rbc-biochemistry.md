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

## Ion Transport

### Na+/K+-ATPase
- 3 Na+ out, 2 K+ in per ATP
- Maintains membrane potential
- Rate: ~3 μmol/hr/mL RBC

### Piezo1 (Mechanosensitive)
- Activated by membrane tension
- Conducts Ca2+, triggers:
  - ATP release
  - Gardos channel activation
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
| NAD+ | 0.05-0.07 | Beutler 1984 |
| NADH | 0.003-0.005 | Beutler 1984 |
| Pi | 1.0-1.5 | Beutler 1984 |

## Key Citations

- Joshi & Palsson (1989) J Theor Biol 141:515 - Core RBC metabolic model
- Mulquiney et al. (1999) Biochem J 342:581 - pH-dependent kinetics
- Imai (1982) "Allosteric Effects in Haemoglobin" - Oxygen binding
- Beutler (1984) "Red Cell Metabolism" - Metabolite reference
