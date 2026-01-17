# /parameter - Biological Parameter Lookup

Look up, add, or modify biological parameters with required citations.

## Usage

```
/parameter <action> [args]
```

## Actions

### lookup
Find a parameter value and its source.

```
/parameter lookup <name>
/parameter lookup shear_modulus
/parameter lookup km_hexokinase
```

### add
Add a new parameter with citation (required).

```
/parameter add <name> <value> <unit> --citation "<Author Year>" --doi "<DOI>"
```

Example:
```
/parameter add km_pfk_f6p 0.1 mM --citation "Mulquiney 1999" --doi "10.1042/bj3420581"
```

### list
List all parameters in a category.

```
/parameter list enzymes
/parameter list membrane
/parameter list hemoglobin
/parameter list ions
```

### validate
Check all parameters have valid citations.

```
/parameter validate
```

## Parameter Categories

| Category | File | Contents |
|----------|------|----------|
| `enzymes` | `data/parameters/enzymes.json` | Km, Vmax, inhibitors |
| `metabolites` | `data/parameters/metabolites.json` | Concentrations |
| `membrane` | `data/parameters/membrane.json` | Mechanical properties |
| `hemoglobin` | `data/parameters/hemoglobin.json` | Oxygen binding |
| `ions` | `data/parameters/ions.json` | Channel conductances |

## Parameter Schema

```json
{
  "name": "shear_modulus",
  "value": 5.5,
  "unit": "μN/m",
  "uncertainty": 1.1,
  "citation": {
    "author": "Dao et al.",
    "year": 2003,
    "journal": "PNAS",
    "doi": "10.1073/pnas.0730616100"
  },
  "notes": "Measured via micropipette aspiration"
}
```

## Citation Requirements

All parameters MUST have:
1. Primary literature citation (not review/textbook)
2. DOI or PubMed ID when available
3. Measurement method noted

## Example Session

```
> /parameter lookup p50_standard

Name: p50_standard
Value: 26.8 mmHg
Uncertainty: ±1.0 mmHg
Citation: Imai (1982) "Allosteric Effects in Haemoglobin"
Method: Oxygen equilibrium at pH 7.4, 37°C

> /parameter add vmax_ldh 115 umol/min/mL --citation "Beutler 1984"

Added parameter:
  vmax_ldh = 115 μmol/min/mL
  Citation: Beutler (1984)
  WARNING: DOI not provided - please add for completeness
```
