# /cite - Add Literature Citation

Add or look up literature citations for the project.

## Usage

```
/cite <action> [args]
```

## Actions

### add
Add a new citation to the project bibliography.

```
/cite add --doi "10.1073/pnas.0730616100"
/cite add --pmid "12456789"
/cite add --manual "Author et al. (Year) Title. Journal Vol:Pages"
```

Auto-fetches metadata from DOI/PubMed when available.

### lookup
Find citations by keyword or author.

```
/cite lookup hemoglobin
/cite lookup "Dao"
/cite lookup 2003
```

### format
Format a citation in project style.

```
/cite format DOI
```

### verify
Verify all citations in codebase are in bibliography.

```
/cite verify
```

### export
Export bibliography in various formats.

```
/cite export bibtex > references.bib
/cite export markdown > REFERENCES.md
```

## Citation Format

Project uses inline citations: `(Author Year)` or `Author et al. Year`

Full citation stored in `docs/bibliography.json`:

```json
{
  "id": "dao2003",
  "authors": ["M. Dao", "C.T. Lim", "S. Suresh"],
  "year": 2003,
  "title": "Mechanics of the human red blood cell deformed by optical tweezers",
  "journal": "PNAS",
  "volume": 100,
  "pages": "7263-7268",
  "doi": "10.1073/pnas.0730616100",
  "tags": ["mechanics", "validation", "optical-tweezers"]
}
```

## In-Code Citation

When adding parameters or equations, cite inline:

```rust
/// Membrane shear modulus (Dao et al. 2003)
const SHEAR_MODULUS_UN_M: f64 = 5.5;

/// Adair equation for O2 binding (Imai 1982)
fn oxygen_saturation(p_o2: f64, params: &HbParams) -> f64 {
    // ...
}
```

## Example Session

```
> /cite add --doi "10.1016/j.bpj.2010.05.008"

Added citation:
  Fedosov et al. (2010)
  "Multiscale Modeling of Red Blood Cell Mechanics and Blood Flow in Malaria"
  Biophysical Journal
  DOI: 10.1016/j.bpj.2010.05.008

  Use as: (Fedosov et al. 2010) or Fedosov2010

> /cite lookup spectrin

Found 3 citations:
  1. Rief et al. (1999) - Spectrin unfolding mechanics
  2. Liu et al. (2006) - Spectrin network structure
  3. Fedosov et al. (2010) - Spectrin network model
```

## Required Fields

For publication, all citations must have:
- DOI (preferred) or PubMed ID
- Full author list
- Journal/Book title
- Year
- Volume/Pages (for journals)

## Key References

See PRD.md Section 13 for core citations.
