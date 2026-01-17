# Cell Simulator X: Human Erythrocyte Simulation Engine

## Product Requirements Document (PRD)

---

## 1. Executive Summary

**Cell Simulator X** is a cutting-edge, GPU-accelerated simulation of the human red blood cell (erythrocyte) that integrates mechanics, metabolism, and molecular biology into a unified computational model. This will be the **first fully integrated mechano-metabolic cell simulation**, surpassing all existing RBC models by coupling:

- Complete metabolic network (~30 pathways, 200+ reactions)
- Cytoskeleton mechanics (spectrin network dynamics)
- Membrane deformation and lipid bilayer physics
- Hemoglobin oxygen transport with full allosteric effects
- Ion homeostasis and osmotic regulation
- Mechanosensitive processes (Piezo1-mediated ATP release)

**Why this matters**: No existing simulation integrates RBC mechanics and metabolism. Current models treat these as separate domains, missing critical feedback loops that govern real cell behavior.

---

## 2. Background & Market Analysis

### 2.1 Current State of RBC Simulation

| Domain | Leading Work | Key Limitation |
|--------|--------------|----------------|
| **Mechanics** | Fedosov et al. (2010) - DPD multiscale model | No metabolism, treats cell as passive material |
| **Metabolism** | Yachie-Kinoshita (2010) - E-Cell model | Only 3 pathways, no spatial resolution, no mechanics |
| **Oxygen Transport** | MWC/Adair models | Not integrated with cellular ATP/2,3-DPG dynamics |
| **Proteomics** | Bryk & Bhavnani (2022) - 1,202 proteins | Static map, no dynamic simulation |

### 2.2 The Gap We Fill

**No model couples deformation → mechanosensing → ATP release → metabolism → 2,3-DPG → oxygen affinity → tissue oxygenation**

This feedback loop is critical for understanding:
- Blood storage lesion (transfusion medicine)
- Malaria pathophysiology (P. falciparum changes RBC mechanics)
- Sickle cell disease (HbS polymerization + metabolism)
- Diabetic microvascular complications

### 2.3 Data Availability Assessment

| Data Type | Source | Coverage | Quality |
|-----------|--------|----------|---------|
| Proteome | Human Protein Atlas, Bryk et al. 2022 | 1,202 proteins | High |
| Metabolic kinetics | BRENDA, SABIO-RK, Joshi-Palsson 1989 | ~80% of enzymes | Medium-High |
| Membrane mechanics | Dao et al., Suresh group | Comprehensive | High |
| Hemoglobin kinetics | Adair, Imai, Weber | Complete | High |
| Ion channels | Bhavnani et al., Bhalla-Bhagwani | Good | Medium |
| Spectrin structure | Cryo-ET (Grigorieff), AFM | Improving | Medium |

---

## 3. Product Vision & Goals

### 3.1 Vision Statement

> Create the most accurate, comprehensive, and usable simulation of a human cell ever built, starting with the erythrocyte as a tractable yet medically important target.

### 3.2 Primary Goals

1. **Scientific Accuracy**: Every parameter derived from peer-reviewed experimental data
2. **Integration**: First model to couple mechanics ↔ metabolism ↔ oxygen transport
3. **Performance**: Real-time visualization on consumer hardware (M4 Max baseline)
4. **Extensibility**: Architecture supports future cell types
5. **Validation**: Quantitative agreement with published experimental data

### 3.3 Success Metrics

| Metric | Target | Validation Method |
|--------|--------|-------------------|
| Metabolic flux accuracy | R² > 0.9 vs experimental | Compare to NMR/MS metabolomics |
| Mechanical deformation | <10% error vs micropipette | Reproduce Dao et al. experiments |
| Oxygen equilibrium curve | <5% P50 error across conditions | Match Imai data (pH, T, DPG) |
| ATP concentration | Within physiological range (1.5-2.5 mM) | Steady-state validation |
| Simulation speed | >10 fps at full resolution | Benchmarking |

---

## 4. Technical Architecture

### 4.1 Core Simulation Engine

```
┌─────────────────────────────────────────────────────────────────┐
│                    CELL SIMULATOR X ENGINE                       │
├─────────────────────────────────────────────────────────────────┤
│  ┌──────────────┐  ┌──────────────┐  ┌──────────────────────┐  │
│  │   GEOMETRY   │  │   PHYSICS    │  │    BIOCHEMISTRY      │  │
│  │    MODULE    │  │    MODULE    │  │       MODULE         │  │
│  │              │  │              │  │                      │  │
│  │ • Mesh gen   │  │ • Membrane   │  │ • Metabolic network  │  │
│  │ • Spectrin   │  │   mechanics  │  │ • Enzyme kinetics    │  │
│  │   network    │  │ • Fluid      │  │ • Hemoglobin binding │  │
│  │ • Protein    │  │   dynamics   │  │ • Ion transport      │  │
│  │   placement  │  │ • DPD solver │  │ • ODE/SSA solver     │  │
│  └──────┬───────┘  └──────┬───────┘  └──────────┬───────────┘  │
│         │                 │                      │              │
│         └─────────────────┼──────────────────────┘              │
│                           │                                     │
│                    ┌──────▼───────┐                             │
│                    │  INTEGRATION │                             │
│                    │    LAYER     │                             │
│                    │              │                             │
│                    │ • State sync │                             │
│                    │ • Event bus  │                             │
│                    │ • Coupling   │                             │
│                    └──────┬───────┘                             │
│                           │                                     │
│         ┌─────────────────┼─────────────────────┐               │
│         │                 │                     │               │
│  ┌──────▼───────┐  ┌──────▼───────┐  ┌─────────▼────────┐      │
│  │     GPU      │  │   RENDER     │  │    ANALYSIS      │      │
│  │   COMPUTE    │  │    ENGINE    │  │     ENGINE       │      │
│  │              │  │              │  │                  │      │
│  │ • Metal/CUDA │  │ • WebGPU     │  │ • Time series    │      │
│  │ • Parallel   │  │ • Real-time  │  │ • Flux analysis  │      │
│  │   solvers    │  │   viz        │  │ • Sensitivity    │      │
│  └──────────────┘  └──────────────┘  └──────────────────┘      │
└─────────────────────────────────────────────────────────────────┘
```

### 4.2 Module Specifications

#### 4.2.1 Geometry Module

**Purpose**: Define the 3D structure of the RBC including membrane, cytoskeleton, and protein localization.

**Data Sources**:
- Biconcave disc geometry: Fung-Tong parametric equations
- Spectrin network topology: Cryo-ET data (Liu et al.)
- Protein localization: Bryk & Bhavnani 2022 interactome

**Key Components**:
- `RBCMesh`: Triangulated surface mesh (~10,000 vertices)
- `SpectrinNetwork`: Graph of spectrin tetramers + junctional complexes
- `ProteinPlacement`: Spatial distribution of 1,202 proteins
- `MembraneCompartments`: Lipid raft domains, Band3 clusters

#### 4.2.2 Physics Module

**Purpose**: Simulate mechanical behavior including deformation, flow, and membrane dynamics.

**Methods**:
- Dissipative Particle Dynamics (DPD) for fluid-membrane interaction
- Worm-like chain (WLC) model for spectrin elasticity
- Skalak strain energy function for membrane area/shear

**Key Parameters** (from literature):
| Parameter | Value | Source |
|-----------|-------|--------|
| Membrane shear modulus | 5.5 μN/m | Dao et al. 2003 |
| Bending modulus | 2×10⁻¹⁹ J | Evans 1983 |
| Spectrin contour length | 194 nm | Rief et al. 1999 |
| Persistence length | 7.5 nm | Rief et al. 1999 |
| Cytoplasmic viscosity | 6 cP | Cokelet & Meiselman |

#### 4.2.3 Biochemistry Module

**Purpose**: Simulate all metabolic reactions, enzyme kinetics, and molecular interactions.

**Metabolic Pathways Included**:
1. **Glycolysis** (Embden-Meyerhof pathway) - 10 reactions
2. **Rapoport-Luebering shunt** (2,3-DPG) - 2 reactions
3. **Pentose Phosphate Pathway** - 7 reactions
4. **Glutathione metabolism** - 4 reactions
5. **Nucleotide salvage** - 8 reactions
6. **Adenylate kinase equilibrium** - 1 reaction
7. **Hexosamine pathway** - 5 reactions
8. **Polyol pathway** - 2 reactions

**Hemoglobin Submodule**:
- 4-site Adair equation with allosteric modifiers
- Bohr effect (pH dependence)
- Temperature dependence
- 2,3-DPG binding to deoxy-Hb
- CO₂ binding (carbamino formation)
- Methemoglobin reduction

**Ion Transport**:
- Na⁺/K⁺-ATPase (3Na out, 2K in, 1 ATP consumed)
- Ca²⁺-ATPase (PMCA)
- Gardos channel (Ca²⁺-activated K⁺)
- Piezo1 (mechanosensitive cation channel)
- Band3 (Cl⁻/HCO₃⁻ exchanger)
- Aquaporin-1 (water)
- GLUT1 (glucose)

**Kinetic Framework**:
- Michaelis-Menten kinetics for most enzymes
- Hill kinetics for cooperative binding (PFK, PK)
- Mass action for reversible reactions
- Stochastic simulation (Gillespie SSA) for low-copy species

#### 4.2.4 Integration Layer

**Critical Coupling Points**:

1. **Deformation → ATP Release**
   - Piezo1 activation by membrane tension
   - ATP release proportional to deformation magnitude
   - Ca²⁺ influx triggering downstream effects

2. **Metabolism → Mechanics**
   - ATP levels affect spectrin phosphorylation
   - 2,3-DPG affects Hb conformation and cell volume
   - GSH levels affect membrane lipid peroxidation

3. **Oxygen → Metabolism**
   - pO₂ affects Hb quaternary state
   - T-state Hb sequesters 2,3-DPG
   - Oxygen levels affect oxidative stress pathways

### 4.3 Compute Architecture

**Primary**: Metal (Apple Silicon) - for M4 Max development
**Secondary**: CUDA (NVIDIA) - for workstation/cloud deployment
**Future**: WebGPU - for browser-based access

**Parallelization Strategy**:
- Mesh vertices: parallel position updates
- Metabolic reactions: batch ODE integration
- Stochastic events: parallel random number generation
- Rendering: GPU instancing for protein visualization

### 4.4 Data Structures

```rust
struct CellState {
    geometry: GeometryState,
    biochemistry: BiochemistryState,
    physics: PhysicsState,
    environment: EnvironmentState,
}

struct GeometryState {
    vertices: Vec<[f32; 3]>,           // Mesh positions (N×3)
    velocities: Vec<[f32; 3]>,         // Mesh velocities (N×3)
    spectrin_nodes: Vec<[f32; 3]>,     // Junction positions (M×3)
    spectrin_edges: Vec<[u32; 2]>,     // Network topology (E×2)
}

struct BiochemistryState {
    metabolites: HashMap<String, f64>,        // Concentrations (mM)
    enzymes: HashMap<String, EnzymeState>,    // Activity states
    hemoglobin: HemoglobinState,
    ions: IonState,
}

struct HemoglobinState {
    total: f64,              // Total Hb (mM heme)
    saturation: f64,         // O₂ saturation (0-1)
    t_state_fraction: f64,   // Tense state fraction
    met_hb_fraction: f64,    // Methemoglobin fraction
}

struct IonState {
    na_in: f64,
    k_in: f64,
    cl_in: f64,
    ca_in: f64,
    mg_in: f64,
    h_in: f64,  // For pH calculation
}

struct PhysicsState {
    volume: f64,                      // Cell volume (fL)
    surface_area: f64,                // Surface area (μm²)
    membrane_tension: Vec<f64>,       // Local tension per vertex
    temperature: f64,                 // Kelvin
}

struct EnvironmentState {
    p_o2: f64,        // O₂ partial pressure (mmHg)
    p_co2: f64,       // CO₂ partial pressure (mmHg)
    glucose: f64,     // External glucose (mM)
    ph_ext: f64,      // External pH
    osmolarity: f64,  // External osmolarity (mOsm)
}
```

---

## 5. Data Sources & Parameters

### 5.1 Primary Data Sources

| Source | Data Type | URL/Reference |
|--------|-----------|---------------|
| Human Protein Atlas | Proteome | proteinatlas.org |
| BRENDA | Enzyme kinetics | brenda-enzymes.org |
| SABIO-RK | Reaction kinetics | sabiork.h-its.org |
| BioModels | Existing models | ebi.ac.uk/biomodels |
| PDB | Protein structures | rcsb.org |
| UniProt | Protein sequences | uniprot.org |
| Joshi-Palsson 1989 | RBC metabolic model | J Theor Biol 141:515 |
| Mulquiney et al. 1999 | pH-dependent kinetics | Biochem J 342:581 |

### 5.2 Key Parameter Values

#### Metabolite Concentrations (mM, normal RBC)
| Metabolite | Concentration | Source |
|------------|---------------|--------|
| ATP | 1.5-2.0 | Beutler 1984 |
| ADP | 0.2-0.4 | Beutler 1984 |
| 2,3-DPG | 4.0-5.0 | Benesch 1969 |
| GSH | 2.0-2.5 | Beutler 1984 |
| NAD⁺ | 0.05-0.07 | Beutler 1984 |
| NADH | 0.003-0.005 | Beutler 1984 |
| Glucose-6-P | 0.03-0.05 | Beutler 1984 |
| Lactate | 1.0-2.0 | Beutler 1984 |

#### Hemoglobin Parameters
| Parameter | Value | Source |
|-----------|-------|--------|
| Total Hb | 5.0 mM (heme) | Imai 1982 |
| P50 (standard) | 26.8 mmHg | Imai 1982 |
| Hill coefficient | 2.7 | Imai 1982 |
| Bohr coefficient | -0.48 | Imai 1982 |
| ΔH (binding) | -14.5 kcal/mol | Imai 1982 |

---

## 6. Use Cases & User Stories

### 6.1 Research Tool

**US-R1**: As a researcher, I want to perturb enzyme activities to understand metabolic control, so I can identify rate-limiting steps.

**US-R2**: As a biophysicist, I want to simulate RBC deformation under shear flow to understand mechanobiology.

**US-R3**: As a hematologist, I want to model how storage time affects RBC function to improve blood banking.

### 6.2 Drug Discovery

**US-D1**: As a pharma researcher, I want to simulate antimalarial effects on RBC metabolism to prioritize drug candidates.

**US-D2**: As a toxicologist, I want to predict hemolytic potential of compounds.

### 6.3 Disease Modeling

**US-DM1**: As a sickle cell researcher, I want to model HbS polymerization effects on cell mechanics.

**US-DM2**: As a diabetologist, I want to simulate how hyperglycemia affects RBC function over time.

### 6.4 Education

**US-E1**: As a student, I want to visualize how oxygen binds to hemoglobin to understand cooperative binding.

**US-E2**: As a teacher, I want to demonstrate metabolic flux changes under hypoxia.

---

## 7. Implementation Phases

### Phase 1: Foundation (Months 1-3)
**Goal**: Basic infrastructure and geometry

- [ ] Project setup (Rust/Metal/WebGPU toolchain)
- [ ] RBC mesh generation with Fung-Tong parametric surface
- [ ] Basic rendering pipeline with camera controls
- [ ] Spectrin network graph structure
- [ ] Core data structures (CellState)
- [ ] Configuration system for parameters

**Deliverable**: Rotating 3D RBC visualization with spectrin network overlay

### Phase 2: Mechanics (Months 4-6)
**Goal**: Accurate mechanical simulation

- [ ] DPD fluid solver (GPU-accelerated)
- [ ] Membrane mechanics (Skalak model)
- [ ] Spectrin elasticity (WLC model)
- [ ] Deformation validation vs micropipette aspiration data
- [ ] Shear flow simulation
- [ ] Osmotic swelling/shrinking

**Deliverable**: Deformable RBC that matches experimental mechanics data

### Phase 3: Core Metabolism (Months 7-9)
**Goal**: Complete glycolytic pathway with validation

- [ ] ODE solver (adaptive Runge-Kutta)
- [ ] Glycolysis implementation (all 10 reactions)
- [ ] Rapoport-Luebering shunt (2,3-DPG)
- [ ] ATP consumption/production balance
- [ ] Steady-state validation
- [ ] Perturbation response validation (glucose step)

**Deliverable**: Validated metabolic model reproducing Joshi-Palsson results

### Phase 4: Oxygen Transport (Months 10-11)
**Goal**: Hemoglobin-oxygen dynamics

- [ ] Adair equation implementation
- [ ] Bohr effect (pH dependence)
- [ ] 2,3-DPG allosteric effect
- [ ] Temperature dependence
- [ ] Oxygen equilibrium curve validation
- [ ] Dynamic oxygen uptake/release

**Deliverable**: Accurate OEC across all physiological conditions

### Phase 5: Integration (Months 12-14)
**Goal**: Couple mechanics and biochemistry

- [ ] Piezo1 mechanosensitive channel
- [ ] ATP release under deformation
- [ ] Ion homeostasis (Na/K-ATPase, etc.)
- [ ] Volume regulation feedback
- [ ] Integrated state synchronization
- [ ] Performance optimization

**Deliverable**: First integrated mechano-metabolic RBC model

### Phase 6: Extended Biochemistry (Months 15-17)
**Goal**: Complete metabolic network

- [ ] Pentose phosphate pathway
- [ ] Glutathione redox cycle
- [ ] Nucleotide metabolism
- [ ] NADPH/NADH coupling
- [ ] Oxidative stress simulation
- [ ] MetHb formation and reduction

**Deliverable**: Complete metabolic coverage

### Phase 7: Disease Models (Months 18-20)
**Goal**: Pathological simulations

- [ ] Malaria (P. falciparum metabolic takeover)
- [ ] Sickle cell (HbS polymerization)
- [ ] Storage lesion (blood banking)
- [ ] Diabetic RBC changes
- [ ] Validation against clinical data

**Deliverable**: Clinically relevant disease simulations

### Phase 8: Polish & Release (Months 21-24)
**Goal**: Production-ready software

- [ ] Comprehensive documentation
- [ ] Parameter database with citations
- [ ] User interface refinement
- [ ] Educational mode
- [ ] API for external tools
- [ ] Publication of methods paper

**Deliverable**: v1.0 release

---

## 8. Technical Requirements

### 8.1 Development Environment
- **Language**: Rust (core engine), TypeScript (UI)
- **GPU**: Metal (primary), CUDA (secondary), WebGPU (web)
- **Build**: Cargo, wasm-pack
- **Testing**: Rust unit tests, integration benchmarks

### 8.2 Performance Targets
| Metric | Target | Notes |
|--------|--------|-------|
| Frame rate | 60 fps | Visualization mode |
| Simulation step | <1ms | 10,000 mesh vertices |
| Memory | <4 GB | Full cell state |
| Startup | <5 sec | Initial mesh generation |

### 8.3 Validation Requirements
- Each module validated independently before integration
- Quantitative comparison to published experimental data
- Sensitivity analysis for uncertain parameters
- Unit tests for all mathematical functions

---

## 9. Risks & Mitigations

| Risk | Likelihood | Impact | Mitigation |
|------|------------|--------|------------|
| Missing kinetic parameters | High | Medium | Use estimation methods, mark uncertain values |
| Numerical instability | Medium | High | Adaptive time-stepping, careful solver choice |
| GPU compatibility issues | Medium | Medium | Abstract compute backend, fallback to CPU |
| Validation data mismatch | Medium | High | Literature review, expert consultation |
| Scope creep | High | Medium | Strict phase gates, MVP focus |

---

## 10. Validation Strategy

### 10.1 Mechanical Validation
- Micropipette aspiration: Match Dao et al. force-extension curves
- Optical tweezers: Reproduce stretching experiments
- Osmotic fragility: Match hemolysis vs osmolarity curves

### 10.2 Metabolic Validation
- Steady-state metabolite concentrations vs Beutler 1984
- Glucose consumption rate: 1.2-1.5 μmol/hr/mL RBC
- Lactate production: ~2× glucose consumption
- ATP turnover: ~3 μmol/hr/mL RBC

### 10.3 Oxygen Transport Validation
- P50 at standard conditions: 26.8 ± 1 mmHg
- Bohr shift: ΔlogP50/ΔpH = -0.48
- Temperature coefficient: ΔlogP50/Δ(1/T) as per Imai
- 2,3-DPG effect: quantitative shift per mM DPG

### 10.4 Integration Validation
- ATP release under deformation (Wan et al. 2008)
- Volume response to osmotic challenge
- Recovery kinetics after perturbation

---

## 11. File Structure

```
cell-simulator-x/
├── Cargo.toml
├── README.md
├── PRD.md
├── docs/
│   ├── architecture.md
│   └── parameters/
│       ├── metabolites.json
│       ├── enzymes.json
│       └── membrane.json
├── src/
│   ├── main.rs
│   ├── lib.rs
│   ├── geometry/
│   │   ├── mod.rs
│   │   ├── mesh.rs
│   │   ├── spectrin.rs
│   │   └── proteins.rs
│   ├── physics/
│   │   ├── mod.rs
│   │   ├── dpd.rs
│   │   ├── membrane.rs
│   │   └── fluid.rs
│   ├── biochemistry/
│   │   ├── mod.rs
│   │   ├── metabolism.rs
│   │   ├── hemoglobin.rs
│   │   ├── ions.rs
│   │   └── enzymes/
│   │       ├── glycolysis.rs
│   │       ├── ppp.rs
│   │       └── ...
│   ├── integration/
│   │   ├── mod.rs
│   │   ├── coupling.rs
│   │   └── state.rs
│   ├── compute/
│   │   ├── mod.rs
│   │   ├── metal.rs
│   │   ├── cuda.rs
│   │   └── cpu_fallback.rs
│   ├── render/
│   │   ├── mod.rs
│   │   ├── camera.rs
│   │   └── shaders/
│   └── analysis/
│       ├── mod.rs
│       ├── timeseries.rs
│       └── validation.rs
├── shaders/
│   ├── cell.wgsl
│   ├── spectrin.wgsl
│   └── compute/
├── data/
│   ├── parameters/
│   └── validation/
└── tests/
    ├── mechanics_tests.rs
    ├── metabolism_tests.rs
    └── integration_tests.rs
```

---

## 12. Immediate Next Steps

1. **Set up Rust project** with Metal compute backend
2. **Implement RBC mesh** using Fung-Tong parametric equations
3. **Create basic renderer** with wgpu
4. **Build parameter loading system** from JSON
5. **Implement first enzyme** (hexokinase) as proof of concept

---

## 13. References

### Core Papers
1. Joshi & Palsson (1989) "Metabolic dynamics in the human red cell" - J Theor Biol
2. Dao et al. (2003) "Mechanics of the human red blood cell" - PNAS
3. Fedosov et al. (2010) "Multiscale RBC model" - Biophys J
4. Bryk & Bhavnani (2022) "Protein organization of a red blood cell" - Cell Reports
5. Imai (1982) "Allosteric Effects in Haemoglobin" - Cambridge University Press

### Data Sources
6. BRENDA enzyme database - brenda-enzymes.org
7. Human Protein Atlas - proteinatlas.org
8. BioModels Database - ebi.ac.uk/biomodels
9. Beutler (1984) "Red Cell Metabolism" - Grune & Stratton

---

*Document Version: 1.0*
*Last Updated: 2026-01-17*
*Author: Cell Simulator X Team*
