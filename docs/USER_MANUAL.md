# User Manual

Complete guide to using Cell Simulator X.

## Table of Contents

1. [GUI Mode](#gui-mode)
2. [CLI Diagnostic Modes](#cli-diagnostic-modes)
3. [Disease Modeling](#disease-modeling)
4. [Parameter Configuration](#parameter-configuration)
5. [Interpreting Output](#interpreting-output)
6. [Common Workflows](#common-workflows)

---

## GUI Mode

Launch the interactive 3D viewer:

```bash
cargo run --release
```

### Controls

| Key | Action | Description |
|-----|--------|-------------|
| **Mouse drag** | Orbit camera | Click and drag to rotate view |
| **S** | Toggle spectrin | Show/hide spectrin network overlay |
| **R** | Reset camera | Return to default view |
| **P** | Toggle physics | Start/stop physics simulation |
| **F** | Apply force | Apply test force to cell |
| **+** / **=** | Increase substeps | More physics substeps per frame |
| **-** | Decrease substeps | Fewer physics substeps per frame |
| **H** | Toggle help | Display/hide keybindings overlay |
| **E** | Toggle export menu | Show/hide export options |
| **M** | Toggle metabolites | Show/hide metabolites panel |
| **D** | Toggle disease | Show/hide disease panel (when active) |
| **Tab** | Toggle HUD | Show/hide entire HUD overlay |
| **F12** | Export state | Save current state to JSON |
| **Esc** | Exit | Close application |

### HUD Panels

The GUI displays real-time simulation data in two main panels:

**SIMULATION Panel (top-left)**:
- Time: Current simulation time
- FPS: Rendering frame rate
- Steps: Physics steps per second
- Mode: Normal / Stressed / Disease
- Physics: ON/OFF status
- Substeps: Physics substeps per frame

**METABOLITES Panel (top-right)**:
- **Energy**: ATP with status bar, 2,3-DPG with status bar
- **Oxygen**: O₂ saturation gauge (circular), pO₂, P50
- **Acid-Base**: Intracellular pH
- **Redox**: NADPH/NADP⁺, GSH/GSSG, H₂O₂
- **Ions**: Na⁺, K⁺, Ca²⁺ concentrations

**DISEASE Panel (bottom-right, when active)**:
- Disease name and severity
- Affected parameters with modifiers

Status indicators use color coding:
- 🟢 Green: Normal (within physiological range)
- 🟡 Yellow: Warning (approaching limits)
- 🔴 Red: Critical (outside safe range)

### 3D View

The central view displays:
- Biconcave RBC mesh with Phong shading
- Spectrin network overlay (toggle with 'S')
- Real-time membrane deformation under physics

---

## CLI Diagnostic Modes

Cell Simulator X provides extensive CLI diagnostics for validation and exploration.

### Physics Diagnostics

```bash
cargo run --release -- --diagnose
cargo run --release -- --diagnose -n 5000 -f 50.0
```

**Options:**
- `-n, --steps <N>` - Number of physics steps (default: 1000)
- `-f, --force <F>` - Applied force in μN (default: 5.0)

**Output:**
- Vertex displacement over time
- Velocity and acceleration
- Elastic energy

### Metabolism Diagnostics

```bash
cargo run --release -- --diagnose-metabolism
cargo run --release -- --diagnose-metabolism -d 60
cargo run --release -- --diagnose-metabolism --glucose-step 10.0
```

**Options:**
- `-d, --duration <SEC>` - Simulation duration in seconds (default: 60)
- `--glucose-step <MM>` - Apply glucose perturbation (mM)

**Output:**
- ATP, ADP concentrations
- 2,3-DPG levels
- Glucose consumption rate
- Glycolysis reaction rates

### Oxygen Transport Diagnostics

```bash
cargo run --release -- --diagnose-oxygen
cargo run --release -- --diagnose-oxygen --ph 7.2 --dpg 5.0 --temp 37
```

**Options:**
- `--ph <PH>` - Intracellular pH (default: 7.4)
- `--dpg <MM>` - 2,3-DPG concentration in mM (default: 5.0)
- `--temp <C>` - Temperature in °C (default: 37)

**Output:**
- Oxygen equilibrium curve (OEC)
- P50 (half-saturation pressure)
- Hill coefficient
- Bohr effect demonstration

### Integrated Metabolism-Oxygen

```bash
cargo run --release -- --diagnose-integrated
cargo run --release -- --diagnose-integrated --po2 40 -d 120
cargo run --release -- --diagnose-integrated --stress 5.0
```

**Options:**
- `--po2 <MMHG>` - Oxygen partial pressure (default: 100)
- `--stress <MULT>` - ATP consumption multiplier (default: 1.0)
- `-d, --duration <SEC>` - Duration in seconds

**Output:**
- Time series of lactate, pH, P50
- Coupling validation (lactate → pH → P50)
- Hemoglobin saturation

### Full Integration (Glycolysis + PPP + Redox)

```bash
cargo run --release -- --diagnose-full
cargo run --release -- --diagnose-full -d 120 --oxidative-stress 3.0
cargo run --release -- --diagnose-full --tension 2.0
```

**Options:**
- `-d, --duration <SEC>` - Duration (default: 120)
- `--oxidative-stress <MULT>` - H₂O₂ production multiplier (default: 1.0)
- `--tension <PN/NM>` - Membrane tension (default: 0.0)
- `--po2 <MMHG>` - Oxygen partial pressure (default: 100)

**Output:**
- 38-metabolite integration
- NADPH/NADP⁺ ratio
- GSH/GSSG ratio
- H₂O₂ steady state
- Ion concentrations (Na⁺, K⁺)
- Performance metrics

### Mechano-Metabolic Coupling

```bash
cargo run --release -- --diagnose-coupled
cargo run --release -- --diagnose-coupled --tension 2.0 -d 60
cargo run --release -- --diagnose-coupled --physics-substeps 500
```

**Options:**
- `--tension <PN/NM>` - Override membrane tension
- `-d, --duration <SEC>` - Duration in seconds
- `--physics-substeps <N>` - Physics steps per biochemistry step

**Output:**
- Membrane tension over time
- Piezo1 open probability
- Cytosolic Ca²⁺
- Spectrin stiffness modifier
- Bidirectional coupling validation

---

## Disease Modeling

### Storage Lesion (Blood Storage Aging)

```bash
# Fresh blood (day 0)
cargo run --release -- --diagnose-disease storage --disease-param 0

# Half-life (day 21)
cargo run --release -- --diagnose-disease storage --disease-param 21

# Expiration (day 42)
cargo run --release -- --diagnose-disease storage --disease-param 42
```

**Parameter**: Storage days (0-42)

**Effects:**
- ATP exponential decay (half-life 21 days)
- 2,3-DPG depletion (complete by day 14)
- Reduced pump efficiency
- Increased membrane leak

### Diabetic RBC

```bash
# Normal glucose
cargo run --release -- --diagnose-disease diabetic --disease-param 5

# Moderate hyperglycemia
cargo run --release -- --diagnose-disease diabetic --disease-param 12

# Severe hyperglycemia
cargo run --release -- --diagnose-disease diabetic --disease-param 15
```

**Parameter**: External glucose concentration (mM)

**Effects:**
- Elevated intracellular glucose
- 1.5× oxidative stress
- HbA1c accumulation

### Malaria (P. falciparum)

```bash
# Mild (1% parasitemia)
cargo run --release -- --diagnose-disease malaria --disease-param 0.01

# Moderate (5% parasitemia)
cargo run --release -- --diagnose-disease malaria --disease-param 0.05

# Severe (10% parasitemia)
cargo run --release -- --diagnose-disease malaria --disease-param 0.10
```

**Parameter**: Parasitemia fraction (0.0-1.0)

**Effects:**
- Glucose competition (parasite consumption)
- Elevated lactate
- 2× oxidative stress
- Stage-dependent metabolism

### Sickle Cell Disease

```bash
# Normal (HbAA)
cargo run --release -- --diagnose-disease sickle --disease-param 0.0

# Trait carrier (HbAS)
cargo run --release -- --diagnose-disease sickle --disease-param 0.4

# Disease (HbSS) with low O2
cargo run --release -- --diagnose-disease sickle --disease-param 1.0 --po2 40
```

**Parameter**: HbS fraction (0.0 = HbAA, 0.4 = HbAS, 1.0 = HbSS)

**Effects:**
- P50 right shift (26.8 → 31 mmHg for HbSS)
- Polymerization below 35% O₂ saturation
- 1.8× chronic oxidative stress

---

## Parameter Configuration

### Editing Parameters

All biological parameters are in `data/parameters/*.json`. To modify:

1. Open the relevant JSON file
2. Change parameter values
3. Rebuild: `cargo build --release`
4. Run diagnostics to validate

### Example: Increasing G6PDH Activity

```json
// data/parameters/pentose_phosphate.json
{
    "g6p_dehydrogenase": {
        "vmax_mM_per_sec": 0.08,  // Increased from 0.06
        ...
    }
}
```

### Creating Custom Configurations

For experiments, copy parameter files and modify:

```bash
cp data/parameters/glycolysis.json data/parameters/glycolysis_experiment.json
# Edit glycolysis_experiment.json
# Update code to load custom file
```

---

## Interpreting Output

### Understanding Diagnostic Output

#### Validation Status

```
✓ All parameters within physiological range
```
All metrics meet biological targets.

```
⚠️  ATP outside target: 1.2 mM (target: 1.5-2.5 mM)
```
Warning: metric outside expected range.

#### Time Series

```
  Time(s)      ATP      G6P  NADPH/NADP+  GSH/GSSG   H2O2(uM)       pH
--------------------------------------------------------------------------------
    0.00    2.000   0.0500        15.0     250.0       1.00    7.210
   12.00    1.870   0.0420        12.3     245.0       0.85    7.205
```

- **ATP**: Should stabilize at 1.5-2.5 mM
- **G6P**: Should remain >0.01 mM for PPP
- **NADPH/NADP⁺**: Target 10-20 at steady state
- **GSH/GSSG**: Target >100 (healthy), <50 indicates stress
- **H₂O₂**: Should be <5 μM
- **pH**: Should be ~7.2

### Key Coupling Validation

```
✓ G6P sustained at 0.039 mM (glycolysis → PPP coupling works)
✓ NADPH/NADP+ ratio sustained at 12.4 (PPP producing NADPH)
✓ GSH/GSSG ratio maintained at 245 (NADPH → GSH regeneration works)
```

These messages indicate successful metabolic coupling.

### Performance Metrics

```
Wall clock time: 1.23s
Simulation time: 120.00 s
Steps: 120000
Steps/second: 97561
```

- **Steps/second**: Higher is better. ~100K steps/s is typical.

---

## Common Workflows

### Workflow 1: Validating Baseline Metabolism

```bash
# Run full integration for 2 minutes
cargo run --release -- --diagnose-full -d 120

# Check key outputs:
# - ATP: 1.5-2.5 mM
# - NADPH/NADP+: 10-20
# - GSH/GSSG: >100
# - H2O2: <5 µM
# - Na+: 5-15 mM
# - K+: 140-150 mM
```

### Workflow 2: Testing Oxidative Stress Response

```bash
# Baseline
cargo run --release -- --diagnose-full -d 60

# 3× oxidative stress
cargo run --release -- --diagnose-full -d 60 --oxidative-stress 3.0

# 5× oxidative stress
cargo run --release -- --diagnose-full -d 60 --oxidative-stress 5.0

# Compare GSH/GSSG ratios and H2O2 levels
```

### Workflow 3: Blood Storage Simulation

```bash
# Simulate storage over time
for day in 0 7 14 21 28 35 42; do
    echo "=== Day $day ==="
    cargo run --release -- --diagnose-disease storage --disease-param $day -d 30
done

# Track ATP and 2,3-DPG depletion
```

### Workflow 4: Hypoxia Response

```bash
# Arterial conditions (high O2)
cargo run --release -- --diagnose-integrated --po2 100 -d 60

# Venous conditions (low O2)
cargo run --release -- --diagnose-integrated --po2 40 -d 60

# Severe hypoxia
cargo run --release -- --diagnose-integrated --po2 20 -d 60

# Compare O2 saturation and metabolic response
```

### Workflow 5: Mechanical Stimulation

```bash
# No tension (resting)
cargo run --release -- --diagnose-coupled -d 30

# Moderate tension
cargo run --release -- --diagnose-coupled --tension 1.5 -d 30

# High tension
cargo run --release -- --diagnose-coupled --tension 3.0 -d 30

# Compare Ca2+ levels and Piezo1 activation
```

---

## Troubleshooting

### Simulation Instabilities

**Symptom**: NaN or Inf values in output

**Solutions**:
- Reduce timestep (modify config)
- Check initial concentrations are positive
- Verify parameter values are reasonable

### ATP Crashes to Zero

**Symptom**: ATP drops below 1.0 mM rapidly

**Solutions**:
- Check ATP consumption rate isn't too high
- Verify glucose transport is enabled
- Check glycolysis enzyme Vmax values

### GSH Depletes Under Normal Conditions

**Symptom**: GSH/GSSG ratio <50 without stress

**Solutions**:
- Verify NADPH production (check PPP)
- Check glutathione reductase parameters
- Ensure H₂O₂ production rate is reasonable

### Performance Issues

**Symptom**: Slow simulation (<10K steps/sec)

**Solutions**:
- Use release build: `cargo run --release`
- Reduce simulation duration
- Close other applications

---

## Getting Help

- **Documentation**: See `docs/` directory
- **Parameters**: See `docs/PARAMETERS.md`
- **API Reference**: See `docs/API.md`
- **Issues**: Report bugs on GitHub
