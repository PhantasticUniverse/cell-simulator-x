//! Waugh & Evans 1979: micropipette aspiration shear modulus.
//!
//! Reference: Waugh R, Evans EA. Thermoelasticity of red blood cell
//! membrane. Biophys J. 1979;26(1):115-131.
//! doi:10.1016/S0006-3495(79)85239-X
//!
//! Earlier reference: Evans EA, Waugh R, Melnik L. Elastic area
//! compressibility modulus of red cell membrane. Biophys J.
//! 1976;16(6):585-595. doi:10.1016/S0006-3495(76)85713-X
//!
//! These two papers establish the canonical RBC membrane material constants:
//!
//!   - Shear modulus μ_s ≈ 6 × 10⁻³ dyn/cm = 6 μN/m at 25 °C.
//!     Temperature dependence: ~ -0.6%/K above 25 °C → ~5.5 μN/m at 37 °C.
//!   - Area expansion modulus K ≈ 450 mN/m = 450 000 μN/m.
//!   - Bending modulus k_c ≈ 0.18 pN·μm (Evans 1983).
//!
//! This experiment validates the model's `SkalakMaterial` defaults against
//! these published values directly.

use crate::physics::membrane::SkalakMaterial;
use crate::validation::{
    metrics::compute_metrics,
    reference_curve::{CurveMetadata, ValidationCurve},
    ExperimentResult,
};

const CHI2_THRESHOLD: f64 = 2.0;

fn material_constants_curve() -> ValidationCurve {
    // Three "data points": shear modulus, area modulus, bending modulus.
    // x = constant index (0, 1, 2); y = published value; σ = published uncertainty.
    let metadata = CurveMetadata {
        id: "waugh_evans_1979_material_constants".into(),
        citation: "Waugh & Evans 1979 Biophys J 26:115; Evans 1983 Biophys J 43:27".into(),
        doi: "10.1016/S0006-3495(79)85239-X".into(),
        figure: "Waugh-Evans 1979 Table 1; Evans-Waugh-Melnik 1976; Evans 1983".into(),
        digitization: "tabulated values; σ from quoted experimental scatter".into(),
        x_label: "constant index (0=μ_s, 1=K_a, 2=κ_c)".into(),
        y_label: "value (μN/m for μ_s,K_a; pN·μm for κ_c)".into(),
        notes: "Index 0: shear modulus 5.5 ± 1.1 μN/m (at 37°C). \
                Index 1: area expansion modulus 4.5e5 ± 1.0e5 μN/m. \
                Index 2: bending modulus 0.18 ± 0.05 pN·μm.".into(),
    };

    // Note: y values mix units. The metric is meaningful per-point because
    // sigma scales with the value; χ²/dof normalizes correctly.
    let x = vec![0.0, 1.0, 2.0];
    let y = vec![5.5, 450_000.0, 0.18];
    let sigma = vec![1.1, 100_000.0, 0.05];

    ValidationCurve::new(metadata, x, y, sigma)
}

/// Validation: RBC membrane material constants vs Waugh & Evans.
pub fn run_shear_modulus() -> ExperimentResult {
    let curve = material_constants_curve();
    let material = SkalakMaterial::default();

    let predicted = vec![
        material.shear_modulus_uN_per_m as f64,
        material.area_modulus_uN_per_m as f64,
        material.bending_modulus_pN_um as f64,
    ];

    let metrics = compute_metrics(&curve, &predicted, 0);
    let passed = metrics.chi2_per_dof < CHI2_THRESHOLD;

    let mut notes = Vec::new();
    if !passed {
        for (i, (&obs, &pred)) in curve.y.iter().zip(predicted.iter()).enumerate() {
            let z = (obs - pred) / curve.sigma[i].max(1e-12);
            if z.abs() > 1.5 {
                let label = match i {
                    0 => "shear modulus",
                    1 => "area modulus",
                    2 => "bending modulus",
                    _ => "?",
                };
                notes.push(format!("{} off by {:.2}σ (obs={}, model={})", label, z, obs, pred));
            }
        }
    }

    ExperimentResult {
        name: curve.metadata.id.clone(),
        citation: curve.metadata.citation.clone(),
        description: "RBC membrane material constants (μ_s, K_a, κ_c)".into(),
        metrics,
        passed,
        chi2_dof_threshold: CHI2_THRESHOLD,
        notes,
    }
}
