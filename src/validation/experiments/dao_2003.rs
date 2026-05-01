//! Dao 2003: optical-tweezers RBC stretching, axial vs transverse strain.
//!
//! Reference: Dao M, Lim CT, Suresh S. Mechanics of the human red blood
//! cell deformed by optical tweezers. J Mech Phys Solids. 2003;51:2259-2280.
//! doi:10.1016/j.jmps.2003.09.019
//!
//! At applied stretching forces between ~70 and ~200 pN, the RBC elongates
//! axially and contracts transversely. Dao et al. report axial diameter
//! gain (D_axial / D_0) and transverse decrease (D_trans / D_0) as
//! functions of applied force, fit to extract the membrane shear modulus.
//!
//! This experiment is the most demanding to reproduce because it requires
//! the full mechanics solver (DPD + Skalak + WLC) running to mechanical
//! equilibrium under prescribed force. For Phase 10 we validate against the
//! *force-vs-strain relationship* extracted by Dao's analytical fit, which
//! gives an effective shear modulus of 5.3 ± 1.0 μN/m at small strain.
//!
//! A full mechanical simulation requires Phase 12 (vessel geometry +
//! immersed boundary). Until then, we report a "configuration-only" check
//! that the model's shear modulus parameter falls in Dao's measured band.

use crate::physics::membrane::SkalakMaterial;
use crate::validation::{
    metrics::compute_metrics,
    reference_curve::{CurveMetadata, ValidationCurve},
    ExperimentResult,
};

const CHI2_THRESHOLD: f64 = 2.0;

fn shear_modulus_extraction_curve() -> ValidationCurve {
    // Dao 2003 reports μ_s = 5.3 ± 1.0 μN/m (small-strain, optical-tweezers
    // fit). Single-point validation against the published value with σ=1.0.
    let metadata = CurveMetadata {
        id: "dao_2003_extracted_shear_modulus".into(),
        citation: "Dao M, Lim CT, Suresh S. J Mech Phys Solids 2003;51:2259".into(),
        doi: "10.1016/j.jmps.2003.09.019".into(),
        figure: "Dao 2003 Table 2 (extracted shear modulus from optical tweezers)".into(),
        digitization: "tabulated value; σ from quoted experimental scatter across cells".into(),
        x_label: "(constant)".into(),
        y_label: "shear modulus μ_s [μN/m]".into(),
        notes: "Optical-tweezers small-strain extraction. PHASE 10: configuration-only check; \
                full force-extension validation requires Phase 12 vessel-geometry mechanics.".into(),
    };

    ValidationCurve::new(metadata, vec![0.0], vec![5.3], vec![1.0])
}

/// Validation: shear modulus matches Dao 2003 optical-tweezers extraction.
pub fn run_axial_extension_curve() -> ExperimentResult {
    let curve = shear_modulus_extraction_curve();
    let material = SkalakMaterial::default();

    let predicted = vec![material.shear_modulus_uN_per_m as f64];
    let metrics = compute_metrics(&curve, &predicted, 0);
    let passed = metrics.chi2_per_dof < CHI2_THRESHOLD;

    let mut notes = vec![
        format!("PHASE 10: only checks shear modulus parameter (model = {} μN/m).",
            material.shear_modulus_uN_per_m),
        "Full optical-tweezers axial/transverse strain reproduction deferred to Phase 12".into(),
    ];
    if !passed {
        notes.push("model shear modulus falls outside Dao 2003 measurement band".into());
    }

    ExperimentResult {
        name: curve.metadata.id.clone(),
        citation: curve.metadata.citation.clone(),
        description: "Dao 2003 optical-tweezers shear modulus (config-only)".into(),
        metrics,
        passed,
        chi2_dof_threshold: CHI2_THRESHOLD,
        notes,
    }
}
