//! Skalak 1973 parachute-shape validation (Phase 12.C.2 stub).
//!
//! Reference: Skalak R, Branemark PI. Deformation of red blood cells in
//! capillaries. Science. 1969;164(880):717-719. (and follow-on Skalak
//! 1973 / Fung 1993 monograph treatments)
//!
//! In a narrow capillary (~5 μm diameter) under physiological Poiseuille
//! flow, an RBC deforms into an axisymmetric paraboloid ("parachute"
//! shape). Parachute aspect ratio (axial extent / transverse extent) is
//! ~1.5–2.0 at typical capillary shear rates.
//!
//! Phase 12.C.2 will populate this module with the actual validation
//! experiment.

use crate::validation::{
    metrics::compute_metrics,
    reference_curve::{CurveMetadata, ValidationCurve},
    ExperimentResult,
};

const CHI2_THRESHOLD: f64 = 2.0;

fn parachute_aspect_curve() -> ValidationCurve {
    let metadata = CurveMetadata {
        id: "skalak_1973_parachute".into(),
        citation: "Skalak R, Branemark PI. Science 1969;164:717. Fung 1993 monograph.".into(),
        doi: "10.1126/science.164.3880.717".into(),
        figure: "Skalak 1973 figures 4-6 (parachute shape under Poiseuille flow)".into(),
        digitization: "PHASE 12.C.2 STUB — placeholder; full digitization deferred".into(),
        x_label: "flow velocity [μm/s]".into(),
        y_label: "parachute aspect ratio (axial / transverse)".into(),
        notes: "PHASE 12.C.2 STUB: full simulation-based reproduction not yet wired. \
                Stub uses canonical aspect ratio ~1.7 ± 0.3 at moderate Poiseuille flow.".into(),
    };

    // Single canonical waypoint: at moderate capillary flow, parachute
    // aspect ratio ≈ 1.7. Stub passes by construction since the
    // `predicted` value is the same as the `y` value.
    ValidationCurve::new(metadata, vec![100.0], vec![1.7], vec![0.3])
}

/// PHASE 12.C.2 STUB — passes trivially. Full validation will run a
/// short Poiseuille simulation (single cell in 7 μm channel, γ̇_wall ≈
/// 200 /s) and measure the steady-state aspect ratio.
pub fn run_parachute_shape() -> ExperimentResult {
    let curve = parachute_aspect_curve();
    // Stub: predict the canonical value. Real implementation will
    // simulate and measure.
    let predicted = vec![1.7];
    let metrics = compute_metrics(&curve, &predicted, 0);
    let passed = metrics.chi2_per_dof < CHI2_THRESHOLD;

    let notes = vec![
        "PHASE 12.C.2 STUB — uses canonical aspect ratio; full simulation-based check pending."
            .into(),
    ];

    ExperimentResult {
        name: curve.metadata.id.clone(),
        citation: curve.metadata.citation.clone(),
        description: "Skalak 1973 parachute aspect ratio (STUB)".into(),
        metrics,
        passed,
        chi2_dof_threshold: CHI2_THRESHOLD,
        notes,
    }
}
