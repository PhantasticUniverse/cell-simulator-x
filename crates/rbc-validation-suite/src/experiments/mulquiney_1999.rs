//! Mulquiney & Kuchel 1999: steady-state RBC metabolite concentrations.
//!
//! Reference: Mulquiney PJ, Kuchel PW. Model of 2,3-bisphosphoglycerate
//! metabolism in the human erythrocyte based on detailed enzyme kinetic
//! equations: computer simulation and metabolic control analysis.
//! Biochem J. 1999;342(Pt 3):567-580. doi:10.1042/0264-6021:3420567
//!
//! This validates two distinct quantities:
//!
//! 1. Steady-state concentrations of glycolytic intermediates and cofactors
//!    (Mulquiney 1999 Table 4 / 5).
//! 2. PPP flux as a fraction of glycolytic flux (the canary the audit flagged
//!    as ~40% in the current model vs ~5–15% physiological).

use cell_simulator_x::biochemistry::{
    FullyIntegratedConfig, FullyIntegratedSolver, MetabolitePool, FullyIntegratedIndices,
};
use crate::{
    metrics::compute_metrics,
    reference_curve::{CurveMetadata, ValidationCurve},
    ExperimentResult,
};

const CHI2_THRESHOLD: f64 = 2.0;
/// Simulation duration to settle to steady state (seconds).
const STEADY_STATE_DURATION_SEC: f64 = 60.0;

fn steady_state_curve() -> ValidationCurve {
    // Steady-state values for human RBC at standard conditions.
    // Reference: Mulquiney & Kuchel 1999 Table 5; cross-checked against
    // Joshi & Palsson 1989 and Beutler 1984.
    //
    // Each "x" is a metabolite index (0..N) — categorical, but expressed
    // numerically so the existing curve infrastructure works. Order is
    // preserved by the implementation, not by x value.
    //
    // σ values: Mulquiney quotes ±20% on most metabolites due to
    // experimental scatter across samples (Beutler 1984). For low-copy
    // metabolites we use absolute σ = 20% × value, floor 1e-4 mM.
    let metabolite_targets: &[(&str, f64)] = &[
        ("Glucose",      4.5),
        ("G6P",          0.039),
        ("F6P",          0.013),
        ("F1,6BP",       0.007),
        ("DHAP",         0.143),
        ("G3P",          0.005),
        ("3PG",          0.061),
        ("2PG",          0.011),
        ("PEP",          0.017),
        ("Pyruvate",     0.085),
        ("Lactate",      1.5),
        ("ATP",          1.85),
        ("ADP",          0.14),
        ("2,3-BPG",      4.7),
        ("NAD+",         0.082),
        ("NADH",         0.012),
        ("NADPH",        0.034),
        ("NADP+",        0.002),
        ("GSH",          2.5),
        ("GSSG",         0.012),
    ];

    let metadata = CurveMetadata {
        id: "mulquiney_1999_steady_state".into(),
        citation: "Mulquiney PJ, Kuchel PW. Biochem J. 1999;342:567-580".into(),
        doi: "10.1042/0264-6021:3420567".into(),
        figure: "Table 5 (canonical steady-state); cross-checked Beutler 1984; Veech 1969 for NADPH".into(),
        digitization: "tabulated values from cited papers".into(),
        x_label: "metabolite index".into(),
        y_label: "[concentration] [mM]".into(),
        notes: format!(
            "RBC at standard conditions (pH 7.2 cytosolic, 37°C, 5 mM glucose). \
             σ = 20% × target value (Beutler scatter), floor 1e-4 mM. \
             Order: {}",
            metabolite_targets.iter().map(|(n, _)| *n).collect::<Vec<_>>().join(", ")
        ),
    };

    let x: Vec<f64> = (0..metabolite_targets.len()).map(|i| i as f64).collect();
    let y: Vec<f64> = metabolite_targets.iter().map(|(_, v)| *v).collect();
    let sigma: Vec<f64> = metabolite_targets.iter()
        .map(|(_, v)| (0.20 * v).max(1e-4))
        .collect();

    ValidationCurve::new(metadata, x, y, sigma)
}

fn ppp_flux_curve() -> ValidationCurve {
    // Single-point validation: PPP flux as fraction of glycolytic flux.
    // Reference: Beutler 1984 (5–15% physiological); Wood 1985 reviews;
    // Kirkman & Gaetani 1984 reports ~5–10% basal in normal RBCs.
    //
    // We use a single-point curve with target 0.10 (10%) ± 0.05 (1σ window
    // covers the full 5–15% range).
    let metadata = CurveMetadata {
        id: "mulquiney_1999_ppp_flux_fraction".into(),
        citation: "Beutler 1984; Kirkman & Gaetani 1984; Wood 1985".into(),
        doi: "".into(),
        figure: "Beutler 1984 § PPP".into(),
        digitization: "tabulated review value (range 5–15%); midpoint 10%, σ=5%".into(),
        x_label: "(constant)".into(),
        y_label: "PPP flux / glycolytic flux".into(),
        notes: "Critical canary: current model exhibits ~40% (audit-flagged structural compromise).".into(),
    };

    ValidationCurve::new(metadata, vec![0.0], vec![0.10], vec![0.05])
}

fn run_to_steady_state(duration_sec: f64) -> (FullyIntegratedSolver, MetabolitePool) {
    let config = FullyIntegratedConfig::default();
    let mut solver = FullyIntegratedSolver::new(config);
    let mut pool = MetabolitePool::default_fully_integrated();
    // Default basal ATP consumption (config.basal_atp_consumption_mM_per_sec);
    // pass 0.0 here for *additional* external load — the basal load is
    // already applied internally.
    solver.run(&mut pool, duration_sec, 0.0);
    (solver, pool)
}

/// Validation: steady-state metabolite concentrations.
pub fn run_steady_state_metabolites() -> ExperimentResult {
    let curve = steady_state_curve();
    let (_solver, pool) = run_to_steady_state(STEADY_STATE_DURATION_SEC);
    let idx = FullyIntegratedIndices::new();

    let predicted: Vec<f64> = vec![
        pool.get(idx.glycolysis.glucose),
        pool.get(idx.glycolysis.glucose_6_phosphate),
        pool.get(idx.glycolysis.fructose_6_phosphate),
        pool.get(idx.glycolysis.fructose_1_6_bisphosphate),
        pool.get(idx.glycolysis.dihydroxyacetone_phosphate),
        pool.get(idx.glycolysis.glyceraldehyde_3_phosphate),
        pool.get(idx.glycolysis.phosphoglycerate_3),
        pool.get(idx.glycolysis.phosphoglycerate_2),
        pool.get(idx.glycolysis.phosphoenolpyruvate),
        pool.get(idx.glycolysis.pyruvate),
        pool.get(idx.glycolysis.lactate),
        pool.get(idx.glycolysis.atp),
        pool.get(idx.glycolysis.adp),
        pool.get(idx.bisphosphoglycerate_2_3),
        pool.get(idx.glycolysis.nad),
        pool.get(idx.glycolysis.nadh),
        pool.get(idx.redox.nadph),
        pool.get(idx.redox.nadp_plus),
        pool.get(idx.redox.gsh),
        pool.get(idx.redox.gssg),
    ];

    let metrics = compute_metrics(&curve, &predicted, 0);
    let passed = metrics.chi2_per_dof < CHI2_THRESHOLD;

    let mut notes = Vec::new();
    // Per-metabolite diagnostic for the worst offender.
    if !passed {
        let worst = curve.y.iter().zip(predicted.iter()).enumerate()
            .map(|(i, (obs, pred))| {
                let r = (obs - pred) / curve.sigma[i].max(1e-12);
                (i, r * r)
            })
            .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
            .unwrap();
        notes.push(format!(
            "worst-fit metabolite at index {}: target {:.4}, model {:.4}",
            worst.0, curve.y[worst.0], predicted[worst.0]
        ));
    }

    ExperimentResult {
        name: curve.metadata.id.clone(),
        citation: curve.metadata.citation.clone(),
        description: "Steady-state metabolite concentrations after 60 s settling".into(),
        metrics,
        passed,
        chi2_dof_threshold: CHI2_THRESHOLD,
        notes,
    }
}

/// Validation: PPP flux as fraction of glycolytic flux (the audit canary).
pub fn run_ppp_flux_fraction() -> ExperimentResult {
    let curve = ppp_flux_curve();
    let (solver, pool) = run_to_steady_state(STEADY_STATE_DURATION_SEC);

    let diag = solver.diagnostics(&pool);
    // PPP fraction = G6PDH flux / HK flux. Already exposed by diagnostics.
    let predicted = vec![diag.ppp_fraction];
    let ppp_fraction = diag.ppp_fraction;
    let metrics = compute_metrics(&curve, &predicted, 0);
    let passed = metrics.chi2_per_dof < CHI2_THRESHOLD;

    let mut notes = Vec::new();
    notes.push(format!(
        "model PPP/glycolysis flux ratio: {:.1}% (target 5–15%)",
        100.0 * ppp_fraction
    ));
    if ppp_fraction > 0.20 {
        notes.push("STRUCTURAL: PPP gain exceeds physiological range — refit required".into());
    }

    ExperimentResult {
        name: curve.metadata.id.clone(),
        citation: curve.metadata.citation.clone(),
        description: "PPP flux fraction (audit canary)".into(),
        metrics,
        passed,
        chi2_dof_threshold: CHI2_THRESHOLD,
        notes,
    }
}
