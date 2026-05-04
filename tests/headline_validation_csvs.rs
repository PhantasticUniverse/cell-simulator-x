//! Phase 17.4: emit headline validation CSVs (figures 3 and 4).
//!
//! These tests call the prediction-only entry points exposed by
//! `rbc_validation_suite::experiments::*` directly — they do **not** run
//! the full validation suite (which would emit ~7 reference curves and
//! consume far more wall-clock than is needed for the figures).
//!
//! - Figure 3 (`tank_treading_emits_csv`): Keller-Skalak 1982 analytic
//!   prediction of tank-treading frequency vs shear rate, with the
//!   Fischer 2007 K(λ) ∈ [0.04, 0.15] band as the published reference.
//! - Figure 4 (`parachute_emits_csv`): Skalak 1973 capillary-number-
//!   driven parachute aspect ratio, with the reported AR ∈ [1.5, 2.0]
//!   physiological band.

use rbc_validation_suite::experiments::fischer_2007::keller_skalak_frequency_hz;
use rbc_validation_suite::experiments::skalak_1973::{
    capillary_number, empirical_aspect_ratio,
};
use std::io::Write;

/// Canonical RBC equatorial semi-axis (μm).
const ALPHA_UM: f64 = 4.0;
/// Canonical RBC polar semi-axis (μm).
const BETA_UM: f64 = 1.0;

/// Fischer 2007 K(λ) range at physiological → high external viscosity.
/// Used to bracket the Keller-Skalak prediction with a published band.
const FISCHER_K_LOW: f64 = 0.04;
const FISCHER_K_HIGH: f64 = 0.15;

/// Skalak 1973 reported parachute aspect-ratio band under physiological
/// capillary flow (Skalak & Branemark 1969, Fung 1993).
const SKALAK_AR_LOW: f64 = 1.5;
const SKALAK_AR_HIGH: f64 = 2.0;

/// Plasma viscosity for capillary number computation (Pa·s).
const PLASMA_VISCOSITY_PA_S: f64 = 1.2e-3;
/// Cell radius for capillary number (μm).
const CELL_RADIUS_UM: f64 = 4.0;
/// Membrane shear modulus for capillary number (μN/m).
const SHEAR_MODULUS_UN_PER_M: f64 = 5.0;

#[test]
fn tank_treading_emits_csv() {
    let target_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("target");
    std::fs::create_dir_all(&target_dir).expect("create target dir");
    let path = target_dir.join("tank_treading_fischer_2007.csv");
    let mut f = std::fs::File::create(&path).expect("create csv");
    writeln!(
        f,
        "shear_rate_per_sec,predicted_freq_hz,fischer_low_hz,fischer_high_hz"
    )
    .expect("write header");

    let shear_rates = [50.0, 100.0, 200.0, 350.0, 500.0, 700.0, 1000.0];
    for g in shear_rates {
        let predicted = keller_skalak_frequency_hz(g, ALPHA_UM, BETA_UM);
        // Fischer 2007: f_TT = K(λ)·γ̇. Bracket with K_low / K_high.
        let low = FISCHER_K_LOW * g;
        let high = FISCHER_K_HIGH * g;
        writeln!(
            f,
            "{:.4},{:.6},{:.6},{:.6}",
            g, predicted, low, high
        )
        .expect("write row");
    }
    drop(f);

    let written = std::fs::read_to_string(&path).expect("read back");
    assert!(written.starts_with("shear_rate_per_sec"));
    // 1 header + 7 shear rates.
    assert_eq!(written.lines().count(), 1 + shear_rates.len());
    println!(
        "Phase 17.4: wrote {} rows to {}",
        shear_rates.len(),
        path.display()
    );
}

#[test]
fn parachute_emits_csv() {
    let target_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("target");
    std::fs::create_dir_all(&target_dir).expect("create target dir");
    let path = target_dir.join("parachute_skalak_1973.csv");
    let mut f = std::fs::File::create(&path).expect("create csv");
    writeln!(
        f,
        "shear_rate_per_sec,capillary_number,predicted_aspect_ratio,skalak_low,skalak_high"
    )
    .expect("write header");

    let shear_rates = [100.0, 200.0, 400.0, 600.0, 800.0, 1000.0];
    for g in shear_rates {
        let ca = capillary_number(
            PLASMA_VISCOSITY_PA_S,
            g,
            CELL_RADIUS_UM,
            SHEAR_MODULUS_UN_PER_M,
        );
        let ar = empirical_aspect_ratio(ca);
        writeln!(
            f,
            "{:.4},{:.6},{:.6},{:.4},{:.4}",
            g, ca, ar, SKALAK_AR_LOW, SKALAK_AR_HIGH
        )
        .expect("write row");
    }
    drop(f);

    let written = std::fs::read_to_string(&path).expect("read back");
    assert!(written.starts_with("shear_rate_per_sec"));
    // 1 header + 6 shear rates.
    assert_eq!(written.lines().count(), 1 + shear_rates.len());
    println!(
        "Phase 17.4: wrote {} rows to {}",
        shear_rates.len(),
        path.display()
    );
}
