//! Validation tests for oxygen transport (Phase 4)
//!
//! These tests validate the hemoglobin oxygen binding model against
//! literature values from Imai 1982 and other sources.
//!
//! Key validation targets:
//! - P50 = 26.8 ± 1 mmHg at standard conditions
//! - Bohr coefficient = -0.48 ± 0.05
//! - Hill coefficient n = 2.7 ± 0.1
//! - 2,3-DPG sensitivity ~2.4 mmHg/mM

use cell_simulator_x::{
    HemoglobinSolver, HemoglobinState,
    STANDARD_PH, STANDARD_TEMPERATURE_K, STANDARD_DPG_MM, STANDARD_PCO2_MMHG,
};

// ============================================================================
// P50 Validation Tests
// ============================================================================

#[test]
fn test_p50_standard_conditions() {
    let solver = HemoglobinSolver::default();
    let p50 = solver.calculate_p50(
        STANDARD_PH,
        STANDARD_DPG_MM,
        STANDARD_TEMPERATURE_K,
        STANDARD_PCO2_MMHG,
    );

    // Target: 26.8 ± 1 mmHg (Imai 1982)
    let target = 26.8;
    let tolerance = 1.0;

    assert!(
        (p50 - target).abs() <= tolerance,
        "P50 at standard conditions: {:.1} mmHg (expected {:.1} ± {:.1})",
        p50, target, tolerance
    );
}

#[test]
fn test_p50_at_50_percent_saturation() {
    let solver = HemoglobinSolver::default();
    let p50 = solver.calculate_p50(
        STANDARD_PH,
        STANDARD_DPG_MM,
        STANDARD_TEMPERATURE_K,
        STANDARD_PCO2_MMHG,
    );

    // Verify that saturation at P50 is approximately 50%
    let sat_at_p50 = solver.calculate_saturation(
        p50,
        STANDARD_PH,
        STANDARD_DPG_MM,
        STANDARD_TEMPERATURE_K,
        STANDARD_PCO2_MMHG,
    );

    assert!(
        (sat_at_p50 - 0.5).abs() < 0.01,
        "Saturation at P50 should be ~50%: {:.1}%",
        sat_at_p50 * 100.0
    );
}

// ============================================================================
// Bohr Effect Tests
// ============================================================================

#[test]
fn test_bohr_coefficient() {
    let solver = HemoglobinSolver::default();
    let bohr = solver.measured_bohr_coefficient(
        STANDARD_DPG_MM,
        STANDARD_TEMPERATURE_K,
        STANDARD_PCO2_MMHG,
    );

    // Target: -0.48 ± 0.05 (Imai 1982)
    let target = -0.48;
    let tolerance = 0.05;

    assert!(
        (bohr - target).abs() <= tolerance,
        "Bohr coefficient: {:.2} (expected {:.2} ± {:.2})",
        bohr, target, tolerance
    );
}

#[test]
fn test_bohr_effect_direction() {
    let solver = HemoglobinSolver::default();

    // Lower pH should increase P50 (right shift, oxygen release)
    let p50_low_ph = solver.calculate_p50(7.2, STANDARD_DPG_MM, STANDARD_TEMPERATURE_K, STANDARD_PCO2_MMHG);
    let p50_high_ph = solver.calculate_p50(7.6, STANDARD_DPG_MM, STANDARD_TEMPERATURE_K, STANDARD_PCO2_MMHG);

    assert!(
        p50_low_ph > p50_high_ph,
        "Lower pH ({:.1} mmHg) should have higher P50 than higher pH ({:.1} mmHg)",
        p50_low_ph, p50_high_ph
    );
}

#[test]
fn test_saturation_at_different_ph() {
    let solver = HemoglobinSolver::default();
    let po2 = 40.0;  // Typical venous pO2

    // At the same pO2, lower pH should give lower saturation
    let sat_low_ph = solver.calculate_saturation(po2, 7.2, STANDARD_DPG_MM, STANDARD_TEMPERATURE_K, STANDARD_PCO2_MMHG);
    let sat_high_ph = solver.calculate_saturation(po2, 7.6, STANDARD_DPG_MM, STANDARD_TEMPERATURE_K, STANDARD_PCO2_MMHG);

    assert!(
        sat_low_ph < sat_high_ph,
        "At pO2 {}: saturation at pH 7.2 ({:.1}%) should be < pH 7.6 ({:.1}%)",
        po2, sat_low_ph * 100.0, sat_high_ph * 100.0
    );
}

// ============================================================================
// Hill Coefficient Tests
// ============================================================================

#[test]
fn test_hill_coefficient() {
    let solver = HemoglobinSolver::default();
    let n = solver.calculate_hill_coefficient(
        STANDARD_PH,
        STANDARD_DPG_MM,
        STANDARD_TEMPERATURE_K,
        STANDARD_PCO2_MMHG,
    );

    // Target: 2.7 ± 0.1 (Imai 1982)
    let target = 2.7;
    let tolerance = 0.1;

    assert!(
        (n - target).abs() <= tolerance,
        "Hill coefficient: {:.2} (expected {:.1} ± {:.1})",
        n, target, tolerance
    );
}

#[test]
fn test_cooperativity() {
    let solver = HemoglobinSolver::default();

    // The OEC should be sigmoidal, not hyperbolic
    // At P50, the slope of the Hill plot should be > 1
    let n = solver.calculate_hill_coefficient(
        STANDARD_PH,
        STANDARD_DPG_MM,
        STANDARD_TEMPERATURE_K,
        STANDARD_PCO2_MMHG,
    );

    assert!(
        n > 2.0,
        "Hill coefficient should indicate cooperativity (n > 2): {:.2}",
        n
    );
}

// ============================================================================
// 2,3-DPG Effect Tests
// ============================================================================

#[test]
fn test_dpg_effect_direction() {
    let solver = HemoglobinSolver::default();

    // Higher 2,3-DPG should increase P50 (right shift)
    let p50_low_dpg = solver.calculate_p50(STANDARD_PH, 3.0, STANDARD_TEMPERATURE_K, STANDARD_PCO2_MMHG);
    let p50_high_dpg = solver.calculate_p50(STANDARD_PH, 7.0, STANDARD_TEMPERATURE_K, STANDARD_PCO2_MMHG);

    assert!(
        p50_high_dpg > p50_low_dpg,
        "Higher 2,3-DPG ({:.1} mmHg) should have higher P50 than lower 2,3-DPG ({:.1} mmHg)",
        p50_high_dpg, p50_low_dpg
    );
}

#[test]
fn test_dpg_sensitivity() {
    let solver = HemoglobinSolver::default();

    // Calculate sensitivity: ΔP50 / Δ[2,3-DPG]
    let p50_low_dpg = solver.calculate_p50(STANDARD_PH, 3.0, STANDARD_TEMPERATURE_K, STANDARD_PCO2_MMHG);
    let p50_high_dpg = solver.calculate_p50(STANDARD_PH, 7.0, STANDARD_TEMPERATURE_K, STANDARD_PCO2_MMHG);

    let delta_p50 = p50_high_dpg - p50_low_dpg;
    let delta_dpg = 7.0 - 3.0;
    let sensitivity = delta_p50 / delta_dpg;

    // Target: ~2.4 mmHg/mM (Benesch 1969), allow 1.5-3.5 range
    assert!(
        sensitivity > 1.5 && sensitivity < 3.5,
        "2,3-DPG sensitivity: {:.1} mmHg/mM (expected ~2.4)",
        sensitivity
    );
}

// ============================================================================
// Temperature Effect Tests
// ============================================================================

#[test]
fn test_temperature_effect_direction() {
    let solver = HemoglobinSolver::default();

    // Higher temperature should increase P50 (right shift, less affinity)
    let p50_30C = solver.calculate_p50(STANDARD_PH, STANDARD_DPG_MM, 303.15, STANDARD_PCO2_MMHG);  // 30°C
    let p50_40C = solver.calculate_p50(STANDARD_PH, STANDARD_DPG_MM, 313.15, STANDARD_PCO2_MMHG);  // 40°C

    assert!(
        p50_40C > p50_30C,
        "Higher temp (40°C: {:.1} mmHg) should have higher P50 than lower temp (30°C: {:.1} mmHg)",
        p50_40C, p50_30C
    );
}

// ============================================================================
// Saturation Limits Tests
// ============================================================================

#[test]
fn test_saturation_at_zero_po2() {
    let solver = HemoglobinSolver::default();

    let sat = solver.calculate_saturation(
        0.0,
        STANDARD_PH,
        STANDARD_DPG_MM,
        STANDARD_TEMPERATURE_K,
        STANDARD_PCO2_MMHG,
    );

    assert_eq!(sat, 0.0, "Saturation at pO2 = 0 should be 0");
}

#[test]
fn test_saturation_at_high_po2() {
    let solver = HemoglobinSolver::default();

    let sat = solver.calculate_saturation(
        500.0,  // Very high pO2
        STANDARD_PH,
        STANDARD_DPG_MM,
        STANDARD_TEMPERATURE_K,
        STANDARD_PCO2_MMHG,
    );

    assert!(
        sat > 0.99,
        "Saturation at pO2 = 500 mmHg should be >99%: {:.1}%",
        sat * 100.0
    );
}

#[test]
fn test_saturation_monotonic() {
    let solver = HemoglobinSolver::default();

    let mut prev_sat = 0.0;
    for po2 in [5.0, 10.0, 20.0, 30.0, 50.0, 80.0, 120.0] {
        let sat = solver.calculate_saturation(
            po2,
            STANDARD_PH,
            STANDARD_DPG_MM,
            STANDARD_TEMPERATURE_K,
            STANDARD_PCO2_MMHG,
        );
        assert!(
            sat >= prev_sat,
            "Saturation should increase monotonically: {:.1}% at {} mmHg vs {:.1}% previously",
            sat * 100.0, po2, prev_sat * 100.0
        );
        prev_sat = sat;
    }
}

// ============================================================================
// Hemoglobin State Tests
// ============================================================================

#[test]
fn test_hemoglobin_state_creation() {
    let state = HemoglobinState::at_saturation(0.75, 5.0);

    assert!((state.saturation - 0.75).abs() < 1e-6);
    assert!((state.total_hb_mM - 5.0).abs() < 1e-6);
    assert!((state.max_oxygen_capacity_mM() - 20.0).abs() < 1e-6);  // 4 × 5.0
    assert!((state.bound_o2_mM - 15.0).abs() < 1e-6);  // 0.75 × 20.0
}

#[test]
fn test_hemoglobin_state_update() {
    let mut state = HemoglobinState::default();

    state.update_from_saturation(0.5);

    assert!((state.saturation - 0.5).abs() < 1e-6);
    assert!((state.bound_o2_mM - 10.0).abs() < 1e-6);  // 0.5 × 20.0
}

// ============================================================================
// Dynamic Binding Tests
// ============================================================================

#[test]
fn test_binding_approaches_equilibrium() {
    let solver = HemoglobinSolver::default();
    let mut state = HemoglobinState::at_saturation(0.5, 5.0);
    let conditions = (STANDARD_PH, STANDARD_DPG_MM, STANDARD_TEMPERATURE_K, STANDARD_PCO2_MMHG);

    // Apply high pO2 and verify saturation increases
    let initial_sat = state.saturation;
    let po2 = 100.0;  // High pO2

    for _ in 0..1000 {
        solver.step(&mut state, po2, conditions, 0.001);
    }

    // Should approach equilibrium saturation
    let eq_sat = solver.calculate_saturation(po2, STANDARD_PH, STANDARD_DPG_MM, STANDARD_TEMPERATURE_K, STANDARD_PCO2_MMHG);

    assert!(
        state.saturation > initial_sat,
        "Saturation should increase at high pO2: {:.1}% -> {:.1}%",
        initial_sat * 100.0, state.saturation * 100.0
    );

    // Should be approaching equilibrium
    assert!(
        (state.saturation - eq_sat).abs() < 0.05,
        "Should approach equilibrium: {:.1}% (equilibrium: {:.1}%)",
        state.saturation * 100.0, eq_sat * 100.0
    );
}

#[test]
fn test_oxygen_release_at_low_po2() {
    let solver = HemoglobinSolver::default();
    let mut state = HemoglobinState::at_saturation(0.97, 5.0);  // Arterial saturation
    let conditions = (STANDARD_PH, STANDARD_DPG_MM, STANDARD_TEMPERATURE_K, STANDARD_PCO2_MMHG);

    // Apply low pO2 (tissue level) and verify saturation decreases
    let initial_sat = state.saturation;
    let po2 = 20.0;  // Low pO2 (tissue)

    for _ in 0..1000 {
        solver.step(&mut state, po2, conditions, 0.001);
    }

    assert!(
        state.saturation < initial_sat,
        "Saturation should decrease at low pO2: {:.1}% -> {:.1}%",
        initial_sat * 100.0, state.saturation * 100.0
    );
}

// ============================================================================
// OEC Generation Tests
// ============================================================================

#[test]
fn test_oec_generation() {
    let solver = HemoglobinSolver::default();
    let conditions = (STANDARD_PH, STANDARD_DPG_MM, STANDARD_TEMPERATURE_K, STANDARD_PCO2_MMHG);

    let oec = solver.generate_oec(conditions, 50);

    assert_eq!(oec.len(), 50, "OEC should have 50 points");

    // First point should have low saturation
    assert!(oec[0].1 < 0.5, "First point saturation should be < 50%");

    // Last point should have high saturation
    assert!(oec[oec.len()-1].1 > 0.95, "Last point saturation should be > 95%");
}

// ============================================================================
// Adair Equation Tests
// ============================================================================

#[test]
fn test_adair_constants_cooperativity() {
    use cell_simulator_x::biochemistry::AdairConstants;

    let constants = AdairConstants::default();

    // Successive K values should increase (positive cooperativity)
    assert!(constants.k1_per_mmHg < constants.k2_per_mmHg);
    assert!(constants.k2_per_mmHg < constants.k3_per_mmHg);
    assert!(constants.k3_per_mmHg < constants.k4_per_mmHg);
}

#[test]
fn test_adair_from_p50_hill() {
    use cell_simulator_x::biochemistry::AdairConstants;

    let target_p50 = 30.0;
    let target_n = 2.5;

    let constants = AdairConstants::from_p50_and_hill(target_p50, target_n);

    // Verify cooperativity is preserved
    assert!(constants.k1_per_mmHg < constants.k4_per_mmHg,
        "K1 ({}) should be less than K4 ({})",
        constants.k1_per_mmHg, constants.k4_per_mmHg
    );
}

// ============================================================================
// Integration Test
// ============================================================================

#[test]
fn test_full_oxygen_transport_validation() {
    let solver = HemoglobinSolver::default();

    // Calculate all validation metrics
    let p50 = solver.calculate_p50(STANDARD_PH, STANDARD_DPG_MM, STANDARD_TEMPERATURE_K, STANDARD_PCO2_MMHG);
    let hill_n = solver.calculate_hill_coefficient(STANDARD_PH, STANDARD_DPG_MM, STANDARD_TEMPERATURE_K, STANDARD_PCO2_MMHG);
    let bohr_coeff = solver.measured_bohr_coefficient(STANDARD_DPG_MM, STANDARD_TEMPERATURE_K, STANDARD_PCO2_MMHG);

    // P50: 26.8 ± 1 mmHg
    assert!(
        (p50 - 26.8).abs() <= 1.0,
        "P50 validation failed: {:.1} (target: 26.8 ± 1.0)",
        p50
    );

    // Hill coefficient: 2.7 ± 0.1
    assert!(
        (hill_n - 2.7).abs() <= 0.1,
        "Hill coefficient validation failed: {:.2} (target: 2.7 ± 0.1)",
        hill_n
    );

    // Bohr coefficient: -0.48 ± 0.05
    assert!(
        (bohr_coeff - (-0.48)).abs() <= 0.05,
        "Bohr coefficient validation failed: {:.2} (target: -0.48 ± 0.05)",
        bohr_coeff
    );

    println!("=== Phase 4 Oxygen Transport Validation ===");
    println!("P50:              {:.1} mmHg (target: 26.8 ± 1.0) ✓", p50);
    println!("Hill coefficient: {:.2} (target: 2.7 ± 0.1) ✓", hill_n);
    println!("Bohr coefficient: {:.2} (target: -0.48 ± 0.05) ✓", bohr_coeff);
    println!("All validation targets met!");
}
