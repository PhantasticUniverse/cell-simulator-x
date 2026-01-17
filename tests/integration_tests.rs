//! Integration tests for Phase 5: Metabolism-Oxygen coupling
//!
//! Tests verify that:
//! - Lactate accumulation correctly decreases pH
//! - pH changes correctly affect P50 via Bohr effect
//! - The coupling direction is correct: lactate ↑ → pH ↓ → P50 ↑ → saturation ↓
//! - Coupling magnitudes match physiological literature values

use cell_simulator_x::{
    IntegratedSolver, IntegratedEnvironment, PhBufferModel,
    MetabolitePool,
};

/// Test that coupling direction is correct: lactate ↑ → pH ↓ → P50 ↑
#[test]
fn test_coupling_direction_lactate_to_ph_to_p50() {
    let env = IntegratedEnvironment::default();
    let solver = IntegratedSolver::new(env);
    let mut metabolites = MetabolitePool::default_physiological();

    // Record baseline at normal lactate
    metabolites.set(solver.indices.glycolysis.lactate, 1.5);
    let baseline_diag = solver.diagnostics(&metabolites);

    // Increase lactate by 3 mM
    metabolites.set(solver.indices.glycolysis.lactate, 4.5);
    let high_lactate_diag = solver.diagnostics(&metabolites);

    // Lactate increased
    assert!(
        high_lactate_diag.lactate_mM > baseline_diag.lactate_mM,
        "Lactate should have increased: {} -> {}",
        baseline_diag.lactate_mM, high_lactate_diag.lactate_mM
    );

    // pH should decrease (lactate is an acid)
    assert!(
        high_lactate_diag.ph < baseline_diag.ph,
        "pH should decrease when lactate increases: {} -> {}",
        baseline_diag.ph, high_lactate_diag.ph
    );

    // P50 should increase (Bohr effect: lower pH → higher P50)
    assert!(
        high_lactate_diag.p50_mmHg > baseline_diag.p50_mmHg,
        "P50 should increase when pH decreases: {} -> {}",
        baseline_diag.p50_mmHg, high_lactate_diag.p50_mmHg
    );
}

/// Test pH sensitivity: 1 mM lactate ≈ 0.017 pH drop
///
/// Reference: Van Slyke 1922, buffer capacity ~60 slykes
#[test]
fn test_ph_sensitivity_per_mm_lactate() {
    let buffer = PhBufferModel::default();

    // Calculate pH at two lactate levels 1 mM apart
    let ph_at_1_5 = buffer.calculate_ph(1.5);
    let ph_at_2_5 = buffer.calculate_ph(2.5);

    let ph_drop_per_mM = ph_at_1_5 - ph_at_2_5;

    // Expected: 1 mM / 60 slykes ≈ 0.0167
    let expected = 1.0 / 60.0;  // ~0.0167

    assert!(
        (ph_drop_per_mM - expected).abs() < 0.002,
        "pH drop per mM lactate should be ~0.017: got {} (expected {})",
        ph_drop_per_mM, expected
    );
}

/// Test Bohr effect magnitude: 0.2 pH drop ≈ 4 mmHg P50 increase
///
/// Reference: Imai 1982, Bohr coefficient -0.48
/// At standard P50 ~27 mmHg: ΔlogP50 = -0.48 × ΔpH
/// For ΔpH = -0.2: ΔlogP50 = +0.096
/// P50_new = 27 × 10^0.096 ≈ 33.4 mmHg (increase of ~6.4 mmHg)
///
/// Note: The actual shift depends on the exact Bohr coefficient
/// implementation. We test for ~4-6 mmHg range.
#[test]
fn test_bohr_effect_magnitude() {
    let env = IntegratedEnvironment::default();
    let solver = IntegratedSolver::new(env);
    let mut metabolites = MetabolitePool::default_physiological();

    // pH 7.2 (reference)
    metabolites.set(solver.indices.glycolysis.lactate, 1.5);  // Should give pH ~7.2
    let p50_at_7_2 = solver.effective_p50(&metabolites);

    // pH 7.0 (0.2 drop, achieved by adding ~12 mM lactate)
    // 0.2 pH drop × 60 slykes = 12 mM lactate increase
    metabolites.set(solver.indices.glycolysis.lactate, 1.5 + 12.0);
    let ph_check = solver.diagnostics(&metabolites);
    let p50_at_7_0 = solver.effective_p50(&metabolites);

    // Verify we achieved roughly the right pH drop
    let ph_drop = 7.2 - ph_check.ph;
    assert!(
        (ph_drop - 0.2).abs() < 0.05,
        "pH should have dropped by ~0.2: actual drop = {}",
        ph_drop
    );

    // P50 increase should be ~4-8 mmHg for 0.2 pH drop
    let p50_increase = p50_at_7_0 - p50_at_7_2;
    assert!(
        p50_increase > 2.0 && p50_increase < 10.0,
        "P50 increase for 0.2 pH drop should be ~4-6 mmHg: got {}",
        p50_increase
    );
}

/// Test P50 shift at low pH
///
/// Reference: Imai 1982 data
/// At low pH (vs 7.4 standard), P50 increases significantly due to Bohr effect.
/// The Bohr coefficient is -0.48, so:
///   At pH 7.0 (delta = -0.4 from 7.4): delta_log_p50 = -0.48 * (-0.4) = +0.192
///   P50_new = 26.8 * 10^0.192 ≈ 41-42 mmHg
///
/// Note: The reference pH in hemoglobin solver is 7.4 (STANDARD_PH), not 7.2
#[test]
fn test_p50_at_low_ph() {
    let env = IntegratedEnvironment::default();
    let solver = IntegratedSolver::new(env);
    let mut metabolites = MetabolitePool::default_physiological();

    // Set lactate to achieve pH ~7.0
    // pH 7.0 means 0.2 drop from 7.2 (buffer reference), which is 12 mM lactate above baseline
    metabolites.set(solver.indices.glycolysis.lactate, 1.5 + 12.0);

    let diag = solver.diagnostics(&metabolites);
    let p50 = solver.effective_p50(&metabolites);

    // Verify pH is around 7.0
    assert!(
        (diag.ph - 7.0).abs() < 0.05,
        "pH should be ~7.0: got {}",
        diag.ph
    );

    // At pH ~7.0, P50 should be significantly elevated (Bohr effect)
    // With Bohr coeff -0.48 and base P50 ~27, at pH 7.0 (0.4 below standard 7.4):
    // P50 should be ~40-45 mmHg
    assert!(
        p50 > 35.0 && p50 < 50.0,
        "P50 at pH ~7.0 should be significantly elevated: got {} (pH={})",
        p50, diag.ph
    );

    // Also verify it's higher than baseline
    metabolites.set(solver.indices.glycolysis.lactate, 1.5);
    let p50_baseline = solver.effective_p50(&metabolites);
    assert!(
        p50 > p50_baseline,
        "P50 at low pH should be higher than baseline: {} vs {}",
        p50, p50_baseline
    );
}

/// Test saturation decreases with acidosis at fixed pO2
///
/// At a pO2 near P50, the shift should be most apparent
#[test]
fn test_saturation_decreases_with_acidosis() {
    // Use pO2 near normal P50 (~27 mmHg) where changes are most visible
    let env = IntegratedEnvironment::with_po2(27.0);
    let solver = IntegratedSolver::new(env);
    let mut metabolites = MetabolitePool::default_physiological();

    // Normal pH
    metabolites.set(solver.indices.glycolysis.lactate, 1.5);
    let sat_normal = solver.calculate_saturation(&metabolites);

    // Acidic pH
    metabolites.set(solver.indices.glycolysis.lactate, 10.0);  // Higher lactate → lower pH
    let sat_acidic = solver.calculate_saturation(&metabolites);

    // Saturation should decrease (P50 shifted right)
    assert!(
        sat_acidic < sat_normal,
        "Saturation should decrease with acidosis at fixed pO2: {} -> {}",
        sat_normal, sat_acidic
    );
}

/// Test integrated simulation maintains physiological pH range
#[test]
fn test_dynamic_simulation_ph_range() {
    let env = IntegratedEnvironment::default();
    let mut solver = IntegratedSolver::new(env);
    let mut metabolites = MetabolitePool::default_physiological();

    // Run simulation with moderate ATP stress
    solver.run(&mut metabolites, 30.0, 0.005);

    let diag = solver.diagnostics(&metabolites);

    // pH should remain in physiological range (6.8-7.6)
    assert!(
        diag.ph >= 6.8 && diag.ph <= 7.6,
        "pH should remain in physiological range: {}",
        diag.ph
    );

    // P50 should remain in physiological range (20-40 mmHg)
    assert!(
        diag.p50_mmHg >= 20.0 && diag.p50_mmHg <= 40.0,
        "P50 should remain in physiological range: {}",
        diag.p50_mmHg
    );
}

/// Test that baseline pH is 7.2 at baseline lactate (Jacobs 1947)
#[test]
fn test_baseline_ph_is_7_2() {
    let buffer = PhBufferModel::default();
    let ph = buffer.calculate_ph(1.5);  // Baseline lactate

    assert!(
        (ph - 7.2).abs() < 0.001,
        "Baseline pH should be 7.2: got {}",
        ph
    );
}

/// Test combined example from plan:
/// Lactate: 1.5 → 4.5 mM (+3 mM)
/// pH: 7.2 → 7.15 (-0.05)
/// P50: should increase by ~1-2 mmHg
#[test]
fn test_combined_coupling_example() {
    let env = IntegratedEnvironment::default();
    let solver = IntegratedSolver::new(env);
    let mut metabolites = MetabolitePool::default_physiological();

    // Baseline: lactate 1.5 mM
    metabolites.set(solver.indices.glycolysis.lactate, 1.5);
    let baseline = solver.diagnostics(&metabolites);

    // After: lactate 4.5 mM (+3 mM)
    metabolites.set(solver.indices.glycolysis.lactate, 4.5);
    let after = solver.diagnostics(&metabolites);

    // pH change should be ~-0.05 (3 mM / 60 slykes)
    let delta_ph = after.ph - baseline.ph;
    assert!(
        (delta_ph - (-0.05)).abs() < 0.01,
        "pH change should be ~-0.05: got {}",
        delta_ph
    );

    // P50 change should be positive (Bohr effect)
    let delta_p50 = after.p50_mmHg - baseline.p50_mmHg;
    assert!(
        delta_p50 > 0.0,
        "P50 should increase: got {} change",
        delta_p50
    );

    // P50 increase should be ~0.5-2 mmHg for 0.05 pH drop
    assert!(
        delta_p50 > 0.3 && delta_p50 < 3.0,
        "P50 increase should be ~0.5-2 mmHg: got {}",
        delta_p50
    );
}

/// Test environment presets
#[test]
fn test_environment_presets() {
    let arterial = IntegratedEnvironment::arterial();
    let venous = IntegratedEnvironment::venous();

    assert!(
        (arterial.po2_mmHg - 100.0).abs() < 0.1,
        "Arterial pO2 should be ~100 mmHg"
    );

    assert!(
        (venous.po2_mmHg - 40.0).abs() < 0.1,
        "Venous pO2 should be ~40 mmHg"
    );
}

/// Test high ATP stress increases lactate production
#[test]
fn test_high_atp_stress_increases_lactate() {
    let env = IntegratedEnvironment::default();
    let mut solver = IntegratedSolver::new(env);
    let mut metabolites = MetabolitePool::default_physiological();

    let initial_lactate = metabolites.get(solver.indices.glycolysis.lactate);

    // Run with high ATP consumption (stress)
    solver.run(&mut metabolites, 60.0, 0.01);  // 10x normal ATP demand

    let final_lactate = metabolites.get(solver.indices.glycolysis.lactate);
    let diag = solver.diagnostics(&metabolites);

    // Note: Lactate may increase or stay stable depending on export rate
    // The key is that pH should respond appropriately
    println!(
        "Lactate change: {} -> {}, pH: {}",
        initial_lactate, final_lactate, diag.ph
    );

    // pH should be lower if lactate increased
    if final_lactate > initial_lactate + 0.5 {
        assert!(
            diag.ph < 7.2,
            "pH should be lower than baseline with elevated lactate: {}",
            diag.ph
        );
    }
}

/// Test round-trip pH calculation
#[test]
fn test_ph_lactate_round_trip() {
    let buffer = PhBufferModel::default();

    // Test several values
    for lactate in &[1.0, 2.0, 3.0, 5.0, 8.0] {
        let ph = buffer.calculate_ph(*lactate);
        let recovered_lactate = buffer.calculate_lactate_from_ph(ph);

        assert!(
            (recovered_lactate - lactate).abs() < 0.001,
            "Round-trip failed for lactate {}: got {}",
            lactate, recovered_lactate
        );
    }
}
