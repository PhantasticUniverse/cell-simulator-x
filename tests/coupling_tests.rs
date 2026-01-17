//! Integration tests for mechano-metabolic coupling (Phase 8).
//!
//! Tests verify:
//! - Forward coupling: Tension → Piezo1 → Ca²⁺ → ATP release
//! - Reverse coupling: ATP → spectrin stiffness
//! - Timestep synchronization between physics and biochemistry

use cell_simulator_x::{
    CoupledSolver, CoupledConfig, SpectrinModulator, TensionComputer,
    MetabolitePool, FullyIntegratedConfig, FullyIntegratedSolver,
};
use cell_simulator_x::config::GeometryParameters;
use cell_simulator_x::geometry::{Mesh, SpectrinNetwork};
use cell_simulator_x::physics::{SkalakSolver, SkalakMaterial};
use cell_simulator_x::state::PhysicsState;

fn create_test_geometry() -> (Mesh, SpectrinNetwork, GeometryParameters) {
    let params = GeometryParameters {
        cell_radius_um: 3.91,
        fung_tong_c0_um: 0.81,
        fung_tong_c2_um: 7.83,
        fung_tong_c4_um: -4.39,
        mesh_resolution: 8,
        spectrin_target_count: 50,
    };
    let mesh = Mesh::generate_rbc(&params);
    let spectrin = SpectrinNetwork::generate(&mesh, &params);
    (mesh, spectrin, params)
}

// ============================================================================
// TensionComputer Tests
// ============================================================================

#[test]
fn test_tension_at_rest() {
    let (mesh, _, _) = create_test_geometry();
    let skalak = SkalakSolver::new(&mesh, SkalakMaterial::default());
    let mut computer = TensionComputer::default();

    // At rest, strain invariants should be ~0, so tension ~0
    let tension = computer.compute_global_tension_pN_per_nm(&skalak, &mesh);

    // Tension should be very small at equilibrium
    assert!(
        tension < 0.5,
        "Tension at rest should be near zero, got: {:.4} pN/nm",
        tension
    );
}

#[test]
fn test_tension_computer_reset() {
    let (mesh, _, _) = create_test_geometry();
    let skalak = SkalakSolver::new(&mesh, SkalakMaterial::default());
    let mut computer = TensionComputer::default();

    // Compute tension to populate history
    computer.compute_global_tension_pN_per_nm(&skalak, &mesh);
    computer.compute_global_tension_pN_per_nm(&skalak, &mesh);

    // Instantaneous tension should be available
    let tension = computer.instantaneous_tension();
    assert!(tension >= 0.0);

    // Reset should clear history
    computer.reset();
    assert!(computer.averaged_tension().abs() < 0.001);
}

// ============================================================================
// SpectrinModulator Tests
// ============================================================================

#[test]
fn test_stiffness_modifier_normal_atp() {
    let modulator = SpectrinModulator::default();

    // At reference ATP (2.0 mM), modifier should be 1.0
    let modifier = modulator.stiffness_modifier(2.0);
    assert!(
        (modifier - 1.0).abs() < 0.01,
        "Modifier at normal ATP should be 1.0, got {}",
        modifier
    );
}

#[test]
fn test_stiffness_modifier_low_atp() {
    let modulator = SpectrinModulator::default();

    // At 0.5 mM ATP (75% depleted), modifier should be ~1.375
    // Formula: 1.0 + 0.5 * (1.0 - 0.5/2.0) = 1.0 + 0.5 * 0.75 = 1.375
    let modifier = modulator.stiffness_modifier(0.5);
    assert!(
        modifier > 1.3,
        "Low ATP should increase stiffness modifier, got {}",
        modifier
    );
    assert!(
        modifier < 1.5,
        "Modifier should not exceed max, got {}",
        modifier
    );
}

#[test]
fn test_stiffness_modifier_zero_atp() {
    let modulator = SpectrinModulator::default();

    // At zero ATP, modifier should be 1.5 (maximum stiffening)
    let modifier = modulator.stiffness_modifier(0.0);
    assert!(
        (modifier - 1.5).abs() < 0.01,
        "Modifier at zero ATP should be 1.5, got {}",
        modifier
    );
}

#[test]
fn test_modified_shear_modulus() {
    let modulator = SpectrinModulator::default();
    let base_gs = 5.5; // μN/m

    // At normal ATP, shear modulus unchanged
    let gs_normal = modulator.modified_shear_modulus(base_gs, 2.0);
    assert!((gs_normal - base_gs).abs() < 0.01);

    // At zero ATP, shear modulus increased by 1.5x
    let gs_depleted = modulator.modified_shear_modulus(base_gs, 0.0);
    let expected = base_gs * 1.5;
    assert!(
        (gs_depleted - expected).abs() < 0.1,
        "Expected {}, got {}",
        expected,
        gs_depleted
    );
}

// ============================================================================
// CoupledSolver Tests
// ============================================================================

#[test]
fn test_coupled_solver_creation() {
    let (mesh, _, _) = create_test_geometry();
    let metabolites = MetabolitePool::default_fully_integrated();
    let config = CoupledConfig::default();
    let solver = CoupledSolver::new(&mesh, config);

    assert!(solver.time_sec == 0.0);
    // Check stiffness modifier via diagnostics
    let diag = solver.diagnostics(&metabolites);
    assert!((diag.stiffness_modifier - 1.0).abs() < 0.01);
}

#[test]
fn test_coupled_solver_no_crash() {
    let (mesh, spectrin, _) = create_test_geometry();
    let mut mesh = mesh;
    let mut physics_state = PhysicsState::new(mesh.vertices.len());
    let mut metabolites = MetabolitePool::default_fully_integrated();

    let mut config = CoupledConfig::default();
    config.physics_substeps = 10; // Reduce for faster test
    let mut solver = CoupledSolver::new(&mesh, config);

    // Run for 10 steps - should not crash
    for _ in 0..10 {
        solver.step(&mut mesh, &spectrin, &mut physics_state, &mut metabolites);
    }

    assert!(solver.time_sec > 0.0);
}

#[test]
fn test_piezo1_responds_to_tension_override() {
    // Test using the biochemistry solver directly for cleaner results
    let mut config = FullyIntegratedConfig::default();
    config.membrane_tension_pN_per_nm = 5.0; // Strong tension to activate Piezo1
    config.enable_piezo1 = true;

    let mut solver = FullyIntegratedSolver::new(config);
    let mut metabolites = MetabolitePool::default_fully_integrated();

    // Get initial Ca²⁺
    let initial_ca = metabolites.get(solver.indices.redox.ca2_plus_cytosolic);

    // Run for 500 ms with tension
    solver.run(&mut metabolites, 0.5, 0.0);

    let final_ca = metabolites.get(solver.indices.redox.ca2_plus_cytosolic);

    // Ca²⁺ should have increased due to Piezo1 activation
    // Even a small increase indicates the pathway is working
    assert!(
        final_ca > initial_ca,
        "Expected Ca²⁺ increase from tension. Initial: {:.4}, Final: {:.4}",
        initial_ca,
        final_ca
    );
}

#[test]
fn test_atp_depletion_increases_stiffness() {
    let (mesh, _, _) = create_test_geometry();
    let mut metabolites = MetabolitePool::default_fully_integrated();
    let config = CoupledConfig::default();
    let solver = CoupledSolver::new(&mesh, config);

    // At normal ATP, modifier should be 1.0
    let diag_normal = solver.diagnostics(&metabolites);
    assert!(
        (diag_normal.stiffness_modifier - 1.0).abs() < 0.1,
        "Normal ATP should have modifier ~1.0"
    );

    // Deplete ATP
    metabolites.set(solver.biochemistry.indices.glycolysis.atp, 0.5);

    // Manually check modifier from modulator
    let atp_mM = metabolites.get(solver.biochemistry.indices.glycolysis.atp);
    let modifier = solver.spectrin_modulator.stiffness_modifier(atp_mM);

    assert!(
        modifier > 1.3,
        "Low ATP ({:.1} mM) should increase stiffness modifier to >1.3, got {}",
        atp_mM,
        modifier
    );
}

// ============================================================================
// Bidirectional Coupling Tests
// ============================================================================

#[test]
fn test_forward_coupling_tension_to_biochemistry() {
    // Test that membrane tension affects biochemistry through Piezo1
    let mut solver = FullyIntegratedSolver::new(FullyIntegratedConfig::default());
    let mut metabolites = MetabolitePool::default_fully_integrated();

    // Run without tension
    solver.config.membrane_tension_pN_per_nm = 0.0;
    solver.run(&mut metabolites, 0.1, 0.0);
    let ca_no_tension = metabolites.get(solver.indices.redox.ca2_plus_cytosolic);

    // Reset and run with tension
    let mut metabolites2 = MetabolitePool::default_fully_integrated();
    let mut solver2 = FullyIntegratedSolver::new(FullyIntegratedConfig::default());
    solver2.config.membrane_tension_pN_per_nm = 3.0; // Significant tension

    solver2.run(&mut metabolites2, 0.1, 0.0);
    let ca_with_tension = metabolites2.get(solver2.indices.redox.ca2_plus_cytosolic);

    // Ca²⁺ should be higher with tension
    assert!(
        ca_with_tension > ca_no_tension,
        "Tension should increase Ca²⁺ via Piezo1. Without: {:.3} μM, With: {:.3} μM",
        ca_no_tension,
        ca_with_tension
    );
}

#[test]
fn test_reverse_coupling_atp_to_stiffness() {
    // Test that ATP level affects mechanical stiffness
    let modulator = SpectrinModulator::default();

    // Create two conditions: normal ATP and depleted ATP
    let normal_atp = 2.0;  // mM
    let low_atp = 0.5;     // mM (storage lesion level)

    let stiffness_normal = modulator.stiffness_modifier(normal_atp);
    let stiffness_depleted = modulator.stiffness_modifier(low_atp);

    // Stiffness should increase as ATP depletes
    assert!(
        stiffness_depleted > stiffness_normal,
        "Low ATP should increase stiffness. Normal: {}, Low: {}",
        stiffness_normal,
        stiffness_depleted
    );

    // Verify magnitudes
    assert!((stiffness_normal - 1.0).abs() < 0.01, "Normal should be 1.0");
    assert!(stiffness_depleted > 1.3, "Low ATP should be >1.3");
}

#[test]
fn test_monotonic_stiffness_with_atp() {
    let modulator = SpectrinModulator::default();

    // As ATP decreases from 2.0 to 0.0, stiffness should monotonically increase
    let atp_levels = [2.0, 1.5, 1.0, 0.5, 0.0];
    let mut prev_modifier = 0.0;

    for atp in atp_levels {
        let modifier = modulator.stiffness_modifier(atp);
        assert!(
            modifier >= prev_modifier,
            "Stiffness should increase monotonically as ATP decreases"
        );
        prev_modifier = modifier;
    }
}

// ============================================================================
// Timestep Synchronization Tests
// ============================================================================

#[test]
fn test_timestep_ratio() {
    let config = CoupledConfig::default();

    // Default should be 1000 physics steps per biochem step
    // (1 μs physics / 1 ms biochem)
    let expected_ratio = (config.biochem_dt_sec / config.physics_dt_sec) as usize;
    assert_eq!(
        config.physics_substeps, expected_ratio,
        "Substeps should match timestep ratio"
    );
}
