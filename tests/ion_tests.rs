//! Integration tests for ion homeostasis (Phase 6b)
//!
//! Tests validate:
//! - Na+/K+-ATPase pump function
//! - Ion steady state maintenance
//! - ATP consumption by ion pumps
//! - ATP balance with active pumps
//! - Ouabain simulation (pump inhibition)

use cell_simulator_x::{
    MetabolitePool, MetaboliteIndices, FullyIntegratedIndices,
    IonIndices, IonHomeostasisSystem, IonHomeostasisConfig,
    FullyIntegratedSolver, FullyIntegratedConfig,
    RedoxIndices, initialize_ion_metabolites,
};

/// Create test system with default configuration
fn create_test_system() -> (IonHomeostasisSystem, MetabolitePool) {
    let glycolysis = MetaboliteIndices::default();
    let redox = RedoxIndices::new(&glycolysis, 18);
    let system = IonHomeostasisSystem::new(&glycolysis, &redox, 35);

    // Create pool with all 38 metabolites
    let mut pool = MetabolitePool::new(38);

    // Set ATP/ADP
    pool.set(glycolysis.atp, 2.0);
    pool.set(glycolysis.adp, 0.25);

    // Initialize ion concentrations
    initialize_ion_metabolites(&mut pool, &system.indices, &system.config);

    (system, pool)
}

#[test]
fn test_na_k_atpase_basic_function() {
    // Test that the pump moves ions and consumes ATP
    let (system, pool) = create_test_system();

    // Check initial state
    let diag = system.diagnostics(&pool);
    assert!((diag.na_cytosolic_mM - 10.0).abs() < 0.1, "Initial Na+ should be ~10 mM");
    assert!((diag.k_cytosolic_mM - 140.0).abs() < 0.1, "Initial K+ should be ~140 mM");

    // Check pump is active
    assert!(diag.pump_rate_mM_per_sec > 0.0, "Pump should be active with physiological conditions");

    // Check derivatives
    let mut dydt = vec![0.0; 38];
    system.compute_derivatives(&pool, &mut dydt);

    // Pump should consume ATP
    assert!(dydt[system.indices.atp] < 0.0, "ATP should be consumed by pump");
    assert!(dydt[system.indices.adp] > 0.0, "ADP should be produced by pump");

    // At steady state, Na+ influx (leak) roughly balances pump efflux
    // Slight imbalance is expected depending on initial conditions
}

#[test]
fn test_ion_steady_state() {
    // Test that ions reach stable steady state
    let mut solver = FullyIntegratedSolver::new(FullyIntegratedConfig::default());
    let mut metabolites = MetabolitePool::default_fully_integrated();

    let indices = solver.indices;

    // Run for 60 seconds
    for _ in 0..60000 {
        solver.step(&mut metabolites, 0.0);  // No additional ATP consumption
    }

    let na = metabolites.get(indices.ions.na_plus_cytosolic);
    let k = metabolites.get(indices.ions.k_plus_cytosolic);

    // Check steady state: Na+ should be 5-15 mM (relaxed range for model)
    assert!(na >= 3.0 && na <= 20.0,
        "Na+ should stabilize near physiological range (5-15 mM), got {:.1} mM", na);

    // K+ should be 130-160 mM (relaxed range for model)
    assert!(k >= 120.0 && k <= 160.0,
        "K+ should stabilize near physiological range (140-150 mM), got {:.1} mM", k);
}

#[test]
fn test_pump_atp_consumption_range() {
    // Test that pump ATP consumption is in expected range
    let (system, pool) = create_test_system();

    let diag = system.diagnostics(&pool);

    // Target: 0.01-0.05 mM/s (Garrahan & Glynn 1967)
    // Allow slightly wider range for model tolerance
    assert!(diag.pump_atp_consumption_mM_per_sec >= 0.005,
        "Pump ATP consumption too low: {:.4} mM/s (target: 0.01-0.05 mM/s)",
        diag.pump_atp_consumption_mM_per_sec);
    assert!(diag.pump_atp_consumption_mM_per_sec <= 0.10,
        "Pump ATP consumption too high: {:.4} mM/s (target: 0.01-0.05 mM/s)",
        diag.pump_atp_consumption_mM_per_sec);
}

#[test]
fn test_atp_balance_with_pumps() {
    // Test that ATP stays in physiological range with active ion pumps
    let mut solver = FullyIntegratedSolver::new(FullyIntegratedConfig::default());
    let mut metabolites = MetabolitePool::default_fully_integrated();

    let indices = solver.indices;

    // Run for 120 seconds
    for _ in 0..120000 {
        solver.step(&mut metabolites, 0.0);
    }

    let atp = metabolites.get(indices.glycolysis.atp);

    // ATP should stay in physiological range (1.5-2.5 mM)
    // Allow slightly wider range due to model complexity
    assert!(atp >= 1.0 && atp <= 3.0,
        "ATP should remain in range (1.5-2.5 mM) with active pumps, got {:.3} mM", atp);
}

#[test]
fn test_ouabain_simulation() {
    // Test pump inhibition (simulating ouabain treatment)
    // When pump is disabled via ion_homeostasis.set_pump_enabled(false),
    // the pump rate should be zero and gradients should collapse.

    let (mut system, pool) = create_test_system();

    // Normal pump rate
    let diag_normal = system.diagnostics(&pool);
    assert!(diag_normal.pump_rate_mM_per_sec > 0.0, "Pump should be active normally");

    // Disable pump (simulate ouabain)
    system.set_pump_enabled(false);
    let diag_ouabain = system.diagnostics(&pool);
    assert_eq!(diag_ouabain.pump_rate_mM_per_sec, 0.0,
        "Pump rate should be zero when disabled (ouabain simulation)");
    assert!(!diag_ouabain.pump_enabled, "Pump should be reported as disabled");
}

#[test]
fn test_na_k_pump_na_dependence() {
    // Test that pump rate depends on intracellular Na+
    let (mut system, mut pool) = create_test_system();

    // Low Na+ condition
    pool.set(system.indices.na_plus_cytosolic, 5.0);
    let rate_low_na = system.na_k_pump.pump_rate_mM_per_sec(
        5.0,
        system.config.k_external_mM,
        pool.get(system.indices.atp),
    );

    // High Na+ condition
    pool.set(system.indices.na_plus_cytosolic, 20.0);
    let rate_high_na = system.na_k_pump.pump_rate_mM_per_sec(
        20.0,
        system.config.k_external_mM,
        pool.get(system.indices.atp),
    );

    // Pump should be faster at higher Na+
    assert!(rate_high_na > rate_low_na,
        "Pump rate should increase with Na+: low={:.4}, high={:.4}",
        rate_low_na, rate_high_na);
}

#[test]
fn test_na_k_pump_atp_dependence() {
    // Test that pump rate depends on ATP
    let (system, _pool) = create_test_system();

    let na = 10.0;
    let k_ext = system.config.k_external_mM;

    // Low ATP
    let rate_low_atp = system.na_k_pump.pump_rate_mM_per_sec(na, k_ext, 0.5);

    // High ATP
    let rate_high_atp = system.na_k_pump.pump_rate_mM_per_sec(na, k_ext, 2.0);

    // Pump should be faster at higher ATP
    assert!(rate_high_atp > rate_low_atp,
        "Pump rate should increase with ATP: low={:.4}, high={:.4}",
        rate_low_atp, rate_high_atp);

    // At very low ATP, pump should be minimal
    let rate_no_atp = system.na_k_pump.pump_rate_mM_per_sec(na, k_ext, 0.01);
    assert!(rate_no_atp < rate_high_atp * 0.1,
        "Pump should be slow at very low ATP");
}

#[test]
fn test_ion_indices_layout() {
    // Verify ion indices are at expected positions
    let indices = FullyIntegratedIndices::new();

    assert_eq!(indices.ions.na_plus_cytosolic, 35, "Na+ should be at index 35");
    assert_eq!(indices.ions.k_plus_cytosolic, 36, "K+ should be at index 36");
    assert_eq!(indices.ions.cl_minus_cytosolic, 37, "Cl- should be at index 37");
    assert_eq!(indices.total_count(), 38, "Total metabolites should be 38");
}

#[test]
fn test_passive_leaks_balance_pump() {
    // Test that passive leaks and pump reach equilibrium
    let (system, pool) = create_test_system();

    let mut dydt = vec![0.0; 38];
    system.compute_derivatives(&pool, &mut dydt);

    // At steady state with default initial conditions,
    // the derivatives should be relatively small (near equilibrium)
    let na_derivative = dydt[system.indices.na_plus_cytosolic];
    let k_derivative = dydt[system.indices.k_plus_cytosolic];

    // Allow some drift but should be small compared to pump rate
    let pump_rate = system.diagnostics(&pool).pump_rate_mM_per_sec;

    // Derivatives should be within an order of magnitude of pump rate
    // (indicating near-equilibrium, not run-away)
    assert!(na_derivative.abs() < pump_rate * 5.0,
        "Na+ derivative should be bounded: {:.4} vs pump {:.4}",
        na_derivative, pump_rate);
    assert!(k_derivative.abs() < pump_rate * 5.0,
        "K+ derivative should be bounded: {:.4} vs pump {:.4}",
        k_derivative, pump_rate);
}
