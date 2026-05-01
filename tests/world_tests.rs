//! Phase 10.5 integration tests for the multi-cell `World` API.
//!
//! Exit criteria from the Phase 10.5 plan:
//! - `World::step` runs N=1, N=10, N=100 cells on CPU.
//! - All Phase 10 validation tests still pass at N=1.
//! - Per-cell results stay consistent across N (no aliasing between cells).

use cell_simulator_x::config::{GeometryParameters, Parameters};
use cell_simulator_x::coupling::CoupledConfig;
use cell_simulator_x::world::World;

fn test_params() -> Parameters {
    Parameters {
        geometry: GeometryParameters {
            cell_radius_um: 3.91,
            fung_tong_c0_um: 0.81,
            fung_tong_c2_um: 7.83,
            fung_tong_c4_um: -4.39,
            // Small mesh keeps per-step cost low so N=100 finishes
            // quickly in CI. 8 → ~250 vertices.
            mesh_resolution: 8,
            spectrin_target_count: 50,
        },
        ..Parameters::default()
    }
}

fn test_config() -> CoupledConfig {
    CoupledConfig {
        // Reduced substeps so N=100 finishes in seconds, not minutes.
        physics_substeps: 5,
        ..CoupledConfig::default()
    }
}

#[test]
fn world_n_equals_one_runs_one_step() {
    let mut world = World::new(test_config());
    let h = world.add_cell(&test_params());
    world.step();
    let cell = world.cell(h);
    assert!(cell.time_sec > 0.0, "cell time should advance after step");
    assert!(world.time_sec > 0.0, "world time should advance with cells");
}

#[test]
fn world_n_equals_ten_steps_in_lockstep() {
    let mut world = World::new(test_config());
    let mut handles = Vec::with_capacity(10);
    for _ in 0..10 {
        handles.push(world.add_cell(&test_params()));
    }
    world.step();
    let t0 = world.cell(handles[0]).time_sec;
    for h in &handles {
        let t = world.cell(*h).time_sec;
        assert!((t - t0).abs() < 1e-12, "all cells must step in lockstep");
    }
}

#[test]
fn world_n_equals_one_hundred_runs_without_panic() {
    // Demonstrates the multi-cell API at the upper end of the Phase 10.5
    // exit criterion. With reduced physics_substeps and small mesh
    // resolution, this completes well under the test timeout.
    let mut world = World::new(test_config());
    for _ in 0..100 {
        world.add_cell(&test_params());
    }
    assert_eq!(world.len(), 100);

    // One step is enough — the goal here is "expressible at scale", not
    // "long sim". Real performance characterization is in
    // `--diagnose-multi-cell`.
    world.step();
    assert!(world.time_sec > 0.0);
}

#[test]
fn world_atp_matches_across_independent_cells() {
    // After running for a while, every cell's ATP should be in the same
    // physiological band — they are independent runs of the same model
    // with the same starting conditions, only diverging through DPD
    // random forces.
    let mut world = World::new(test_config());
    for _ in 0..5 {
        world.add_cell(&test_params());
    }

    // Run 0.1 s of simulated time (100 biochem steps).
    world.run(0.1);

    let atp_idx = world.cells()[0].biochemistry.indices.glycolysis.atp;
    let atps: Vec<f64> = world
        .cells()
        .iter()
        .map(|c| c.metabolites.get(atp_idx))
        .collect();
    let mean = atps.iter().sum::<f64>() / atps.len() as f64;
    let max_dev = atps.iter().map(|v| (v - mean).abs()).fold(0.0f64, f64::max);

    assert!(
        max_dev < 0.05,
        "ATP across cells should agree within 0.05 mM, got max deviation {} mM (atps: {:?})",
        max_dev, atps,
    );
    assert!(
        mean >= 1.0 && mean <= 2.5,
        "Mean ATP across cells should be physiological (1.0–2.5 mM), got {}",
        mean
    );
}

#[test]
fn world_handles_remain_stable() {
    let mut world = World::new(test_config());
    let h0 = world.add_cell(&test_params());
    let h1 = world.add_cell(&test_params());
    let h2 = world.add_cell(&test_params());

    assert_ne!(h0, h1);
    assert_ne!(h1, h2);
    assert_ne!(h0, h2);
    assert_eq!(h0.index(), 0);
    assert_eq!(h1.index(), 1);
    assert_eq!(h2.index(), 2);
}

#[test]
fn world_reset_clears_time() {
    let mut world = World::new(test_config());
    world.add_cell(&test_params());
    world.step();
    assert!(world.time_sec > 0.0);

    world.reset();
    assert_eq!(world.time_sec, 0.0);
    for cell in world.cells() {
        assert_eq!(cell.time_sec, 0.0);
    }
}

#[test]
fn world_tension_override_applies_per_cell() {
    let mut world = World::new(test_config());
    world.add_cell(&test_params());
    world.add_cell(&test_params());

    world.set_tension_override(2.5);
    for cell in world.cells() {
        assert!((cell.current_tension_pN_per_nm - 2.5).abs() < 1e-9);
    }
}
