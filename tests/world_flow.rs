//! Phase 12.B.2: `World::apply_poiseuille_drag` integration test.
//!
//! Verifies that calling `apply_poiseuille_drag` followed by `step()`
//! produces a directional cell drift along the flow axis, matching the
//! single-cell behavior validated in `flow_deformation.rs`.

#![allow(deprecated)]

use cell_simulator_x::config::{GeometryParameters, Parameters};
use cell_simulator_x::coupling::CoupledConfig;
use cell_simulator_x::flow::{CylindricalChannel, Poiseuille};
use cell_simulator_x::world::World;
use glam::Vec3;

fn small_params() -> Parameters {
    Parameters {
        geometry: GeometryParameters {
            cell_radius_um: 3.91,
            fung_tong_c0_um: 0.81,
            fung_tong_c2_um: 7.83,
            fung_tong_c4_um: -4.39,
            mesh_resolution: 6,
            spectrin_target_count: 30,
        },
        ..Parameters::default()
    }
}

fn small_config() -> CoupledConfig {
    CoupledConfig {
        physics_substeps: 5,
        ..CoupledConfig::default()
    }
}

fn cell_centroid(cell: &cell_simulator_x::world::Cell) -> Vec3 {
    let mut sum = Vec3::ZERO;
    for v in &cell.mesh.vertices {
        sum += v.position_vec3();
    }
    sum / cell.mesh.vertices.len() as f32
}

#[test]
fn apply_poiseuille_drag_writes_to_all_cells() {
    let mut world = World::new(small_config());
    let params = small_params();
    for _ in 0..3 {
        world.add_cell(&params);
    }
    let flow = Poiseuille::new(CylindricalChannel::default(), 100.0);
    world.apply_poiseuille_drag(&flow, 0.5);
    for cell in world.cells() {
        // Each cell should have non-empty external-forces buffer.
        assert_eq!(
            cell.physics_state.external_forces_uN.len(),
            cell.mesh.vertices.len()
        );
        // At least some forces should be non-zero (cells are inside
        // the 5 µm channel by default).
        let any_nonzero = cell.physics_state.external_forces_uN
            .iter()
            .any(|f| f.length() > 1e-6);
        assert!(any_nonzero, "expected non-zero drag forces on each cell");
    }
}

#[test]
fn cells_drift_along_flow_axis_under_drag() {
    let mut world = World::new(small_config());
    let params = small_params();
    for _ in 0..3 {
        world.add_cell(&params);
    }
    let flow = Poiseuille::new(CylindricalChannel::default(), 100.0);

    // Record initial centroids.
    let initial: Vec<Vec3> = world.cells().iter().map(cell_centroid).collect();

    // Apply drag, step a few times.
    for _ in 0..30 {
        world.apply_poiseuille_drag(&flow, 0.5);
        world.step();
    }

    // Each cell's centroid must have shifted in +z (the flow axis).
    for (i, cell) in world.cells().iter().enumerate() {
        let final_pos = cell_centroid(cell);
        let dz = final_pos.z - initial[i].z;
        println!("cell {}: dz = {:.6} µm", i, dz);
        assert!(dz > 0.0, "cell {} did not drift along flow axis (dz = {})", i, dz);
    }
}

#[test]
fn no_drift_with_zero_velocity_flow() {
    let mut world = World::new(small_config());
    let params = small_params();
    world.add_cell(&params);
    // Zero-velocity flow → no drag.
    let flow = Poiseuille::new(CylindricalChannel::default(), 0.0);

    let initial = cell_centroid(&world.cells()[0]);
    for _ in 0..30 {
        world.apply_poiseuille_drag(&flow, 0.5);
        world.step();
    }
    let final_pos = cell_centroid(&world.cells()[0]);
    let dz = final_pos.z - initial.z;
    println!("zero-flow dz = {:.6} µm", dz);
    // Allow microns of thermal noise drift.
    assert!(dz.abs() < 1.0, "unexpected drift {} under zero flow", dz);
}
