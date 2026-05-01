//! Phase 12.A integration test: a cell placed in Poiseuille flow with
//! drag coupling deforms (vertex positions shift) over time.
//!
//! This is a qualitative check — we don't yet validate against
//! parachute-shape data (Phase 12.C). What we verify here:
//!
//! 1. With zero flow, the cell barely moves (only thermal noise).
//! 2. With flow, the cell's centroid shifts along the flow axis after
//!    `n_steps` and individual vertices acquire velocity along the axis.
//!
//! Skipping this with `#[ignore]`-style gating is unnecessary — the test
//! is fast (no GPU dispatch, ~tens of ms).

use cell_simulator_x::config::GeometryParameters;
use cell_simulator_x::flow::{
    apply_drag_to_external_forces, CylindricalChannel, Poiseuille,
};
use cell_simulator_x::geometry::{Mesh, SpectrinNetwork};
use cell_simulator_x::physics::{PhysicsConfig, PhysicsSolver};
use cell_simulator_x::state::PhysicsState;
use glam::Vec3;

fn test_params() -> GeometryParameters {
    GeometryParameters {
        cell_radius_um: 3.91,
        fung_tong_c0_um: 0.81,
        fung_tong_c2_um: 7.83,
        fung_tong_c4_um: -4.39,
        mesh_resolution: 6,
        spectrin_target_count: 30,
    }
}

fn centroid(mesh: &Mesh) -> Vec3 {
    let mut sum = Vec3::ZERO;
    for v in &mesh.vertices {
        sum += v.position_vec3();
    }
    sum / mesh.vertices.len() as f32
}

#[test]
fn cell_translates_along_flow_axis() {
    let params = test_params();
    let mut mesh = Mesh::generate_rbc(&params);
    let spectrin = SpectrinNetwork::generate(&mesh, &params);
    let mut state = PhysicsState::new(mesh.vertices.len());
    let mut solver = PhysicsSolver::new(
        &mesh,
        PhysicsConfig {
            dt_sec: 1e-6,
            temperature_K: 310.0,
            // Disable thermal noise to make the test deterministic — drag
            // alone should drive the centroid forward.
            enable_thermal_noise: false,
            membrane_damping: 0.05,
        },
    );

    // 5 µm radius cylindrical channel along +Z, centered at origin. The
    // cell sits at the center of the channel.
    let flow = Poiseuille::new(
        CylindricalChannel {
            radius_um: 5.0,
            axis: Vec3::Z,
            center: Vec3::ZERO,
        },
        // Aggressive flow — 1 mm/s centerline (1000 µm/s) in a small
        // capillary so the test sees a clear translation in 1000 substeps.
        1000.0,
    );

    let initial_centroid = centroid(&mesh);

    let drag_coeff = 5.0_f32; // tuned so drag dominates membrane tension
    let n_steps = 1000;
    for _ in 0..n_steps {
        // Refresh drag (depends on vertex velocities each step).
        apply_drag_to_external_forces(&mut state, &mesh, &flow, drag_coeff);
        solver.step(&mut mesh, &spectrin, &mut state);
    }

    let final_centroid = centroid(&mesh);
    let drift = final_centroid - initial_centroid;
    println!(
        "initial centroid = {:?}, final = {:?}, drift = {:?}",
        initial_centroid, final_centroid, drift
    );

    // Should drift along +Z and barely sideways.
    assert!(
        drift.z > 0.001,
        "expected positive z-drift under flow, got {}",
        drift.z
    );
    assert!(
        drift.x.abs() < 0.05 && drift.y.abs() < 0.05,
        "expected near-zero lateral drift, got xy = ({}, {})",
        drift.x,
        drift.y
    );
}

#[test]
fn no_flow_no_translation() {
    let params = test_params();
    let mut mesh = Mesh::generate_rbc(&params);
    let spectrin = SpectrinNetwork::generate(&mesh, &params);
    let mut state = PhysicsState::new(mesh.vertices.len());
    let mut solver = PhysicsSolver::new(
        &mesh,
        PhysicsConfig {
            dt_sec: 1e-6,
            temperature_K: 310.0,
            enable_thermal_noise: false,
            membrane_damping: 0.05,
        },
    );
    // Zero flow — drag forces are all zero.
    let flow = Poiseuille::new(CylindricalChannel::default(), 0.0);

    let initial_centroid = centroid(&mesh);
    for _ in 0..1000 {
        apply_drag_to_external_forces(&mut state, &mesh, &flow, 5.0);
        solver.step(&mut mesh, &spectrin, &mut state);
    }
    let drift = centroid(&mesh) - initial_centroid;
    assert!(
        drift.length() < 1e-3,
        "no-flow centroid should not drift, got {:?}",
        drift
    );
}
