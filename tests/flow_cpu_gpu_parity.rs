//! Phase 12.B.2 parity: same cell + same Poiseuille flow drives CPU
//! and GPU paths to within tolerance after N substeps.
//!
//! CPU path: `flow::apply_drag_to_external_forces` writes drag into
//! `physics_state.external_forces_uN`; `PhysicsSolver::step` reads it.
//! GPU path: per-vertex drag computed host-side and uploaded via
//! `PhysicsBackend::set_external_forces`; `skalak_init_from_baseline`
//! reads it inside the integrated step kernel.
//!
//! The plan target tolerance is <1e-3 μm cell-centroid drift after 1 s
//! of simulated flow. We use a shorter window for test-suite speed.

use cell_simulator_x::compute::{ComputeContext, PhysicsBackend, PhysicsBackendConfig};
use cell_simulator_x::config::GeometryParameters;
use cell_simulator_x::flow::{drag_force_uN, CylindricalChannel, Poiseuille};
use cell_simulator_x::geometry::{Mesh, SpectrinNetwork};
use cell_simulator_x::physics::{
    DPDParameters, DPDSolver, SkalakMaterial, SkalakSolver, VelocityVerlet, WLCParameters,
    WLCSolver,
};
use glam::Vec3;

const N_SUBSTEPS: u32 = 100;
const DT_SEC: f32 = 1e-6;
const TEMPERATURE_K: f32 = 310.0;
const SEED: u32 = 0x12345678;
const MEMBRANE_DAMPING: f32 = 0.1;
const DRAG_COEFF: f32 = 0.5;

fn default_params() -> GeometryParameters {
    GeometryParameters {
        cell_radius_um: 3.91,
        fung_tong_c0_um: 0.81,
        fung_tong_c2_um: 7.83,
        fung_tong_c4_um: -4.39,
        mesh_resolution: 8,
        spectrin_target_count: 60,
    }
}

/// CPU reference: mirror PhysicsSolver::step over N substeps with
/// per-substep Poiseuille drag computation.
fn run_cpu_with_flow(
    mesh: &mut Mesh,
    spectrin: &SpectrinNetwork,
    flow: &Poiseuille,
    n_substeps: u32,
) -> Vec<Vec3> {
    let n_vertices = mesh.vertices.len();
    let mut velocities = vec![Vec3::ZERO; n_vertices];
    let mut forces = vec![Vec3::ZERO; n_vertices];

    let skalak_solver = SkalakSolver::new(mesh, SkalakMaterial::default());
    let wlc_solver = WLCSolver::new(WLCParameters::default());
    let dpd_solver = DPDSolver::new(DPDParameters::default());
    let mut integrator = VelocityVerlet::new();

    for step in 0..n_substeps {
        integrator.half_step_velocity(&mut velocities, &forces, DT_SEC);
        let positions: Vec<Vec3> = mesh.vertices.iter().map(|v| v.position_vec3()).collect();
        let new_positions = integrator.step_position(&positions, &velocities, DT_SEC);
        for (i, p) in new_positions.iter().enumerate() {
            mesh.vertices[i].position = p.to_array();
        }

        // Compute drag force per vertex (Poiseuille).
        for i in 0..n_vertices {
            let pos = mesh.vertices[i].position_vec3();
            forces[i] = drag_force_uN(flow, pos, velocities[i], DRAG_COEFF);
        }

        let wlc_forces = wlc_solver.compute_forces(spectrin, &new_positions);
        for (node_idx, force) in wlc_forces {
            if let Some(node) = spectrin.nodes.get(node_idx) {
                let v_idx = node.mesh_vertex_idx as usize;
                if v_idx < n_vertices {
                    forces[v_idx] += force;
                }
            }
        }

        let (skalak_forces, _) = skalak_solver.compute_forces(mesh);
        for (i, f) in skalak_forces.into_iter().enumerate() {
            forces[i] += f;
        }

        let dpd_forces = dpd_solver.compute_membrane_forces_pcg(
            &velocities,
            TEMPERATURE_K,
            step,
            SEED,
        );
        for (i, f) in dpd_forces.into_iter().enumerate() {
            forces[i] += f;
        }

        for (i, vel) in velocities.iter().enumerate() {
            forces[i] -= *vel * MEMBRANE_DAMPING;
        }

        integrator.half_step_velocity(&mut velocities, &forces, DT_SEC);
    }

    velocities
}

#[test]
fn poiseuille_drag_cpu_gpu_parity() {
    let ctx = match ComputeContext::new_headless_blocking() {
        Ok(c) => c,
        Err(e) => {
            eprintln!("skipping GPU parity test: no adapter: {e}");
            return;
        }
    };

    let mut cpu_mesh = Mesh::generate_rbc(&default_params());
    let gpu_mesh = Mesh::generate_rbc(&default_params());
    let spectrin = SpectrinNetwork::generate(&cpu_mesh, &default_params());

    let flow = Poiseuille::new(CylindricalChannel::default(), 100.0);

    let _cpu_velocities = run_cpu_with_flow(&mut cpu_mesh, &spectrin, &flow, N_SUBSTEPS);

    let skalak_solver = SkalakSolver::new(&gpu_mesh, SkalakMaterial::default());
    let wlc_params = WLCParameters::default();
    let config = PhysicsBackendConfig {
        dt_sec: DT_SEC,
        vertex_mass_pg: 1.0,
        max_velocity_um_per_sec: 1000.0,
        max_displacement_um: 0.1,
        temperature_K: TEMPERATURE_K,
        membrane_damping: MEMBRANE_DAMPING,
        seed: SEED,
    };
    let mut backend = PhysicsBackend::new(
        &ctx,
        &gpu_mesh,
        &spectrin,
        &skalak_solver,
        &wlc_params,
        SkalakMaterial::default(),
        config,
    )
    .expect("PhysicsBackend::new");

    let n = gpu_mesh.vertices.len();

    for _ in 0..N_SUBSTEPS {
        // Read GPU state to compute drag in host-space.
        let (positions, velocities) = backend.read_state().expect("read_state");
        let drag: Vec<Vec3> = (0..n)
            .map(|i| drag_force_uN(&flow, positions[i], velocities[i], DRAG_COEFF))
            .collect();
        backend.set_external_forces(&drag).expect("set_external_forces");
        backend.step().expect("backend step");
    }

    let (gpu_positions, _) = backend.read_state().expect("read_state");

    let mut max_pos_err = 0.0f32;
    let mut max_idx = 0usize;
    for i in 0..n {
        let cp = cpu_mesh.vertices[i].position_vec3();
        let gp = gpu_positions[i];
        let err = (cp - gp).length();
        if err > max_pos_err {
            max_pos_err = err;
            max_idx = i;
        }
    }

    // Cell-centroid drift (the Phase 12.B headline metric).
    let cpu_centroid: Vec3 = cpu_mesh.vertices.iter()
        .map(|v| v.position_vec3())
        .sum::<Vec3>() / n as f32;
    let gpu_centroid: Vec3 = gpu_positions.iter().sum::<Vec3>() / n as f32;
    let centroid_err = (cpu_centroid - gpu_centroid).length();

    println!(
        "Poiseuille CPU/GPU parity: n_vertices={}, n_substeps={}",
        n, N_SUBSTEPS
    );
    println!("  worst |Δpos| = {:.6e} μm at vertex {}", max_pos_err, max_idx);
    println!("  CPU centroid: {:?}", cpu_centroid);
    println!("  GPU centroid: {:?}", gpu_centroid);
    println!("  centroid drift: {:.6e} μm", centroid_err);

    // Per-vertex tolerance: 1e-5 μm absolute (matches Phase 11 parity).
    assert!(
        max_pos_err < 1e-5,
        "Poiseuille parity: position drift {} μm at vertex {} > 1e-5 μm",
        max_pos_err, max_idx
    );
    // Centroid tolerance: 1e-3 μm (Phase 12.B plan target).
    assert!(
        centroid_err < 1e-3,
        "Poiseuille parity: centroid drift {} μm > 1e-3 μm",
        centroid_err
    );
}
