//! Phase 12.B.1 parity test: PhysicsBackend external_forces buffer.
//!
//! With external_forces set to a uniform constant force per vertex, GPU
//! `step()` must produce vertex displacement that matches a CPU
//! reference loop adding the same constant force inside the integrator.
//!
//! Tolerance target (per Phase 12.B plan): <1e-6 μm after 100 substeps.

use cell_simulator_x::compute::{ComputeContext, PhysicsBackend, PhysicsBackendConfig};
use cell_simulator_x::config::GeometryParameters;
use cell_simulator_x::geometry::{Mesh, SpectrinNetwork};
use cell_simulator_x::physics::{
    DPDParameters, DPDSolver, SkalakMaterial, SkalakSolver, VelocityVerlet, WLCParameters,
    WLCSolver,
};
use glam::Vec3;

const N_SUBSTEPS: u32 = 100;
const DT_SEC: f32 = 1e-6;
const TEMPERATURE_K: f32 = 310.0;
const SEED: u32 = 0xCAFEBABE;
const MEMBRANE_DAMPING: f32 = 0.1;

/// Constant external force per vertex (in μN). Magnitude chosen so the
/// resulting displacement is observable but well within max_displacement.
fn external_force_per_vertex() -> Vec3 {
    Vec3::new(0.0, 0.0, 1.0e-3) // +z drag, 1 nN per vertex
}

fn default_params() -> GeometryParameters {
    GeometryParameters {
        cell_radius_um: 3.91,
        fung_tong_c0_um: 0.81,
        fung_tong_c2_um: 7.83,
        fung_tong_c4_um: -4.39,
        mesh_resolution: 10,
        spectrin_target_count: 100,
    }
}

fn run_cpu_reference_with_external(
    mesh: &mut Mesh,
    spectrin: &SpectrinNetwork,
    external: Vec3,
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

        // Reset forces, then add external + WLC + Skalak + DPD.
        for f in forces.iter_mut() {
            *f = external;
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
fn external_forces_uniform_constant_parity() {
    let ctx = match ComputeContext::new_headless_blocking() {
        Ok(c) => c,
        Err(e) => {
            eprintln!("skipping GPU parity test: no adapter: {e}");
            return;
        }
    };

    let mut cpu_mesh = Mesh::generate_rbc(&default_params());
    let mut gpu_mesh = Mesh::generate_rbc(&default_params());
    let spectrin = SpectrinNetwork::generate(&cpu_mesh, &default_params());

    let external = external_force_per_vertex();
    let n = cpu_mesh.vertices.len();

    let _cpu_velocities = run_cpu_reference_with_external(&mut cpu_mesh, &spectrin, external, N_SUBSTEPS);

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

    // Upload constant external force, then run substeps. Set once before
    // the loop — buffer persists across substeps.
    let external_per_vertex = vec![external; n];
    backend.set_external_forces(&external_per_vertex).expect("set_external_forces");
    for _ in 0..N_SUBSTEPS {
        backend.step().expect("backend step");
    }
    let (gpu_positions, _gpu_velocities) = backend.read_state().expect("read_state");

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
    println!(
        "external_forces uniform-constant parity: n_vertices = {}, n_substeps = {}",
        n, N_SUBSTEPS
    );
    println!(
        "  worst |Δpos| = {:.6e} μm at vertex {} (cpu = {:?} gpu = {:?})",
        max_pos_err, max_idx,
        cpu_mesh.vertices[max_idx].position_vec3(),
        gpu_positions[max_idx]
    );

    assert!(
        max_pos_err < 1e-5,
        "external_forces parity: position drift {} μm at vertex {} > 1e-5 μm",
        max_pos_err, max_idx
    );
}

#[test]
fn external_forces_zero_matches_baseline_parity() {
    // With external_forces = 0, the backend must produce IDENTICAL
    // results to the existing physics_backend_parity test (since the
    // new buffer just adds 0). This is a regression-safety check.
    let ctx = match ComputeContext::new_headless_blocking() {
        Ok(c) => c,
        Err(e) => {
            eprintln!("skipping GPU parity test: no adapter: {e}");
            return;
        }
    };

    let mut cpu_mesh = Mesh::generate_rbc(&default_params());
    let mut gpu_mesh = Mesh::generate_rbc(&default_params());
    let spectrin = SpectrinNetwork::generate(&cpu_mesh, &default_params());

    let _cpu_velocities = run_cpu_reference_with_external(&mut cpu_mesh, &spectrin, Vec3::ZERO, N_SUBSTEPS);

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
    // Don't set external_forces — buffer is zero-initialized.
    for _ in 0..N_SUBSTEPS {
        backend.step().expect("backend step");
    }
    let (gpu_positions, _) = backend.read_state().expect("read_state");

    let n = cpu_mesh.vertices.len();
    let mut max_err = 0.0f32;
    for i in 0..n {
        let cp = cpu_mesh.vertices[i].position_vec3();
        let gp = gpu_positions[i];
        let err = (cp - gp).length();
        if err > max_err {
            max_err = err;
        }
    }
    println!("zero-external parity worst |Δpos| = {:.6e} μm", max_err);
    assert!(max_err < 1e-5, "zero-external parity drift {} μm > 1e-5", max_err);
}
