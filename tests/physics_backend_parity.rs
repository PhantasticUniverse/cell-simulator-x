//! Phase 11.3.E parity test: GPU PhysicsBackend vs a CPU equivalent loop
//! that mirrors `PhysicsSolver::step` (with PCG-DPD instead of `StdRng`-DPD).
//!
//! Asserts per-vertex position and velocity match within tolerance after
//! N substeps. Skipped if no GPU adapter is available.

use cell_simulator_x::compute::{ComputeContext, PhysicsBackend, PhysicsBackendConfig};
use cell_simulator_x::config::GeometryParameters;
use cell_simulator_x::geometry::{Mesh, SpectrinNetwork};
use cell_simulator_x::physics::{
    DPDParameters, DPDSolver, SkalakMaterial, SkalakSolver, VelocityVerlet, WLCParameters,
    WLCSolver,
};
use glam::Vec3;

const N_SUBSTEPS: u32 = 10;
const DT_SEC: f32 = 1e-6;
const TEMPERATURE_K: f32 = 310.0;
const SEED: u32 = 0xDEADBEEF;
const MEMBRANE_DAMPING: f32 = 0.1;

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

/// Apply a small perturbation so forces are non-trivial.
fn perturb_mesh(mesh: &mut Mesh) {
    for (i, v) in mesh.vertices.iter_mut().enumerate() {
        let phase = (i as f32) * 0.7345;
        let scale = 0.05 * phase.sin();
        let normal = v.normal_vec3();
        let p = v.position_vec3() + normal * scale;
        v.position = p.to_array();
    }
}

/// CPU reference: mirror `PhysicsSolver::step` over N substeps but with
/// `compute_membrane_forces_pcg` for DPD and a fixed seed for parity.
fn run_cpu_reference(mesh: &mut Mesh, spectrin: &SpectrinNetwork, n_substeps: u32) -> Vec<Vec3> {
    let n_vertices = mesh.vertices.len();
    let mut velocities = vec![Vec3::ZERO; n_vertices];
    let mut forces = vec![Vec3::ZERO; n_vertices];

    let skalak_solver = SkalakSolver::new(mesh, SkalakMaterial::default());
    let wlc_solver = WLCSolver::new(WLCParameters::default());
    let dpd_solver = DPDSolver::new(DPDParameters::default());
    let mut integrator = VelocityVerlet::new();

    for step in 0..n_substeps {
        // 1. Half-step velocity using current forces.
        integrator.half_step_velocity(&mut velocities, &forces, DT_SEC);

        // 2. Update positions.
        let positions: Vec<Vec3> = mesh.vertices.iter().map(|v| v.position_vec3()).collect();
        let new_positions = integrator.step_position(&positions, &velocities, DT_SEC);
        for (i, p) in new_positions.iter().enumerate() {
            mesh.vertices[i].position = p.to_array();
        }

        // 3. Reset forces (no external forces in this test).
        for f in forces.iter_mut() {
            *f = Vec3::ZERO;
        }

        // 4. WLC forces — mapped from spectrin nodes to mesh vertices.
        let wlc_forces = wlc_solver.compute_forces(spectrin, &new_positions);
        for (node_idx, force) in wlc_forces {
            if let Some(node) = spectrin.nodes.get(node_idx) {
                let v_idx = node.mesh_vertex_idx as usize;
                if v_idx < n_vertices {
                    forces[v_idx] += force;
                }
            }
        }

        // 5. Skalak forces.
        let (skalak_forces, _energy) = skalak_solver.compute_forces(mesh);
        for (i, f) in skalak_forces.into_iter().enumerate() {
            forces[i] += f;
        }

        // 6. DPD forces (PCG path for parity).
        let dpd_forces = dpd_solver.compute_membrane_forces_pcg(
            &velocities,
            TEMPERATURE_K,
            step,
            SEED,
        );
        for (i, f) in dpd_forces.into_iter().enumerate() {
            forces[i] += f;
        }

        // 7. Membrane damping.
        for (i, vel) in velocities.iter().enumerate() {
            forces[i] -= *vel * MEMBRANE_DAMPING;
        }

        // 8. Second half-step velocity.
        integrator.half_step_velocity(&mut velocities, &forces, DT_SEC);
    }

    velocities
}

#[test]
fn physics_backend_cpu_gpu_parity() {
    let ctx = match ComputeContext::new_headless_blocking() {
        Ok(c) => c,
        Err(e) => {
            eprintln!("skipping GPU parity test: no adapter: {e}");
            return;
        }
    };

    // Build matching CPU and GPU initial state.
    let mut cpu_mesh = Mesh::generate_rbc(&default_params());
    perturb_mesh(&mut cpu_mesh);
    let mut gpu_mesh = Mesh::generate_rbc(&default_params());
    perturb_mesh(&mut gpu_mesh);
    let spectrin = SpectrinNetwork::generate(&cpu_mesh, &default_params());

    // === CPU run ===
    let cpu_velocities = run_cpu_reference(&mut cpu_mesh, &spectrin, N_SUBSTEPS);

    // === GPU run ===
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

    for _ in 0..N_SUBSTEPS {
        backend.step().expect("backend step");
    }
    let (gpu_positions, gpu_velocities) = backend.read_state().expect("read_state");

    // === Compare positions ===
    let n = cpu_mesh.vertices.len();
    assert_eq!(n, gpu_positions.len());
    let mut max_pos_err = 0.0f32;
    let mut max_pos_idx = 0usize;
    for i in 0..n {
        let cp = cpu_mesh.vertices[i].position_vec3();
        let gp = gpu_positions[i];
        let err = (cp - gp).length();
        if err > max_pos_err {
            max_pos_err = err;
            max_pos_idx = i;
        }
    }

    let mut max_vel_err = 0.0f32;
    let mut max_vel_idx = 0usize;
    for i in 0..n {
        let cv = cpu_velocities[i];
        let gv = gpu_velocities[i];
        let err = (cv - gv).length();
        if err > max_vel_err {
            max_vel_err = err;
            max_vel_idx = i;
        }
    }

    println!(
        "n_vertices = {}, n_substeps = {}",
        n, N_SUBSTEPS
    );
    println!(
        "  worst |Δpos| = {:.6e} at vertex {} (cpu = {:?} gpu = {:?})",
        max_pos_err, max_pos_idx,
        cpu_mesh.vertices[max_pos_idx].position_vec3(),
        gpu_positions[max_pos_idx]
    );
    println!(
        "  worst |Δvel| = {:.6e} at vertex {} (cpu = {:?} gpu = {:?})",
        max_vel_err, max_vel_idx,
        cpu_velocities[max_vel_idx], gpu_velocities[max_vel_idx]
    );

    // Tolerances: positions are tiny because dt is tiny, so allow 1e-5 μm
    // absolute. Velocities can be ~50 μm/s, so allow 1% relative or 1e-3 μm/s.
    assert!(
        max_pos_err < 1e-5,
        "position drift too large: {} at vertex {}",
        max_pos_err, max_pos_idx
    );
    let cpu_speed = cpu_velocities[max_vel_idx].length();
    let vel_pass = max_vel_err < 1e-3 || (cpu_speed > 1e-6 && max_vel_err / cpu_speed < 0.01);
    assert!(
        vel_pass,
        "velocity drift too large: {} at vertex {} (cpu speed {})",
        max_vel_err, max_vel_idx, cpu_speed
    );
}
