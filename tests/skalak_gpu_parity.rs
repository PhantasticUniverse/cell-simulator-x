//! Phase 11.3.A parity test: GPU vs CPU `SkalakSolver::compute_forces`.
//!
//! Builds a default RBC mesh, perturbs the vertices slightly so forces are
//! non-trivially nonzero, and asserts that GPU and CPU produce per-vertex
//! forces that match within 1e-3 μN absolute or 1% relative.
//!
//! Skipped if no GPU adapter is available.

use cell_simulator_x::compute::{run_skalak_forces, ComputeContext, SkalakBackendData};
use cell_simulator_x::config::GeometryParameters;
use cell_simulator_x::geometry::Mesh;
use cell_simulator_x::physics::{SkalakMaterial, SkalakSolver};
use glam::Vec3;

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

/// Build a mesh and perturb its vertices so forces are non-zero.
fn build_perturbed_mesh() -> Mesh {
    let mut mesh = Mesh::generate_rbc(&default_params());
    // Perturb each vertex slightly along its current normal so we get a
    // mix of positive and negative strains across the mesh. Use a simple
    // deterministic hash so the test is reproducible.
    for (i, v) in mesh.vertices.iter_mut().enumerate() {
        let phase = (i as f32) * 0.7345;
        let scale = 0.05 * phase.sin();
        let normal = v.normal_vec3();
        let p = v.position_vec3() + normal * scale;
        v.position = p.to_array();
    }
    mesh
}

#[test]
fn skalak_cpu_gpu_parity() {
    let ctx = match ComputeContext::new_headless_blocking() {
        Ok(c) => c,
        Err(e) => {
            eprintln!("skipping GPU parity test: no adapter: {e}");
            return;
        }
    };

    // Build matching CPU and GPU state.
    let cpu_mesh = build_perturbed_mesh();

    // The Skalak solver caches reference geometry from the *unperturbed*
    // mesh, so build the solver on a reference mesh and apply forces to
    // the perturbed positions. Keep both meshes' topology identical.
    let ref_mesh = Mesh::generate_rbc(&default_params());
    let material = SkalakMaterial::default();
    let solver = SkalakSolver::new(&ref_mesh, material.clone());

    // CPU run.
    let (cpu_forces, _energy) = solver.compute_forces(&cpu_mesh);

    // GPU run.
    let backend = SkalakBackendData::from_solver(&solver, cpu_mesh.vertices.len());
    let gpu_forces = run_skalak_forces(&ctx, &cpu_mesh, &backend, &material)
        .expect("GPU dispatch");

    // Compare per-vertex.
    assert_eq!(cpu_forces.len(), gpu_forces.len());

    let mut max_abs_err = 0.0f32;
    let mut max_abs_idx = 0usize;
    let mut max_rel_err = 0.0f32;
    let mut max_rel_idx = 0usize;
    let mut failures: Vec<(usize, Vec3, Vec3, f32)> = Vec::new();

    for (i, (cpu_f, gpu_f)) in cpu_forces.iter().zip(gpu_forces.iter()).enumerate() {
        let diff = *cpu_f - *gpu_f;
        let abs_err = diff.length();
        let cpu_mag = cpu_f.length();
        let rel_err = if cpu_mag > 1e-6 { abs_err / cpu_mag } else { 0.0 };

        if abs_err > max_abs_err {
            max_abs_err = abs_err;
            max_abs_idx = i;
        }
        if rel_err > max_rel_err {
            max_rel_err = rel_err;
            max_rel_idx = i;
        }

        // Tolerance: 1e-3 μN absolute OR 1% relative.
        let pass = abs_err < 1e-3 || rel_err < 0.01;
        if !pass {
            failures.push((i, *cpu_f, *gpu_f, abs_err));
        }
    }

    println!(
        "n_vertices = {}, n_elements = {}",
        backend.n_vertices, backend.n_elements
    );
    println!(
        "  worst absolute: vertex {} cpu = {:?} gpu = {:?} |err| = {:.6e}",
        max_abs_idx, cpu_forces[max_abs_idx], gpu_forces[max_abs_idx], max_abs_err
    );
    println!(
        "  worst relative: vertex {} cpu = {:?} gpu = {:?} rel-err = {:.4}%",
        max_rel_idx,
        cpu_forces[max_rel_idx],
        gpu_forces[max_rel_idx],
        max_rel_err * 100.0
    );

    assert!(
        failures.is_empty(),
        "{} vertices exceeded tolerance (first 5: {:?})",
        failures.len(),
        &failures[..failures.len().min(5)]
    );
}

#[test]
fn skalak_at_rest_is_small() {
    let ctx = match ComputeContext::new_headless_blocking() {
        Ok(c) => c,
        Err(e) => {
            eprintln!("skipping GPU parity test: no adapter: {e}");
            return;
        }
    };

    // Unperturbed mesh — forces should all be ~zero on both sides.
    let mesh = Mesh::generate_rbc(&default_params());
    let material = SkalakMaterial::default();
    let solver = SkalakSolver::new(&mesh, material.clone());
    let backend = SkalakBackendData::from_solver(&solver, mesh.vertices.len());
    let gpu_forces = run_skalak_forces(&ctx, &mesh, &backend, &material).unwrap();
    let max_force = gpu_forces.iter().map(|f| f.length()).fold(0f32, f32::max);
    assert!(max_force < 1.0, "GPU rest-state max force {} should be small", max_force);
}
