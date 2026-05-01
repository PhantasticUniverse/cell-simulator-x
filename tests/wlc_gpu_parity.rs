//! Phase 11.3.B parity test: GPU vs CPU `WLCSolver::compute_forces`.
//!
//! Builds a default RBC mesh + spectrin network, perturbs spectrin node
//! positions slightly, then asserts that GPU per-node forces match CPU
//! within 1e-3 μN absolute or 1% relative.
//!
//! Skipped if no GPU adapter is available.

use cell_simulator_x::compute::{run_wlc_forces, ComputeContext, WlcBackendData};
use cell_simulator_x::config::GeometryParameters;
use cell_simulator_x::geometry::{Mesh, SpectrinNetwork};
use cell_simulator_x::physics::{WLCParameters, WLCSolver};
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

fn build_perturbed_network() -> SpectrinNetwork {
    let mesh = Mesh::generate_rbc(&default_params());
    let mut net = SpectrinNetwork::generate(&mesh, &default_params());
    // Perturb each node deterministically so spectrin tetramers are
    // stretched a non-trivial amount.
    for (i, n) in net.nodes.iter_mut().enumerate() {
        let phase_x = (i as f32) * 0.314;
        let phase_y = (i as f32) * 0.6739;
        let phase_z = (i as f32) * 1.0731;
        let dx = 0.02 * phase_x.sin();
        let dy = 0.02 * phase_y.sin();
        let dz = 0.02 * phase_z.sin();
        let p = n.position_vec3() + Vec3::new(dx, dy, dz);
        n.position = p.to_array();
    }
    net
}

#[test]
fn wlc_cpu_gpu_parity() {
    let ctx = match ComputeContext::new_headless_blocking() {
        Ok(c) => c,
        Err(e) => {
            eprintln!("skipping GPU parity test: no adapter: {e}");
            return;
        }
    };

    let net = build_perturbed_network();
    let params = WLCParameters::default();
    let solver = WLCSolver::new(params.clone());

    // CPU: returns sparse Vec<(idx, Vec3)>; densify for comparison.
    let cpu_sparse = solver.compute_forces(&net, &[]);
    let mut cpu_dense = vec![Vec3::ZERO; net.nodes.len()];
    for (idx, f) in cpu_sparse {
        cpu_dense[idx] = f;
    }

    let backend = WlcBackendData::from_network(&net);
    let gpu_dense = run_wlc_forces(&ctx, &net, &backend, &params).expect("GPU dispatch");

    assert_eq!(cpu_dense.len(), gpu_dense.len());

    let mut max_abs_err = 0.0f32;
    let mut max_abs_idx = 0usize;
    let mut max_rel_err = 0.0f32;
    let mut max_rel_idx = 0usize;
    let mut failures: Vec<(usize, Vec3, Vec3, f32)> = Vec::new();

    for (i, (cpu_f, gpu_f)) in cpu_dense.iter().zip(gpu_dense.iter()).enumerate() {
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

        let pass = abs_err < 1e-3 || rel_err < 0.01;
        if !pass {
            failures.push((i, *cpu_f, *gpu_f, abs_err));
        }
    }

    println!(
        "n_nodes = {}, n_edges = {}",
        backend.n_nodes, backend.n_edges
    );
    println!(
        "  worst absolute: node {} cpu = {:?} gpu = {:?} |err| = {:.6e}",
        max_abs_idx, cpu_dense[max_abs_idx], gpu_dense[max_abs_idx], max_abs_err
    );
    println!(
        "  worst relative: node {} cpu = {:?} gpu = {:?} rel-err = {:.4}%",
        max_rel_idx,
        cpu_dense[max_rel_idx],
        gpu_dense[max_rel_idx],
        max_rel_err * 100.0
    );

    assert!(
        failures.is_empty(),
        "{} nodes exceeded tolerance (first 5: {:?})",
        failures.len(),
        &failures[..failures.len().min(5)]
    );
}
