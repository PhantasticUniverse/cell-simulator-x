//! Phase 11.3.C parity test: GPU vs CPU `DPDSolver::compute_membrane_forces_pcg`.
//!
//! Both backends use the same stateless PCG hash + Box-Muller, so given
//! identical (seed, step_count) the output is bit-identical.
//!
//! Skipped if no GPU adapter is available.

use cell_simulator_x::compute::{run_dpd_membrane_forces, ComputeContext};
use cell_simulator_x::physics::{DPDParameters, DPDSolver};
use glam::Vec3;

const TEMPERATURE_K: f32 = 310.15;
const SEED: u32 = 0xDEADBEEF;

#[test]
fn dpd_cpu_gpu_parity() {
    let ctx = match ComputeContext::new_headless_blocking() {
        Ok(c) => c,
        Err(e) => {
            eprintln!("skipping GPU parity test: no adapter: {e}");
            return;
        }
    };

    // Synthetic velocities — non-trivial mix.
    let n = 401;
    let mut velocities: Vec<Vec3> = Vec::with_capacity(n);
    for i in 0..n {
        let f = i as f32;
        velocities.push(Vec3::new(
            (f * 0.317).sin() * 100.0,
            (f * 0.429).cos() * 100.0,
            (f * 0.591).sin() * 100.0,
        ));
    }

    let params = DPDParameters::default();
    let solver = DPDSolver::new(params.clone());

    for step in [0u32, 1, 17, 1234] {
        let cpu_forces = solver.compute_membrane_forces_pcg(&velocities, TEMPERATURE_K, step, SEED);
        let gpu_forces = run_dpd_membrane_forces(
            &ctx,
            &velocities,
            &params,
            TEMPERATURE_K,
            step,
            SEED,
        )
        .expect("GPU dispatch");

        assert_eq!(cpu_forces.len(), gpu_forces.len());

        let mut max_abs_err = 0.0f32;
        let mut max_idx = 0usize;
        for i in 0..n {
            let diff = cpu_forces[i] - gpu_forces[i];
            let abs = diff.length();
            if abs > max_abs_err {
                max_abs_err = abs;
                max_idx = i;
            }
        }

        println!(
            "step={} worst |Δf|={:.6e} at vertex {} (cpu={:?} gpu={:?})",
            step, max_abs_err, max_idx, cpu_forces[max_idx], gpu_forces[max_idx]
        );
        // Bit-near-perfect — kernel and CPU code use identical math.
        assert!(
            max_abs_err < 1e-5,
            "step {}: worst |Δf|={} at vertex {}",
            step,
            max_abs_err,
            max_idx
        );
    }
}
