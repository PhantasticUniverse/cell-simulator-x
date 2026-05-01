//! Phase 11.3.D parity test: GPU vs CPU `VelocityVerlet::half_step_velocity`
//! + `step_position` for one combined step.
//!
//! The GPU kernel runs both operations in sequence; the CPU side does the
//! same. After one step, positions and velocities should match
//! bit-near-perfectly (within f32 noise from the speed/displacement clamps).
//!
//! Skipped if no GPU adapter is available.

use cell_simulator_x::compute::{run_verlet_step, ComputeContext, VerletConfig};
use cell_simulator_x::physics::VelocityVerlet;
use glam::Vec3;

#[test]
fn verlet_cpu_gpu_parity_one_step() {
    let ctx = match ComputeContext::new_headless_blocking() {
        Ok(c) => c,
        Err(e) => {
            eprintln!("skipping GPU parity test: no adapter: {e}");
            return;
        }
    };

    let n = 200;
    let mut positions: Vec<Vec3> = Vec::with_capacity(n);
    let mut velocities: Vec<Vec3> = Vec::with_capacity(n);
    let mut forces: Vec<Vec3> = Vec::with_capacity(n);
    for i in 0..n {
        let f = i as f32;
        positions.push(Vec3::new(
            (f * 0.317).sin() * 5.0,
            (f * 0.429).cos() * 5.0,
            (f * 0.591).sin() * 5.0,
        ));
        velocities.push(Vec3::new(
            (f * 0.713).sin() * 50.0,
            (f * 0.451).cos() * 50.0,
            (f * 0.234).sin() * 50.0,
        ));
        forces.push(Vec3::new(
            (f * 0.111).sin() * 10.0,
            (f * 0.222).cos() * 10.0,
            (f * 0.333).sin() * 10.0,
        ));
    }

    // CPU run.
    let mut cpu_pos = positions.clone();
    let mut cpu_vel = velocities.clone();
    let mut cpu_int = VelocityVerlet::new();
    let dt = 1e-5_f32;
    cpu_int.half_step_velocity(&mut cpu_vel, &forces, dt);
    let new_pos = cpu_int.step_position(&cpu_pos, &cpu_vel, dt);
    cpu_pos = new_pos;

    // GPU run.
    let mut gpu_pos = positions.clone();
    let mut gpu_vel = velocities.clone();
    let cfg = VerletConfig {
        dt_sec: dt,
        vertex_mass_pg: 1.0,
        max_velocity_um_per_sec: 1000.0,
        max_displacement_um: 0.1,
    };
    run_verlet_step(&ctx, &mut gpu_pos, &mut gpu_vel, &forces, &cfg).expect("GPU dispatch");

    let mut max_pos_err = 0.0f32;
    let mut max_vel_err = 0.0f32;
    for i in 0..n {
        let pe = (cpu_pos[i] - gpu_pos[i]).length();
        let ve = (cpu_vel[i] - gpu_vel[i]).length();
        max_pos_err = max_pos_err.max(pe);
        max_vel_err = max_vel_err.max(ve);
    }
    println!("max |Δpos| = {:.6e}, max |Δvel| = {:.6e}", max_pos_err, max_vel_err);

    assert!(max_pos_err < 1e-5, "max |Δpos| = {} exceeds 1e-5", max_pos_err);
    assert!(max_vel_err < 1e-3, "max |Δvel| = {} exceeds 1e-3", max_vel_err);
}
