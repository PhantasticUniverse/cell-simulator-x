// Velocity-Verlet integrator on GPU — Phase 11.3.D.
//
// Two compute kernels mirroring `VelocityVerlet::half_step_velocity` and
// `VelocityVerlet::step_position` from src/physics/integrator.rs:
//
//   1. `verlet_half_step_velocity` — per-vertex update:
//        v += (F/m) * accel_scale * (dt/2)
//      with subsequent speed clamp to `max_velocity`.
//
//   2. `verlet_step_position` — per-vertex update:
//        x += clamp(v * dt, max_disp)
//
// Same `accel_scale = 10.0` magic constant from the CPU side, exposed as a
// uniform so callers and tests can verify it's set correctly.

struct U {
    n_vertices: u32,
    half_dt: f32,
    dt: f32,
    inv_mass: f32,
    accel_scale: f32,
    max_velocity: f32,
    max_displacement: f32,
    _pad0: f32,
}

@group(0) @binding(0) var<storage, read>       forces:     array<f32>;        // n_vertices * 3
@group(0) @binding(1) var<storage, read_write> velocities: array<f32>;        // n_vertices * 3
@group(0) @binding(2) var<storage, read_write> positions:  array<f32>;        // n_vertices * 3
@group(0) @binding(3) var<uniform>             u: U;

@compute @workgroup_size(64)
fn verlet_half_step_velocity(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    if (i >= u.n_vertices) { return; }

    let base = i * 3u;
    let f = vec3<f32>(forces[base + 0u], forces[base + 1u], forces[base + 2u]);
    var v = vec3<f32>(velocities[base + 0u], velocities[base + 1u], velocities[base + 2u]);

    let accel = f * u.inv_mass * u.accel_scale;
    v = v + accel * u.half_dt;

    let speed = length(v);
    if (speed > u.max_velocity) {
        v = v * (u.max_velocity / speed);
    }

    velocities[base + 0u] = v.x;
    velocities[base + 1u] = v.y;
    velocities[base + 2u] = v.z;
}

@compute @workgroup_size(64)
fn verlet_step_position(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    if (i >= u.n_vertices) { return; }

    let base = i * 3u;
    let v = vec3<f32>(velocities[base + 0u], velocities[base + 1u], velocities[base + 2u]);
    var p = vec3<f32>(positions[base + 0u], positions[base + 1u], positions[base + 2u]);

    var disp = v * u.dt;
    let mag = length(disp);
    if (mag > u.max_displacement) {
        disp = disp * (u.max_displacement / mag);
    }
    p = p + disp;

    positions[base + 0u] = p.x;
    positions[base + 1u] = p.y;
    positions[base + 2u] = p.z;
}
