// DPD membrane forces on GPU — Phase 11.3.C.
//
// Ports `DPDSolver::compute_membrane_forces_pcg` (Rust) to a single
// per-vertex compute kernel. Each vertex computes:
//
//   F_d = -γ * 0.001 * v
//   F_r = σ * 0.0001 * (gx, gy, gz)
//
// where σ = sqrt(2*γ*kbt) and gx/gy/gz are independent Gaussians produced
// from a stateless PCG hash keyed on (vertex_id, step_count, seed, axis).
// The hash + Box-Muller implementation matches `pcg_hash_u32` /
// `pcg_uniform` / `pcg_gaussian` in src/physics/dpd.rs verbatim, so CPU
// and GPU produce bit-identical noise sequences.

const TAU: f32 = 6.28318530718;

struct U {
    n_vertices: u32,
    step_count: u32,
    seed: u32,
    gamma: f32,
    sigma: f32,
    visc_scale: f32,        // 0.001 — matches CPU constant
    rand_scale: f32,        // 0.0001 — matches CPU constant
    _pad0: f32,
}

@group(0) @binding(0) var<storage, read>       velocities: array<f32>;        // n_vertices * 3
@group(0) @binding(1) var<storage, read_write> forces:     array<f32>;        // n_vertices * 3
@group(0) @binding(2) var<uniform>             u: U;

fn pcg_hash_u32(seed: u32) -> u32 {
    let state = seed * 747796405u + 2891336453u;
    let shift = ((state >> 28u) + 4u) & 31u;
    let word = ((state >> shift) ^ state) * 277803737u;
    return (word >> 22u) ^ word;
}

fn pcg_uniform(seed: u32) -> f32 {
    return f32(pcg_hash_u32(seed)) / 4294967296.0;
}

fn pcg_gaussian(seed1: u32, seed2: u32) -> f32 {
    let u1 = max(pcg_uniform(seed1), 1e-30);
    let u2 = pcg_uniform(seed2);
    let r = sqrt(-2.0 * log(u1));
    let theta = TAU * u2;
    return r * cos(theta);
}

@compute @workgroup_size(64)
fn dpd_membrane(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    if (i >= u.n_vertices) { return; }

    let base = i * 3u;
    let vel = vec3<f32>(velocities[base + 0u], velocities[base + 1u], velocities[base + 2u]);

    // Viscous force.
    let f_d = -vel * u.gamma * u.visc_scale;

    // Random force from PCG hash.
    let key = (i * 0x9E3779B1u) ^ (u.step_count * 0x85EBCA77u) ^ u.seed;
    let gx = pcg_gaussian(key, key ^ 0x00000001u);
    let gy = pcg_gaussian(key ^ 0x00010000u, key ^ 0x00010001u);
    let gz = pcg_gaussian(key ^ 0x00020000u, key ^ 0x00020001u);
    let f_r = vec3<f32>(gx, gy, gz) * u.sigma * u.rand_scale;

    let total = f_d + f_r;
    forces[base + 0u] = total.x;
    forces[base + 1u] = total.y;
    forces[base + 2u] = total.z;
}
