// Integrated physics substep on GPU — Phase 11.3.E.
//
// Six compute kernels share the same bind group layout and run in sequence
// to advance one substep, mirroring `PhysicsSolver::step` from
// src/physics/mod.rs:
//
//   1. `verlet_half_step_velocity` — v += (F/m)*accel_scale*(dt/2), clamp.
//   2. `verlet_step_position`      — x += clamp(v*dt, max_disp).
//   3. `skalak_per_element`        — per-element edge + area forces → scratch.
//   4. `skalak_init_from_baseline` — per-vertex: F = wlc_baseline + sum_skalak.
//   5. `dpd_add`                   — per-vertex: F += dpd(v, vel).
//   6. `damping_and_half_step`     — per-vertex: F -= vel*damping; v += (F/m)…
//
// WLC is baked into the per-vertex `wlc_baseline` buffer once at backend
// creation, since the CPU `WLCSolver::compute_forces` reads from
// `network.nodes` (which never updates as the mesh moves), making WLC
// forces a static contribution per substep — see docs/phase_11_3_notes.md.

const LOCAL_IDX_MASK: u32 = 3u;

struct Element {
    v0: u32,
    v1: u32,
    v2: u32,
    ref_area: f32,
    ref_len_01: f32,
    ref_len_02: f32,
    ref_len_12: f32,
    _pad: f32,
}

struct U {
    n_vertices: u32,
    n_elements: u32,
    half_dt: f32,
    dt: f32,
    inv_mass: f32,
    accel_scale: f32,
    max_velocity: f32,
    max_displacement: f32,
    k_spring: f32,
    area_modulus: f32,
    area_force_scale: f32,
    membrane_damping: f32,
    gamma: f32,
    sigma: f32,
    visc_scale: f32,
    rand_scale: f32,
    step_count: u32,
    seed: u32,
    _pad0: f32,
    _pad1: f32,
}

@group(0) @binding(0)  var<storage, read_write> positions:     array<f32>;        // n_vertices * 3
@group(0) @binding(1)  var<storage, read_write> velocities:    array<f32>;        // n_vertices * 3
@group(0) @binding(2)  var<storage, read_write> forces:        array<f32>;        // n_vertices * 3
@group(0) @binding(3)  var<storage, read>       wlc_baseline:  array<f32>;        // n_vertices * 3
@group(0) @binding(4)  var<storage, read>       elements:      array<Element>;    // n_elements
@group(0) @binding(5)  var<storage, read>       csr_offsets:   array<u32>;        // n_vertices + 1
@group(0) @binding(6)  var<storage, read>       csr_data:      array<u32>;
@group(0) @binding(7)  var<storage, read_write> elem_forces:   array<f32>;        // n_elements * 9
@group(0) @binding(8)  var<uniform>             u: U;

// === Verlet kernels =================================================

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

// === Skalak ========================================================

@compute @workgroup_size(64)
fn skalak_per_element(@builtin(global_invocation_id) gid: vec3<u32>) {
    let elem_id = gid.x;
    if (elem_id >= u.n_elements) { return; }

    let elem = elements[elem_id];
    let p0_b = elem.v0 * 3u;
    let p1_b = elem.v1 * 3u;
    let p2_b = elem.v2 * 3u;
    let p0 = vec3<f32>(positions[p0_b + 0u], positions[p0_b + 1u], positions[p0_b + 2u]);
    let p1 = vec3<f32>(positions[p1_b + 0u], positions[p1_b + 1u], positions[p1_b + 2u]);
    let p2 = vec3<f32>(positions[p2_b + 0u], positions[p2_b + 1u], positions[p2_b + 2u]);

    var f0 = vec3<f32>(0.0, 0.0, 0.0);
    var f1 = vec3<f32>(0.0, 0.0, 0.0);
    var f2 = vec3<f32>(0.0, 0.0, 0.0);

    let e01 = p1 - p0;
    let len01 = length(e01);
    if (len01 > 1e-10) {
        let strain = (len01 - elem.ref_len_01) / max(elem.ref_len_01, 1e-10);
        let mag = u.k_spring * strain * elem.ref_len_01;
        let f = (e01 / len01) * mag;
        f0 = f0 + f;
        f1 = f1 - f;
    }

    let e02 = p2 - p0;
    let len02 = length(e02);
    if (len02 > 1e-10) {
        let strain = (len02 - elem.ref_len_02) / max(elem.ref_len_02, 1e-10);
        let mag = u.k_spring * strain * elem.ref_len_02;
        let f = (e02 / len02) * mag;
        f0 = f0 + f;
        f2 = f2 - f;
    }

    let e12 = p2 - p1;
    let len12 = length(e12);
    if (len12 > 1e-10) {
        let strain = (len12 - elem.ref_len_12) / max(elem.ref_len_12, 1e-10);
        let mag = u.k_spring * strain * elem.ref_len_12;
        let f = (e12 / len12) * mag;
        f1 = f1 + f;
        f2 = f2 - f;
    }

    let cross_v = cross(e01, e02);
    let current_area = 0.5 * length(cross_v);
    let area_strain = (current_area - elem.ref_area) / max(elem.ref_area, 1e-10);
    if (abs(area_strain) > 0.01) {
        let center = (p0 + p1 + p2) / 3.0;
        let scale = u.area_modulus * area_strain * u.area_force_scale;
        let to_c0 = center - p0;
        let d0 = length(to_c0);
        if (d0 > 1e-10) { f0 = f0 + (to_c0 / d0) * scale; }
        let to_c1 = center - p1;
        let d1 = length(to_c1);
        if (d1 > 1e-10) { f1 = f1 + (to_c1 / d1) * scale; }
        let to_c2 = center - p2;
        let d2 = length(to_c2);
        if (d2 > 1e-10) { f2 = f2 + (to_c2 / d2) * scale; }
    }

    let base = elem_id * 9u;
    elem_forces[base + 0u] = f0.x;
    elem_forces[base + 1u] = f0.y;
    elem_forces[base + 2u] = f0.z;
    elem_forces[base + 3u] = f1.x;
    elem_forces[base + 4u] = f1.y;
    elem_forces[base + 5u] = f1.z;
    elem_forces[base + 6u] = f2.x;
    elem_forces[base + 7u] = f2.y;
    elem_forces[base + 8u] = f2.z;
}

@compute @workgroup_size(64)
fn skalak_init_from_baseline(@builtin(global_invocation_id) gid: vec3<u32>) {
    let v = gid.x;
    if (v >= u.n_vertices) { return; }

    // Start from WLC baseline (constant per substep — see file header).
    let bb = v * 3u;
    var sum = vec3<f32>(wlc_baseline[bb + 0u], wlc_baseline[bb + 1u], wlc_baseline[bb + 2u]);

    let start = csr_offsets[v];
    let end = csr_offsets[v + 1u];
    for (var i = start; i < end; i = i + 1u) {
        let packed = csr_data[i];
        let elem_id = packed >> 2u;
        let local_idx = packed & LOCAL_IDX_MASK;
        let base = elem_id * 9u + local_idx * 3u;
        sum.x = sum.x + elem_forces[base + 0u];
        sum.y = sum.y + elem_forces[base + 1u];
        sum.z = sum.z + elem_forces[base + 2u];
    }

    let out = v * 3u;
    forces[out + 0u] = sum.x;
    forces[out + 1u] = sum.y;
    forces[out + 2u] = sum.z;
}

// === DPD (PCG hash) — ADD semantics =================================

const TAU: f32 = 6.28318530718;

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
fn dpd_add(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    if (i >= u.n_vertices) { return; }
    let base = i * 3u;
    let vel = vec3<f32>(velocities[base + 0u], velocities[base + 1u], velocities[base + 2u]);

    let f_d = -vel * u.gamma * u.visc_scale;

    let key = (i * 0x9E3779B1u) ^ (u.step_count * 0x85EBCA77u) ^ u.seed;
    let gx = pcg_gaussian(key, key ^ 0x00000001u);
    let gy = pcg_gaussian(key ^ 0x00010000u, key ^ 0x00010001u);
    let gz = pcg_gaussian(key ^ 0x00020000u, key ^ 0x00020001u);
    let f_r = vec3<f32>(gx, gy, gz) * u.sigma * u.rand_scale;

    let total = f_d + f_r;
    forces[base + 0u] = forces[base + 0u] + total.x;
    forces[base + 1u] = forces[base + 1u] + total.y;
    forces[base + 2u] = forces[base + 2u] + total.z;
}

// === Damping + final half-step (combined) ===========================

@compute @workgroup_size(64)
fn damping_and_half_step(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    if (i >= u.n_vertices) { return; }
    let base = i * 3u;

    var v = vec3<f32>(velocities[base + 0u], velocities[base + 1u], velocities[base + 2u]);
    var f = vec3<f32>(forces[base + 0u], forces[base + 1u], forces[base + 2u]);

    // Membrane damping.
    f = f - v * u.membrane_damping;

    // Half-step velocity (using post-damping forces).
    let accel = f * u.inv_mass * u.accel_scale;
    v = v + accel * u.half_dt;
    let speed = length(v);
    if (speed > u.max_velocity) {
        v = v * (u.max_velocity / speed);
    }

    velocities[base + 0u] = v.x;
    velocities[base + 1u] = v.y;
    velocities[base + 2u] = v.z;
    // Write back the damped force so the host can read it back if needed.
    forces[base + 0u] = f.x;
    forces[base + 1u] = f.y;
    forces[base + 2u] = f.z;
}
