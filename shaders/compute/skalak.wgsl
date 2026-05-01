// Skalak membrane forces on GPU — Phase 11.3.A.
//
// Ports `SkalakSolver::compute_forces` from src/physics/membrane.rs to a
// two-kernel compute pipeline:
//
//   1. `skalak_per_element` — one thread per triangle. Reads the three
//      vertex positions, computes the three edge-spring forces and the
//      per-vertex area-preservation force (mirroring the CPU code line
//      244-313), writes 3 force vectors (one per local vertex) to a
//      per-element output buffer of length n_elements*9 f32.
//
//   2. `skalak_aggregate` — one thread per vertex. Reads its incident
//      elements via a precomputed CSR (offsets + packed (element_id,
//      local_vertex_idx)) and sums the per-element-vertex contributions
//      into the per-vertex force buffer.
//
// The CSR aggregation pattern replaces the CPU's serial scatter pattern
// (`forces[i0] += f`) with a deterministic, atomic-free GPU layout. WGSL
// does not have native f32 atomic add, so the two-pass shape is preferred
// over a per-element atomicAdd loop on i32 bit-reinterpretations.
//
// Mixed precision: f32 throughout, matching `SkalakSolver` (which is also
// f32 on the CPU side via glam's Vec3/Mat2). No conversion needed.

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
    n_elements: u32,
    n_vertices: u32,
    k_spring: f32,
    area_modulus: f32,
    area_force_scale: f32,
    _pad0: f32,
    _pad1: f32,
    _pad2: f32,
}

@group(0) @binding(0) var<storage, read>       vertex_pos:    array<f32>;          // n_vertices * 3
@group(0) @binding(1) var<storage, read>       elements:      array<Element>;     // n_elements
@group(0) @binding(2) var<storage, read_write> elem_forces:   array<f32>;         // n_elements * 9
@group(0) @binding(3) var<storage, read>       csr_offsets:   array<u32>;         // n_vertices + 1
@group(0) @binding(4) var<storage, read>       csr_data:      array<u32>;         // packed (element_id<<2 | local_idx)
@group(0) @binding(5) var<storage, read_write> vertex_forces: array<f32>;          // n_vertices * 3
@group(0) @binding(6) var<uniform>             u: U;

fn read_vertex(v: u32) -> vec3<f32> {
    let base = v * 3u;
    return vec3<f32>(vertex_pos[base + 0u], vertex_pos[base + 1u], vertex_pos[base + 2u]);
}

@compute @workgroup_size(64)
fn skalak_per_element(@builtin(global_invocation_id) gid: vec3<u32>) {
    let elem_id = gid.x;
    if (elem_id >= u.n_elements) { return; }

    let elem = elements[elem_id];
    let p0 = read_vertex(elem.v0);
    let p1 = read_vertex(elem.v1);
    let p2 = read_vertex(elem.v2);

    var f0 = vec3<f32>(0.0, 0.0, 0.0);
    var f1 = vec3<f32>(0.0, 0.0, 0.0);
    var f2 = vec3<f32>(0.0, 0.0, 0.0);

    // Edge 0-1.
    let e01 = p1 - p0;
    let len01 = length(e01);
    if (len01 > 1e-10) {
        let strain = (len01 - elem.ref_len_01) / max(elem.ref_len_01, 1e-10);
        let mag = u.k_spring * strain * elem.ref_len_01;
        let f = (e01 / len01) * mag;
        f0 = f0 + f;
        f1 = f1 - f;
    }

    // Edge 0-2.
    let e02 = p2 - p0;
    let len02 = length(e02);
    if (len02 > 1e-10) {
        let strain = (len02 - elem.ref_len_02) / max(elem.ref_len_02, 1e-10);
        let mag = u.k_spring * strain * elem.ref_len_02;
        let f = (e02 / len02) * mag;
        f0 = f0 + f;
        f2 = f2 - f;
    }

    // Edge 1-2.
    let e12 = p2 - p1;
    let len12 = length(e12);
    if (len12 > 1e-10) {
        let strain = (len12 - elem.ref_len_12) / max(elem.ref_len_12, 1e-10);
        let mag = u.k_spring * strain * elem.ref_len_12;
        let f = (e12 / len12) * mag;
        f1 = f1 + f;
        f2 = f2 - f;
    }

    // Area-preservation force (penalty if |area_strain| > 0.01).
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

    // Write per-element output. Layout: [f0_x f0_y f0_z f1_x f1_y f1_z f2_x f2_y f2_z].
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
fn skalak_aggregate(@builtin(global_invocation_id) gid: vec3<u32>) {
    let v = gid.x;
    if (v >= u.n_vertices) { return; }

    var sum = vec3<f32>(0.0, 0.0, 0.0);
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

    let out_base = v * 3u;
    vertex_forces[out_base + 0u] = sum.x;
    vertex_forces[out_base + 1u] = sum.y;
    vertex_forces[out_base + 2u] = sum.z;
}
