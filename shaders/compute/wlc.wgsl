// WLC spectrin network forces on GPU — Phase 11.3.B.
//
// Ports `WLCSolver::compute_forces` from src/physics/wlc.rs to a two-kernel
// pipeline that mirrors the Skalak pattern:
//
//   1. `wlc_per_edge` — one thread per spectrin tetramer. Computes the
//      Marko-Siggia force from the two endpoints' positions and the
//      edge's contour length, writes 6 floats per edge (force on node_a
//      + force on node_b — opposite signs).
//
//   2. `wlc_aggregate` — one thread per spectrin node. Reads its
//      incident edges via CSR (offsets + packed (edge_id, side_bit)) and
//      sums into a per-node force buffer. Replaces the CPU-side HashMap
//      aggregation in `WLCSolver::compute_forces`.
//
// `side_bit = 0` means "this node is node_a of the edge" (read +force_a);
// `side_bit = 1` means "this node is node_b" (read +force_b).
//
// Mixed precision: f32 throughout (CPU side is also f32 via glam Vec3).

const SIDE_BIT_MASK: u32 = 1u;

struct Edge {
    node_a: u32,
    node_b: u32,
    contour_length_um: f32,
    _pad: f32,
}

struct U {
    n_edges: u32,
    n_nodes: u32,
    persistence_length_um: f32,
    kbt_uN_um: f32,             // 4.11e-9 by default
    max_relative_extension: f32, // 0.95 by default
    _pad0: f32,
    _pad1: f32,
    _pad2: f32,
}

@group(0) @binding(0) var<storage, read>       node_pos:    array<f32>;        // n_nodes * 3
@group(0) @binding(1) var<storage, read>       edges:       array<Edge>;
@group(0) @binding(2) var<storage, read_write> edge_forces: array<f32>;        // n_edges * 6 (a then b)
@group(0) @binding(3) var<storage, read>       csr_offsets: array<u32>;        // n_nodes + 1
@group(0) @binding(4) var<storage, read>       csr_data:    array<u32>;        // packed (edge_id<<1 | side_bit)
@group(0) @binding(5) var<storage, read_write> node_forces: array<f32>;        // n_nodes * 3
@group(0) @binding(6) var<uniform>             u: U;

fn read_node(n: u32) -> vec3<f32> {
    let base = n * 3u;
    return vec3<f32>(node_pos[base + 0u], node_pos[base + 1u], node_pos[base + 2u]);
}

fn marko_siggia_force(extension: f32, contour: f32) -> f32 {
    // Mirror src/physics/wlc.rs::marko_siggia_force.
    let max_ext = contour * u.max_relative_extension;
    let x = clamp(extension, 0.001 * contour, max_ext);
    let xi = x / contour;
    let one_minus_xi = 1.0 - xi;
    let term1 = 1.0 / (4.0 * one_minus_xi * one_minus_xi);
    return (u.kbt_uN_um / u.persistence_length_um) * (term1 - 0.25 + xi);
}

@compute @workgroup_size(64)
fn wlc_per_edge(@builtin(global_invocation_id) gid: vec3<u32>) {
    let edge_id = gid.x;
    if (edge_id >= u.n_edges) { return; }

    let edge = edges[edge_id];
    let pos_a = read_node(edge.node_a);
    let pos_b = read_node(edge.node_b);
    let diff = pos_b - pos_a;
    let len = length(diff);

    var force_a = vec3<f32>(0.0, 0.0, 0.0);
    var force_b = vec3<f32>(0.0, 0.0, 0.0);

    if (len >= 1e-10) {
        let dir = diff / len;
        let mag = marko_siggia_force(len, edge.contour_length_um);
        force_a = dir * mag;
        force_b = -force_a;
    }

    let base = edge_id * 6u;
    edge_forces[base + 0u] = force_a.x;
    edge_forces[base + 1u] = force_a.y;
    edge_forces[base + 2u] = force_a.z;
    edge_forces[base + 3u] = force_b.x;
    edge_forces[base + 4u] = force_b.y;
    edge_forces[base + 5u] = force_b.z;
}

@compute @workgroup_size(64)
fn wlc_aggregate(@builtin(global_invocation_id) gid: vec3<u32>) {
    let n = gid.x;
    if (n >= u.n_nodes) { return; }

    var sum = vec3<f32>(0.0, 0.0, 0.0);
    let start = csr_offsets[n];
    let end = csr_offsets[n + 1u];

    for (var i = start; i < end; i = i + 1u) {
        let packed = csr_data[i];
        let edge_id = packed >> 1u;
        let side = packed & SIDE_BIT_MASK;
        let base = edge_id * 6u + side * 3u;  // 0 → +force_a, 3 → +force_b
        sum.x = sum.x + edge_forces[base + 0u];
        sum.y = sum.y + edge_forces[base + 1u];
        sum.z = sum.z + edge_forces[base + 2u];
    }

    let out_base = n * 3u;
    node_forces[out_base + 0u] = sum.x;
    node_forces[out_base + 1u] = sum.y;
    node_forces[out_base + 2u] = sum.z;
}
