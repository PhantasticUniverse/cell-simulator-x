//! GPU physics kernels (Phase 11.3).
//!
//! Currently ports `SkalakSolver::compute_forces` to a two-kernel compute
//! pipeline: per-element force compute + per-vertex CSR aggregation. The
//! WLC, DPD, and Velocity-Verlet ports follow in 11.3.B/C/D.

use std::sync::mpsc;

use anyhow::{Context as _, Result};
use bytemuck::{Pod, Zeroable};
use glam::Vec3;
use wgpu::util::DeviceExt;

use crate::geometry::{Mesh, SpectrinNetwork};
use crate::physics::{DPDParameters, SkalakMaterial, SkalakSolver, WLCParameters};

use super::ComputeContext;

const KERNEL_WGSL: &str = include_str!("../../shaders/compute/skalak.wgsl");
const KERNEL_WLC_WGSL: &str = include_str!("../../shaders/compute/wlc.wgsl");
const KERNEL_DPD_WGSL: &str = include_str!("../../shaders/compute/dpd.wgsl");
const KERNEL_INTEGRATOR_WGSL: &str =
    include_str!("../../shaders/compute/integrator.wgsl");
const KERNEL_STEP_WGSL: &str = include_str!("../../shaders/compute/physics_step.wgsl");

#[repr(C)]
#[derive(Copy, Clone, Debug, Pod, Zeroable)]
struct ElementGpu {
    v0: u32,
    v1: u32,
    v2: u32,
    ref_area: f32,
    ref_len_01: f32,
    ref_len_02: f32,
    ref_len_12: f32,
    _pad: f32,
}

#[repr(C)]
#[derive(Copy, Clone, Debug, Pod, Zeroable)]
struct SkalakUniforms {
    n_elements: u32,
    n_vertices: u32,
    k_spring: f32,
    area_modulus: f32,
    area_force_scale: f32,
    _pad0: f32,
    _pad1: f32,
    _pad2: f32,
}

/// Precomputed per-mesh data needed to dispatch Skalak forces on the GPU.
///
/// Built once on the host from the mesh + solver and reused across many
/// dispatches (each substep re-uploads only the vertex positions).
#[derive(Debug, Clone)]
pub struct SkalakBackendData {
    /// Per-element packed data.
    pub elements: Vec<ElementGpuPub>,
    /// CSR offsets — `csr_offsets[v+1] - csr_offsets[v]` is the number of
    /// triangles incident at vertex `v`.
    pub csr_offsets: Vec<u32>,
    /// CSR data — packed `(element_id << 2) | local_idx`. `local_idx` ∈ {0,1,2}
    /// indicates which corner of the triangle is `v`.
    pub csr_data: Vec<u32>,
    /// Number of vertices.
    pub n_vertices: usize,
    /// Number of elements.
    pub n_elements: usize,
}

/// Public re-export of the GPU element struct so callers can keep their
/// own backend data without depending on the (private) crate-internal Pod.
#[repr(C)]
#[derive(Copy, Clone, Debug, Pod, Zeroable)]
pub struct ElementGpuPub {
    pub v0: u32,
    pub v1: u32,
    pub v2: u32,
    pub ref_area: f32,
    pub ref_len_01: f32,
    pub ref_len_02: f32,
    pub ref_len_12: f32,
    pub _pad: f32,
}

impl SkalakBackendData {
    /// Build backend data from a mesh + Skalak solver.
    pub fn from_solver(solver: &SkalakSolver, n_vertices: usize) -> Self {
        let mut elements = Vec::with_capacity(solver.elements.len());
        let mut incidence: Vec<Vec<u32>> = vec![Vec::new(); n_vertices];

        for (eid, e) in solver.elements.iter().enumerate() {
            // Mirror the CPU's reference-edge length computation.
            let ref_e0 = e.reference_edges[0];
            let ref_e1 = e.reference_edges[1];
            let ref_e2 = e.reference_positions_2d[2] - e.reference_positions_2d[1];
            let ref_len_01 = ref_e0.length();
            let ref_len_02 = ref_e1.length();
            let ref_len_12 = ref_e2.length();

            let v0 = e.vertex_indices[0];
            let v1 = e.vertex_indices[1];
            let v2 = e.vertex_indices[2];

            elements.push(ElementGpuPub {
                v0,
                v1,
                v2,
                ref_area: e.reference_area_um2,
                ref_len_01,
                ref_len_02,
                ref_len_12,
                _pad: 0.0,
            });

            // Build CSR incidence: vertex → list of (element_id, local_idx).
            let eid = eid as u32;
            incidence[v0 as usize].push((eid << 2) | 0);
            incidence[v1 as usize].push((eid << 2) | 1);
            incidence[v2 as usize].push((eid << 2) | 2);
        }

        // Flatten incidence to CSR.
        let mut csr_offsets = Vec::with_capacity(n_vertices + 1);
        let mut csr_data = Vec::new();
        csr_offsets.push(0);
        for inc in &incidence {
            csr_data.extend_from_slice(inc);
            csr_offsets.push(csr_data.len() as u32);
        }

        Self {
            n_elements: elements.len(),
            elements,
            csr_offsets,
            csr_data,
            n_vertices,
        }
    }
}

/// Run Skalak forces on the GPU.
///
/// One-shot dispatch: uploads vertex positions, dispatches both kernels,
/// reads back per-vertex forces. Pipeline + bind group are recreated per
/// call. A persistent `SkalakBackend` will land in 11.3.D once the kernel
/// surface is frozen.
///
/// Returns the per-vertex force vector (μN), in the same order as
/// `mesh.vertices`.
pub fn run_skalak_forces(
    ctx: &ComputeContext,
    mesh: &Mesh,
    backend: &SkalakBackendData,
    material: &SkalakMaterial,
) -> Result<Vec<Vec3>> {
    let n_vertices = backend.n_vertices;
    let n_elements = backend.n_elements;
    if n_vertices != mesh.vertices.len() {
        anyhow::bail!(
            "vertex count mismatch: backend has {}, mesh has {}",
            n_vertices,
            mesh.vertices.len()
        );
    }

    // Upload vertex positions (SoA: x,y,z per vertex, interleaved as
    // 3-tuple but kept scalar to avoid the WGSL vec3 alignment trap).
    let mut pos_f32 = vec![0f32; n_vertices * 3];
    for (i, v) in mesh.vertices.iter().enumerate() {
        pos_f32[i * 3 + 0] = v.position[0];
        pos_f32[i * 3 + 1] = v.position[1];
        pos_f32[i * 3 + 2] = v.position[2];
    }

    let buf_pos = ctx.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("skalak vertex_pos"),
        contents: bytemuck::cast_slice(&pos_f32),
        usage: wgpu::BufferUsages::STORAGE,
    });

    let elements_pod: Vec<ElementGpu> = backend
        .elements
        .iter()
        .map(|e| ElementGpu {
            v0: e.v0,
            v1: e.v1,
            v2: e.v2,
            ref_area: e.ref_area,
            ref_len_01: e.ref_len_01,
            ref_len_02: e.ref_len_02,
            ref_len_12: e.ref_len_12,
            _pad: 0.0,
        })
        .collect();

    let buf_elements = ctx.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("skalak elements"),
        contents: bytemuck::cast_slice(&elements_pod),
        usage: wgpu::BufferUsages::STORAGE,
    });

    let elem_forces_bytes = (n_elements * 9 * std::mem::size_of::<f32>()) as u64;
    let buf_elem_forces = ctx.device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("skalak elem_forces"),
        size: elem_forces_bytes.max(16),
        usage: wgpu::BufferUsages::STORAGE,
        mapped_at_creation: false,
    });

    let buf_csr_offsets = ctx.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("skalak csr_offsets"),
        contents: bytemuck::cast_slice(&backend.csr_offsets),
        usage: wgpu::BufferUsages::STORAGE,
    });
    let buf_csr_data = ctx.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("skalak csr_data"),
        contents: bytemuck::cast_slice(&backend.csr_data),
        usage: wgpu::BufferUsages::STORAGE,
    });

    let vertex_forces_bytes = (n_vertices * 3 * std::mem::size_of::<f32>()) as u64;
    let buf_vertex_forces = ctx.device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("skalak vertex_forces"),
        size: vertex_forces_bytes,
        usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
        mapped_at_creation: false,
    });

    // CPU spring stiffness mirrors `SkalakSolver::compute_forces` line 242.
    let k_spring = material.shear_modulus_uN_per_m * 1.732;
    // Area force scale mirrors `SkalakSolver::compute_forces` line 303.
    let area_force_scale = 0.001f32;

    let uniforms = SkalakUniforms {
        n_elements: n_elements as u32,
        n_vertices: n_vertices as u32,
        k_spring,
        area_modulus: material.area_modulus_uN_per_m,
        area_force_scale,
        _pad0: 0.0,
        _pad1: 0.0,
        _pad2: 0.0,
    };
    let buf_uniforms = ctx.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("skalak uniforms"),
        contents: bytemuck::bytes_of(&uniforms),
        usage: wgpu::BufferUsages::UNIFORM,
    });

    let buf_staging = ctx.device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("skalak staging"),
        size: vertex_forces_bytes,
        usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
        mapped_at_creation: false,
    });

    // Bind group layout — six storage buffers + one uniform.
    let storage_entry = |binding: u32, read_only: bool| wgpu::BindGroupLayoutEntry {
        binding,
        visibility: wgpu::ShaderStages::COMPUTE,
        ty: wgpu::BindingType::Buffer {
            ty: wgpu::BufferBindingType::Storage { read_only },
            has_dynamic_offset: false,
            min_binding_size: None,
        },
        count: None,
    };

    let bgl = ctx.device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
        label: Some("skalak bgl"),
        entries: &[
            storage_entry(0, true),  // vertex_pos
            storage_entry(1, true),  // elements
            storage_entry(2, false), // elem_forces (read_write)
            storage_entry(3, true),  // csr_offsets
            storage_entry(4, true),  // csr_data
            storage_entry(5, false), // vertex_forces (read_write)
            wgpu::BindGroupLayoutEntry {
                binding: 6,
                visibility: wgpu::ShaderStages::COMPUTE,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Uniform,
                    has_dynamic_offset: false,
                    min_binding_size: None,
                },
                count: None,
            },
        ],
    });

    let bg = ctx.device.create_bind_group(&wgpu::BindGroupDescriptor {
        label: Some("skalak bg"),
        layout: &bgl,
        entries: &[
            wgpu::BindGroupEntry { binding: 0, resource: buf_pos.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 1, resource: buf_elements.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 2, resource: buf_elem_forces.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 3, resource: buf_csr_offsets.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 4, resource: buf_csr_data.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 5, resource: buf_vertex_forces.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 6, resource: buf_uniforms.as_entire_binding() },
        ],
    });

    let shader = ctx.device.create_shader_module(wgpu::ShaderModuleDescriptor {
        label: Some("skalak shader"),
        source: wgpu::ShaderSource::Wgsl(KERNEL_WGSL.into()),
    });

    let pipeline_layout = ctx.device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
        label: Some("skalak pl"),
        bind_group_layouts: &[&bgl],
        push_constant_ranges: &[],
    });

    let pipe_per_element = ctx.device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
        label: Some("skalak per_element pipeline"),
        layout: Some(&pipeline_layout),
        module: &shader,
        entry_point: "skalak_per_element",
        compilation_options: Default::default(),
    });
    let pipe_aggregate = ctx.device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
        label: Some("skalak aggregate pipeline"),
        layout: Some(&pipeline_layout),
        module: &shader,
        entry_point: "skalak_aggregate",
        compilation_options: Default::default(),
    });

    let mut encoder = ctx.device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
        label: Some("skalak encoder"),
    });

    {
        let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
            label: Some("skalak pass"),
            timestamp_writes: None,
        });
        pass.set_bind_group(0, &bg, &[]);

        // Per-element kernel.
        pass.set_pipeline(&pipe_per_element);
        let groups = ((n_elements as u32) + 63) / 64;
        pass.dispatch_workgroups(groups.max(1), 1, 1);

        // Per-vertex aggregation kernel.
        pass.set_pipeline(&pipe_aggregate);
        let groups = ((n_vertices as u32) + 63) / 64;
        pass.dispatch_workgroups(groups.max(1), 1, 1);
    }

    encoder.copy_buffer_to_buffer(&buf_vertex_forces, 0, &buf_staging, 0, vertex_forces_bytes);
    ctx.queue.submit(std::iter::once(encoder.finish()));

    let slice = buf_staging.slice(..);
    let (tx, rx) = mpsc::channel();
    slice.map_async(wgpu::MapMode::Read, move |result| {
        let _ = tx.send(result);
    });
    ctx.device.poll(wgpu::Maintain::Wait);
    rx.recv()
        .context("skalak: map_async channel closed")?
        .map_err(|e| anyhow::anyhow!("skalak: map_async failed: {e:?}"))?;

    let mapped = slice.get_mapped_range();
    let result_f32: &[f32] = bytemuck::cast_slice(&mapped);

    let mut out = Vec::with_capacity(n_vertices);
    for i in 0..n_vertices {
        out.push(Vec3::new(
            result_f32[i * 3 + 0],
            result_f32[i * 3 + 1],
            result_f32[i * 3 + 2],
        ));
    }

    drop(mapped);
    buf_staging.unmap();

    Ok(out)
}

// =============================================================================
// Phase 11.3.B — WLC spectrin forces
// =============================================================================

#[repr(C)]
#[derive(Copy, Clone, Debug, Pod, Zeroable)]
struct EdgeGpu {
    node_a: u32,
    node_b: u32,
    contour_length_um: f32,
    _pad: f32,
}

#[repr(C)]
#[derive(Copy, Clone, Debug, Pod, Zeroable)]
struct WlcUniforms {
    n_edges: u32,
    n_nodes: u32,
    persistence_length_um: f32,
    kbt_uN_um: f32,
    max_relative_extension: f32,
    _pad0: f32,
    _pad1: f32,
    _pad2: f32,
}

/// Precomputed per-network data for WLC dispatch (analogue of
/// `SkalakBackendData`).
#[derive(Debug, Clone)]
pub struct WlcBackendData {
    /// Per-edge data uploaded to the GPU.
    pub edges: Vec<EdgeGpuPub>,
    /// CSR offsets indexed by node.
    pub csr_offsets: Vec<u32>,
    /// Packed `(edge_id << 1) | side_bit` (0 = node_a, 1 = node_b).
    pub csr_data: Vec<u32>,
    /// Number of nodes.
    pub n_nodes: usize,
    /// Number of edges.
    pub n_edges: usize,
}

#[repr(C)]
#[derive(Copy, Clone, Debug, Pod, Zeroable)]
pub struct EdgeGpuPub {
    pub node_a: u32,
    pub node_b: u32,
    pub contour_length_um: f32,
    pub _pad: f32,
}

impl WlcBackendData {
    pub fn from_network(network: &SpectrinNetwork) -> Self {
        let n_nodes = network.nodes.len();
        let n_edges = network.edges.len();

        let mut edges = Vec::with_capacity(n_edges);
        let mut incidence: Vec<Vec<u32>> = vec![Vec::new(); n_nodes];

        for (eid, e) in network.edges.iter().enumerate() {
            edges.push(EdgeGpuPub {
                node_a: e.node_a,
                node_b: e.node_b,
                contour_length_um: e.contour_length_um,
                _pad: 0.0,
            });
            let eid = eid as u32;
            incidence[e.node_a as usize].push((eid << 1) | 0);
            incidence[e.node_b as usize].push((eid << 1) | 1);
        }

        let mut csr_offsets = Vec::with_capacity(n_nodes + 1);
        let mut csr_data = Vec::new();
        csr_offsets.push(0);
        for inc in &incidence {
            csr_data.extend_from_slice(inc);
            csr_offsets.push(csr_data.len() as u32);
        }

        Self {
            edges,
            csr_offsets,
            csr_data,
            n_nodes,
            n_edges,
        }
    }
}

/// Run WLC spectrin forces on the GPU. Returns dense per-node forces in
/// μN. Unlike the CPU counterpart, which returns `Vec<(idx, Vec3)>` only
/// for nodes touched by some edge, this returns a dense vector of length
/// `network.nodes.len()` (zero for unconnected nodes).
pub fn run_wlc_forces(
    ctx: &ComputeContext,
    network: &SpectrinNetwork,
    backend: &WlcBackendData,
    params: &WLCParameters,
) -> Result<Vec<Vec3>> {
    let n_nodes = backend.n_nodes;
    let n_edges = backend.n_edges;
    if n_nodes != network.nodes.len() {
        anyhow::bail!(
            "node count mismatch: backend has {}, network has {}",
            n_nodes,
            network.nodes.len()
        );
    }

    // Upload node positions.
    let mut pos_f32 = vec![0f32; n_nodes * 3];
    for (i, n) in network.nodes.iter().enumerate() {
        pos_f32[i * 3 + 0] = n.position[0];
        pos_f32[i * 3 + 1] = n.position[1];
        pos_f32[i * 3 + 2] = n.position[2];
    }
    let buf_pos = ctx.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("wlc node_pos"),
        contents: bytemuck::cast_slice(&pos_f32),
        usage: wgpu::BufferUsages::STORAGE,
    });

    let edges_pod: Vec<EdgeGpu> = backend
        .edges
        .iter()
        .map(|e| EdgeGpu {
            node_a: e.node_a,
            node_b: e.node_b,
            contour_length_um: e.contour_length_um,
            _pad: 0.0,
        })
        .collect();
    let buf_edges = ctx.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("wlc edges"),
        contents: bytemuck::cast_slice(&edges_pod),
        usage: wgpu::BufferUsages::STORAGE,
    });

    let edge_forces_bytes = (n_edges * 6 * std::mem::size_of::<f32>()) as u64;
    let buf_edge_forces = ctx.device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("wlc edge_forces"),
        size: edge_forces_bytes.max(16),
        usage: wgpu::BufferUsages::STORAGE,
        mapped_at_creation: false,
    });

    let buf_csr_offsets = ctx.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("wlc csr_offsets"),
        contents: bytemuck::cast_slice(&backend.csr_offsets),
        usage: wgpu::BufferUsages::STORAGE,
    });
    let buf_csr_data = ctx.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("wlc csr_data"),
        contents: bytemuck::cast_slice(&backend.csr_data),
        usage: wgpu::BufferUsages::STORAGE,
    });

    let node_forces_bytes = (n_nodes * 3 * std::mem::size_of::<f32>()) as u64;
    let buf_node_forces = ctx.device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("wlc node_forces"),
        size: node_forces_bytes,
        usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
        mapped_at_creation: false,
    });

    // kbt = 4.11 pN·nm = 4.11e-9 μN·μm.
    let kbt_uN_um = 4.11_f32 * 1e-9;

    let uniforms = WlcUniforms {
        n_edges: n_edges as u32,
        n_nodes: n_nodes as u32,
        persistence_length_um: params.persistence_length_um,
        kbt_uN_um,
        max_relative_extension: params.max_relative_extension,
        _pad0: 0.0,
        _pad1: 0.0,
        _pad2: 0.0,
    };
    let buf_uniforms = ctx.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("wlc uniforms"),
        contents: bytemuck::bytes_of(&uniforms),
        usage: wgpu::BufferUsages::UNIFORM,
    });

    let buf_staging = ctx.device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("wlc staging"),
        size: node_forces_bytes,
        usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
        mapped_at_creation: false,
    });

    let storage_entry = |binding: u32, read_only: bool| wgpu::BindGroupLayoutEntry {
        binding,
        visibility: wgpu::ShaderStages::COMPUTE,
        ty: wgpu::BindingType::Buffer {
            ty: wgpu::BufferBindingType::Storage { read_only },
            has_dynamic_offset: false,
            min_binding_size: None,
        },
        count: None,
    };

    let bgl = ctx.device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
        label: Some("wlc bgl"),
        entries: &[
            storage_entry(0, true),
            storage_entry(1, true),
            storage_entry(2, false),
            storage_entry(3, true),
            storage_entry(4, true),
            storage_entry(5, false),
            wgpu::BindGroupLayoutEntry {
                binding: 6,
                visibility: wgpu::ShaderStages::COMPUTE,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Uniform,
                    has_dynamic_offset: false,
                    min_binding_size: None,
                },
                count: None,
            },
        ],
    });

    let bg = ctx.device.create_bind_group(&wgpu::BindGroupDescriptor {
        label: Some("wlc bg"),
        layout: &bgl,
        entries: &[
            wgpu::BindGroupEntry { binding: 0, resource: buf_pos.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 1, resource: buf_edges.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 2, resource: buf_edge_forces.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 3, resource: buf_csr_offsets.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 4, resource: buf_csr_data.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 5, resource: buf_node_forces.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 6, resource: buf_uniforms.as_entire_binding() },
        ],
    });

    let shader = ctx.device.create_shader_module(wgpu::ShaderModuleDescriptor {
        label: Some("wlc shader"),
        source: wgpu::ShaderSource::Wgsl(KERNEL_WLC_WGSL.into()),
    });

    let pipeline_layout = ctx.device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
        label: Some("wlc pl"),
        bind_group_layouts: &[&bgl],
        push_constant_ranges: &[],
    });

    let pipe_per_edge = ctx.device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
        label: Some("wlc per_edge pipeline"),
        layout: Some(&pipeline_layout),
        module: &shader,
        entry_point: "wlc_per_edge",
        compilation_options: Default::default(),
    });
    let pipe_aggregate = ctx.device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
        label: Some("wlc aggregate pipeline"),
        layout: Some(&pipeline_layout),
        module: &shader,
        entry_point: "wlc_aggregate",
        compilation_options: Default::default(),
    });

    let mut encoder = ctx.device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
        label: Some("wlc encoder"),
    });

    {
        let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
            label: Some("wlc pass"),
            timestamp_writes: None,
        });
        pass.set_bind_group(0, &bg, &[]);

        pass.set_pipeline(&pipe_per_edge);
        let groups = ((n_edges as u32) + 63) / 64;
        pass.dispatch_workgroups(groups.max(1), 1, 1);

        pass.set_pipeline(&pipe_aggregate);
        let groups = ((n_nodes as u32) + 63) / 64;
        pass.dispatch_workgroups(groups.max(1), 1, 1);
    }

    encoder.copy_buffer_to_buffer(&buf_node_forces, 0, &buf_staging, 0, node_forces_bytes);
    ctx.queue.submit(std::iter::once(encoder.finish()));

    let slice = buf_staging.slice(..);
    let (tx, rx) = mpsc::channel();
    slice.map_async(wgpu::MapMode::Read, move |result| {
        let _ = tx.send(result);
    });
    ctx.device.poll(wgpu::Maintain::Wait);
    rx.recv()
        .context("wlc: map_async channel closed")?
        .map_err(|e| anyhow::anyhow!("wlc: map_async failed: {e:?}"))?;

    let mapped = slice.get_mapped_range();
    let result_f32: &[f32] = bytemuck::cast_slice(&mapped);

    let mut out = Vec::with_capacity(n_nodes);
    for i in 0..n_nodes {
        out.push(Vec3::new(
            result_f32[i * 3 + 0],
            result_f32[i * 3 + 1],
            result_f32[i * 3 + 2],
        ));
    }

    drop(mapped);
    buf_staging.unmap();

    Ok(out)
}

// =============================================================================
// Phase 11.3.C — DPD membrane forces (stateless PCG hash RNG)
// =============================================================================

#[repr(C)]
#[derive(Copy, Clone, Debug, Pod, Zeroable)]
struct DpdUniforms {
    n_vertices: u32,
    step_count: u32,
    seed: u32,
    gamma: f32,
    sigma: f32,
    visc_scale: f32,
    rand_scale: f32,
    _pad0: f32,
}

/// Run the per-vertex DPD membrane force kernel on the GPU. Uses the same
/// stateless PCG hash + Box-Muller as `DPDSolver::compute_membrane_forces_pcg`
/// on the CPU side, so given identical (seed, step_count) the output is
/// bit-identical between the two backends.
pub fn run_dpd_membrane_forces(
    ctx: &ComputeContext,
    velocities: &[Vec3],
    params: &DPDParameters,
    temperature_K: f32,
    step_count: u32,
    seed: u32,
) -> Result<Vec<Vec3>> {
    let n = velocities.len();
    if n == 0 {
        return Ok(Vec::new());
    }

    // Compute σ for the given temperature (mirrors the CPU code).
    let kbt_J = 1.38e-23_f32 * temperature_K;
    let kbt = kbt_J * 1e12_f32;
    let sigma = (2.0_f32 * params.dissipation_gamma * kbt).sqrt();

    // Upload velocities (SoA scalar).
    let mut vel_f32 = vec![0f32; n * 3];
    for (i, v) in velocities.iter().enumerate() {
        vel_f32[i * 3 + 0] = v.x;
        vel_f32[i * 3 + 1] = v.y;
        vel_f32[i * 3 + 2] = v.z;
    }
    let buf_vel = ctx.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("dpd velocities"),
        contents: bytemuck::cast_slice(&vel_f32),
        usage: wgpu::BufferUsages::STORAGE,
    });

    let forces_bytes = (n * 3 * std::mem::size_of::<f32>()) as u64;
    let buf_forces = ctx.device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("dpd forces"),
        size: forces_bytes,
        usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
        mapped_at_creation: false,
    });

    let uniforms = DpdUniforms {
        n_vertices: n as u32,
        step_count,
        seed,
        gamma: params.dissipation_gamma,
        sigma,
        visc_scale: 0.001,
        rand_scale: 0.0001,
        _pad0: 0.0,
    };
    let buf_u = ctx.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("dpd uniforms"),
        contents: bytemuck::bytes_of(&uniforms),
        usage: wgpu::BufferUsages::UNIFORM,
    });

    let buf_staging = ctx.device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("dpd staging"),
        size: forces_bytes,
        usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
        mapped_at_creation: false,
    });

    let bgl = ctx.device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
        label: Some("dpd bgl"),
        entries: &[
            wgpu::BindGroupLayoutEntry {
                binding: 0,
                visibility: wgpu::ShaderStages::COMPUTE,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Storage { read_only: true },
                    has_dynamic_offset: false,
                    min_binding_size: None,
                },
                count: None,
            },
            wgpu::BindGroupLayoutEntry {
                binding: 1,
                visibility: wgpu::ShaderStages::COMPUTE,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Storage { read_only: false },
                    has_dynamic_offset: false,
                    min_binding_size: None,
                },
                count: None,
            },
            wgpu::BindGroupLayoutEntry {
                binding: 2,
                visibility: wgpu::ShaderStages::COMPUTE,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Uniform,
                    has_dynamic_offset: false,
                    min_binding_size: None,
                },
                count: None,
            },
        ],
    });

    let bg = ctx.device.create_bind_group(&wgpu::BindGroupDescriptor {
        label: Some("dpd bg"),
        layout: &bgl,
        entries: &[
            wgpu::BindGroupEntry { binding: 0, resource: buf_vel.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 1, resource: buf_forces.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 2, resource: buf_u.as_entire_binding() },
        ],
    });

    let shader = ctx.device.create_shader_module(wgpu::ShaderModuleDescriptor {
        label: Some("dpd shader"),
        source: wgpu::ShaderSource::Wgsl(KERNEL_DPD_WGSL.into()),
    });

    let pipeline_layout = ctx.device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
        label: Some("dpd pl"),
        bind_group_layouts: &[&bgl],
        push_constant_ranges: &[],
    });

    let pipeline = ctx.device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
        label: Some("dpd pipeline"),
        layout: Some(&pipeline_layout),
        module: &shader,
        entry_point: "dpd_membrane",
        compilation_options: Default::default(),
    });

    let mut encoder = ctx.device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
        label: Some("dpd encoder"),
    });

    {
        let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
            label: Some("dpd pass"),
            timestamp_writes: None,
        });
        pass.set_pipeline(&pipeline);
        pass.set_bind_group(0, &bg, &[]);
        let groups = ((n as u32) + 63) / 64;
        pass.dispatch_workgroups(groups.max(1), 1, 1);
    }

    encoder.copy_buffer_to_buffer(&buf_forces, 0, &buf_staging, 0, forces_bytes);
    ctx.queue.submit(std::iter::once(encoder.finish()));

    let slice = buf_staging.slice(..);
    let (tx, rx) = mpsc::channel();
    slice.map_async(wgpu::MapMode::Read, move |result| {
        let _ = tx.send(result);
    });
    ctx.device.poll(wgpu::Maintain::Wait);
    rx.recv()
        .context("dpd: map_async channel closed")?
        .map_err(|e| anyhow::anyhow!("dpd: map_async failed: {e:?}"))?;

    let mapped = slice.get_mapped_range();
    let result_f32: &[f32] = bytemuck::cast_slice(&mapped);

    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        out.push(Vec3::new(
            result_f32[i * 3 + 0],
            result_f32[i * 3 + 1],
            result_f32[i * 3 + 2],
        ));
    }

    drop(mapped);
    buf_staging.unmap();

    Ok(out)
}

// =============================================================================
// Phase 11.3.D — Velocity-Verlet integrator (half-step + position update)
// =============================================================================

#[repr(C)]
#[derive(Copy, Clone, Debug, Pod, Zeroable)]
struct VerletUniforms {
    n_vertices: u32,
    half_dt: f32,
    dt: f32,
    inv_mass: f32,
    accel_scale: f32,
    max_velocity: f32,
    max_displacement: f32,
    _pad0: f32,
}

/// Configuration for the GPU Velocity-Verlet step.
#[derive(Clone, Debug)]
pub struct VerletConfig {
    pub dt_sec: f32,
    pub vertex_mass_pg: f32,
    pub max_velocity_um_per_sec: f32,
    pub max_displacement_um: f32,
}

impl Default for VerletConfig {
    fn default() -> Self {
        Self {
            dt_sec: 1e-6,
            vertex_mass_pg: 1.0,
            max_velocity_um_per_sec: 1000.0,
            max_displacement_um: 0.1,
        }
    }
}

/// Run one Velocity-Verlet step on the GPU. Mirrors
/// `VelocityVerlet::half_step_velocity` followed by `step_position`. The
/// `forces` slice is read-only; `positions` and `velocities` are updated
/// in place.
///
/// The CPU side has a magic `accel_scale = 10.0` constant — same value is
/// used here so unit conversions stay aligned.
pub fn run_verlet_step(
    ctx: &ComputeContext,
    positions: &mut [Vec3],
    velocities: &mut [Vec3],
    forces: &[Vec3],
    config: &VerletConfig,
) -> Result<()> {
    let n = positions.len();
    if n != velocities.len() || n != forces.len() {
        anyhow::bail!(
            "verlet length mismatch: pos={}, vel={}, forces={}",
            n,
            velocities.len(),
            forces.len()
        );
    }
    if n == 0 {
        return Ok(());
    }

    let mut pos_f32 = vec![0f32; n * 3];
    let mut vel_f32 = vec![0f32; n * 3];
    let mut force_f32 = vec![0f32; n * 3];
    for i in 0..n {
        pos_f32[i * 3 + 0] = positions[i].x;
        pos_f32[i * 3 + 1] = positions[i].y;
        pos_f32[i * 3 + 2] = positions[i].z;
        vel_f32[i * 3 + 0] = velocities[i].x;
        vel_f32[i * 3 + 1] = velocities[i].y;
        vel_f32[i * 3 + 2] = velocities[i].z;
        force_f32[i * 3 + 0] = forces[i].x;
        force_f32[i * 3 + 1] = forces[i].y;
        force_f32[i * 3 + 2] = forces[i].z;
    }

    let buf_forces = ctx.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("verlet forces"),
        contents: bytemuck::cast_slice(&force_f32),
        usage: wgpu::BufferUsages::STORAGE,
    });
    let buf_vel = ctx.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("verlet velocities"),
        contents: bytemuck::cast_slice(&vel_f32),
        usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
    });
    let buf_pos = ctx.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("verlet positions"),
        contents: bytemuck::cast_slice(&pos_f32),
        usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
    });

    let uniforms = VerletUniforms {
        n_vertices: n as u32,
        half_dt: config.dt_sec * 0.5,
        dt: config.dt_sec,
        inv_mass: 1.0 / config.vertex_mass_pg,
        accel_scale: 10.0,
        max_velocity: config.max_velocity_um_per_sec,
        max_displacement: config.max_displacement_um,
        _pad0: 0.0,
    };
    let buf_u = ctx.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("verlet uniforms"),
        contents: bytemuck::bytes_of(&uniforms),
        usage: wgpu::BufferUsages::UNIFORM,
    });

    let bytes = (n * 3 * std::mem::size_of::<f32>()) as u64;
    let buf_pos_staging = ctx.device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("verlet pos staging"),
        size: bytes,
        usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
        mapped_at_creation: false,
    });
    let buf_vel_staging = ctx.device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("verlet vel staging"),
        size: bytes,
        usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
        mapped_at_creation: false,
    });

    let storage_entry = |binding: u32, read_only: bool| wgpu::BindGroupLayoutEntry {
        binding,
        visibility: wgpu::ShaderStages::COMPUTE,
        ty: wgpu::BindingType::Buffer {
            ty: wgpu::BufferBindingType::Storage { read_only },
            has_dynamic_offset: false,
            min_binding_size: None,
        },
        count: None,
    };

    let bgl = ctx.device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
        label: Some("verlet bgl"),
        entries: &[
            storage_entry(0, true),  // forces (read-only)
            storage_entry(1, false), // velocities (read-write)
            storage_entry(2, false), // positions (read-write)
            wgpu::BindGroupLayoutEntry {
                binding: 3,
                visibility: wgpu::ShaderStages::COMPUTE,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Uniform,
                    has_dynamic_offset: false,
                    min_binding_size: None,
                },
                count: None,
            },
        ],
    });

    let bg = ctx.device.create_bind_group(&wgpu::BindGroupDescriptor {
        label: Some("verlet bg"),
        layout: &bgl,
        entries: &[
            wgpu::BindGroupEntry { binding: 0, resource: buf_forces.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 1, resource: buf_vel.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 2, resource: buf_pos.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 3, resource: buf_u.as_entire_binding() },
        ],
    });

    let shader = ctx.device.create_shader_module(wgpu::ShaderModuleDescriptor {
        label: Some("verlet shader"),
        source: wgpu::ShaderSource::Wgsl(KERNEL_INTEGRATOR_WGSL.into()),
    });

    let pl = ctx.device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
        label: Some("verlet pl"),
        bind_group_layouts: &[&bgl],
        push_constant_ranges: &[],
    });
    let pipe_half = ctx.device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
        label: Some("verlet half_step pipeline"),
        layout: Some(&pl),
        module: &shader,
        entry_point: "verlet_half_step_velocity",
        compilation_options: Default::default(),
    });
    let pipe_pos = ctx.device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
        label: Some("verlet step_position pipeline"),
        layout: Some(&pl),
        module: &shader,
        entry_point: "verlet_step_position",
        compilation_options: Default::default(),
    });

    let mut encoder = ctx.device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
        label: Some("verlet encoder"),
    });
    {
        let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
            label: Some("verlet pass"),
            timestamp_writes: None,
        });
        pass.set_bind_group(0, &bg, &[]);
        let groups = ((n as u32) + 63) / 64;
        pass.set_pipeline(&pipe_half);
        pass.dispatch_workgroups(groups.max(1), 1, 1);
        pass.set_pipeline(&pipe_pos);
        pass.dispatch_workgroups(groups.max(1), 1, 1);
    }

    encoder.copy_buffer_to_buffer(&buf_pos, 0, &buf_pos_staging, 0, bytes);
    encoder.copy_buffer_to_buffer(&buf_vel, 0, &buf_vel_staging, 0, bytes);
    ctx.queue.submit(std::iter::once(encoder.finish()));

    let pos_slice = buf_pos_staging.slice(..);
    let (tx_pos, rx_pos) = mpsc::channel();
    pos_slice.map_async(wgpu::MapMode::Read, move |r| { let _ = tx_pos.send(r); });
    let vel_slice = buf_vel_staging.slice(..);
    let (tx_vel, rx_vel) = mpsc::channel();
    vel_slice.map_async(wgpu::MapMode::Read, move |r| { let _ = tx_vel.send(r); });

    ctx.device.poll(wgpu::Maintain::Wait);
    rx_pos.recv()
        .context("verlet: pos map_async channel closed")?
        .map_err(|e| anyhow::anyhow!("verlet: pos map_async failed: {e:?}"))?;
    rx_vel.recv()
        .context("verlet: vel map_async channel closed")?
        .map_err(|e| anyhow::anyhow!("verlet: vel map_async failed: {e:?}"))?;

    let mapped_pos = pos_slice.get_mapped_range();
    let pos_out: &[f32] = bytemuck::cast_slice(&mapped_pos);
    let mapped_vel = vel_slice.get_mapped_range();
    let vel_out: &[f32] = bytemuck::cast_slice(&mapped_vel);

    for i in 0..n {
        positions[i] = Vec3::new(
            pos_out[i * 3 + 0],
            pos_out[i * 3 + 1],
            pos_out[i * 3 + 2],
        );
        velocities[i] = Vec3::new(
            vel_out[i * 3 + 0],
            vel_out[i * 3 + 1],
            vel_out[i * 3 + 2],
        );
    }

    drop(mapped_pos);
    drop(mapped_vel);
    buf_pos_staging.unmap();
    buf_vel_staging.unmap();

    Ok(())
}

// =============================================================================
// Phase 11.3.E — Integrated single-cell physics backend (persistent buffers)
// =============================================================================

#[repr(C)]
#[derive(Copy, Clone, Debug, Pod, Zeroable)]
struct StepUniforms {
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

/// Configuration for the integrated physics backend.
#[derive(Clone, Debug)]
pub struct PhysicsBackendConfig {
    pub dt_sec: f32,
    pub vertex_mass_pg: f32,
    pub max_velocity_um_per_sec: f32,
    pub max_displacement_um: f32,
    pub temperature_K: f32,
    pub membrane_damping: f32,
    pub seed: u32,
}

impl Default for PhysicsBackendConfig {
    fn default() -> Self {
        Self {
            dt_sec: 1e-6,
            vertex_mass_pg: 1.0,
            max_velocity_um_per_sec: 1000.0,
            max_displacement_um: 0.1,
            temperature_K: 310.0,
            membrane_damping: 0.1,
            seed: 0xDEADBEEF,
        }
    }
}

/// Persistent GPU backend for the integrated physics step.
///
/// Owns all persistent buffers (positions, velocities, forces, mesh
/// topology, WLC baseline) and the compiled compute pipelines. Each call
/// to `step()` runs one substep on the GPU; `read_state()` downloads
/// positions and velocities back to the host.
///
/// WLC forces are baked into a `wlc_baseline` buffer at construction
/// because the CPU `WLCSolver::compute_forces` reads from
/// `network.nodes` (which never updates as the mesh moves) — the contribution
/// is therefore static across substeps. See `docs/phase_11_3_notes.md`.
pub struct PhysicsBackend {
    n_vertices: u32,
    n_elements: u32,
    skalak_material: SkalakMaterial,
    config: PhysicsBackendConfig,
    step_count: u32,
    #[allow(dead_code)]
    pipeline_layout: wgpu::PipelineLayout,
    #[allow(dead_code)]
    bgl: wgpu::BindGroupLayout,
    bg: wgpu::BindGroup,
    pipe_half_step: wgpu::ComputePipeline,
    pipe_step_position: wgpu::ComputePipeline,
    pipe_skalak_per_element: wgpu::ComputePipeline,
    pipe_skalak_init: wgpu::ComputePipeline,
    pipe_dpd_add: wgpu::ComputePipeline,
    pipe_damping_half: wgpu::ComputePipeline,
    buf_positions: wgpu::Buffer,
    buf_velocities: wgpu::Buffer,
    buf_forces: wgpu::Buffer,
    buf_uniforms: wgpu::Buffer,
    #[allow(dead_code)]
    buf_wlc_baseline: wgpu::Buffer,
    #[allow(dead_code)]
    buf_elements: wgpu::Buffer,
    #[allow(dead_code)]
    buf_csr_offsets: wgpu::Buffer,
    #[allow(dead_code)]
    buf_csr_data: wgpu::Buffer,
    #[allow(dead_code)]
    buf_elem_forces: wgpu::Buffer,
    /// Phase 12.B.1: per-vertex external forces (e.g., flow drag).
    /// Default-zero; host updates via `set_external_forces`.
    buf_external_forces: wgpu::Buffer,
    pos_bytes: u64,
    device: std::sync::Arc<wgpu::Device>,
    queue: std::sync::Arc<wgpu::Queue>,
}

impl PhysicsBackend {
    /// Construct a new persistent backend for a given mesh + spectrin
    /// network. WLC forces are computed once on the GPU using the existing
    /// `run_wlc_forces` path and uploaded as a per-mesh-vertex baseline.
    pub fn new(
        ctx: &ComputeContext,
        mesh: &Mesh,
        spectrin: &SpectrinNetwork,
        skalak_solver: &SkalakSolver,
        wlc_params: &WLCParameters,
        skalak_material: SkalakMaterial,
        config: PhysicsBackendConfig,
    ) -> Result<Self> {
        let n_vertices = mesh.vertices.len() as u32;

        // === WLC baseline: per-mesh-vertex forces from spectrin (static) ===
        let wlc_backend = WlcBackendData::from_network(spectrin);
        let wlc_per_node = run_wlc_forces(ctx, spectrin, &wlc_backend, wlc_params)?;
        let mut wlc_per_vertex = vec![Vec3::ZERO; n_vertices as usize];
        for (node_idx, node) in spectrin.nodes.iter().enumerate() {
            let v = node.mesh_vertex_idx as usize;
            if v < wlc_per_vertex.len() {
                wlc_per_vertex[v] += wlc_per_node[node_idx];
            }
        }
        let mut wlc_f32 = vec![0f32; n_vertices as usize * 3];
        for (i, f) in wlc_per_vertex.iter().enumerate() {
            wlc_f32[i * 3 + 0] = f.x;
            wlc_f32[i * 3 + 1] = f.y;
            wlc_f32[i * 3 + 2] = f.z;
        }

        // === Skalak topology ===
        let skalak_backend = SkalakBackendData::from_solver(skalak_solver, n_vertices as usize);

        // === Initial positions (from mesh) ===
        let mut pos_f32 = vec![0f32; n_vertices as usize * 3];
        for (i, v) in mesh.vertices.iter().enumerate() {
            pos_f32[i * 3 + 0] = v.position[0];
            pos_f32[i * 3 + 1] = v.position[1];
            pos_f32[i * 3 + 2] = v.position[2];
        }

        // === Persistent buffer allocation ===
        let pos_bytes = (n_vertices as u64) * 3 * std::mem::size_of::<f32>() as u64;
        let buf_positions = ctx.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("backend positions"),
            contents: bytemuck::cast_slice(&pos_f32),
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
        });
        let zero_pos: Vec<f32> = vec![0.0; n_vertices as usize * 3];
        let buf_velocities = ctx.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("backend velocities"),
            contents: bytemuck::cast_slice(&zero_pos),
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
        });
        let buf_forces = ctx.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("backend forces"),
            contents: bytemuck::cast_slice(&zero_pos),
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
        });
        let buf_wlc_baseline = ctx.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("backend wlc_baseline"),
            contents: bytemuck::cast_slice(&wlc_f32),
            usage: wgpu::BufferUsages::STORAGE,
        });
        // Phase 12.B.1: zero-initialized external forces buffer; host
        // writes Poiseuille drag (or any other per-vertex force field)
        // via `PhysicsBackend::set_external_forces`.
        let buf_external_forces = ctx.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("backend external_forces"),
            contents: bytemuck::cast_slice(&zero_pos),
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
        });

        let elements_pod: Vec<ElementGpu> = skalak_backend
            .elements
            .iter()
            .map(|e| ElementGpu {
                v0: e.v0,
                v1: e.v1,
                v2: e.v2,
                ref_area: e.ref_area,
                ref_len_01: e.ref_len_01,
                ref_len_02: e.ref_len_02,
                ref_len_12: e.ref_len_12,
                _pad: 0.0,
            })
            .collect();
        let buf_elements = ctx.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("backend elements"),
            contents: bytemuck::cast_slice(&elements_pod),
            usage: wgpu::BufferUsages::STORAGE,
        });
        let buf_csr_offsets = ctx.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("backend csr_offsets"),
            contents: bytemuck::cast_slice(&skalak_backend.csr_offsets),
            usage: wgpu::BufferUsages::STORAGE,
        });
        let buf_csr_data = ctx.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("backend csr_data"),
            contents: bytemuck::cast_slice(&skalak_backend.csr_data),
            usage: wgpu::BufferUsages::STORAGE,
        });

        let n_elements = skalak_backend.n_elements as u32;
        let elem_forces_bytes = (n_elements as u64) * 9 * std::mem::size_of::<f32>() as u64;
        let buf_elem_forces = ctx.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("backend elem_forces"),
            size: elem_forces_bytes.max(16),
            usage: wgpu::BufferUsages::STORAGE,
            mapped_at_creation: false,
        });

        let buf_uniforms = ctx.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("backend uniforms"),
            size: std::mem::size_of::<StepUniforms>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        // === Bind group layout ===
        let storage = |b: u32, ro: bool| wgpu::BindGroupLayoutEntry {
            binding: b,
            visibility: wgpu::ShaderStages::COMPUTE,
            ty: wgpu::BindingType::Buffer {
                ty: wgpu::BufferBindingType::Storage { read_only: ro },
                has_dynamic_offset: false,
                min_binding_size: None,
            },
            count: None,
        };
        let bgl = ctx.device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("backend bgl"),
            entries: &[
                storage(0, false), // positions
                storage(1, false), // velocities
                storage(2, false), // forces
                storage(3, true),  // wlc_baseline
                storage(4, true),  // elements
                storage(5, true),  // csr_offsets
                storage(6, true),  // csr_data
                storage(7, false), // elem_forces
                wgpu::BindGroupLayoutEntry {
                    binding: 8,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Uniform,
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
                storage(9, true),  // external_forces (Phase 12.B.1)
            ],
        });

        let bg = ctx.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("backend bg"),
            layout: &bgl,
            entries: &[
                wgpu::BindGroupEntry { binding: 0, resource: buf_positions.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 1, resource: buf_velocities.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 2, resource: buf_forces.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 3, resource: buf_wlc_baseline.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 4, resource: buf_elements.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 5, resource: buf_csr_offsets.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 6, resource: buf_csr_data.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 7, resource: buf_elem_forces.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 8, resource: buf_uniforms.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 9, resource: buf_external_forces.as_entire_binding() },
            ],
        });

        let shader = ctx.device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("backend shader"),
            source: wgpu::ShaderSource::Wgsl(KERNEL_STEP_WGSL.into()),
        });
        let pipeline_layout = ctx.device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("backend pl"),
            bind_group_layouts: &[&bgl],
            push_constant_ranges: &[],
        });
        let make_pipe = |entry: &str| {
            ctx.device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                label: Some(entry),
                layout: Some(&pipeline_layout),
                module: &shader,
                entry_point: entry,
                compilation_options: Default::default(),
            })
        };
        let pipe_half_step = make_pipe("verlet_half_step_velocity");
        let pipe_step_position = make_pipe("verlet_step_position");
        let pipe_skalak_per_element = make_pipe("skalak_per_element");
        let pipe_skalak_init = make_pipe("skalak_init_from_baseline");
        let pipe_dpd_add = make_pipe("dpd_add");
        let pipe_damping_half = make_pipe("damping_and_half_step");

        Ok(Self {
            n_vertices,
            n_elements,
            skalak_material,
            config,
            step_count: 0,
            pipeline_layout,
            bgl,
            bg,
            pipe_half_step,
            pipe_step_position,
            pipe_skalak_per_element,
            pipe_skalak_init,
            pipe_dpd_add,
            pipe_damping_half,
            buf_positions,
            buf_velocities,
            buf_forces,
            buf_uniforms,
            buf_wlc_baseline,
            buf_elements,
            buf_csr_offsets,
            buf_csr_data,
            buf_elem_forces,
            buf_external_forces,
            pos_bytes,
            device: ctx.device.clone(),
            queue: ctx.queue.clone(),
        })
    }

    /// Phase 12.B.1: upload per-vertex external forces (Poiseuille drag,
    /// any other force field). Picked up by the next call to `step()`
    /// inside `skalak_init_from_baseline`. The buffer persists across
    /// substeps; call this once per substep with the latest drag (or
    /// once per outer flow update if drag is quasi-stationary).
    pub fn set_external_forces(&self, forces: &[Vec3]) -> Result<()> {
        anyhow::ensure!(
            forces.len() == self.n_vertices as usize,
            "external_forces length {} does not match n_vertices {}",
            forces.len(), self.n_vertices
        );
        let mut buf = vec![0f32; forces.len() * 3];
        for (i, f) in forces.iter().enumerate() {
            buf[i * 3 + 0] = f.x;
            buf[i * 3 + 1] = f.y;
            buf[i * 3 + 2] = f.z;
        }
        self.queue
            .write_buffer(&self.buf_external_forces, 0, bytemuck::cast_slice(&buf));
        Ok(())
    }

    /// Reset external forces to zero. Equivalent to
    /// `set_external_forces(&vec![Vec3::ZERO; n_vertices])` but cheaper.
    pub fn clear_external_forces(&self) {
        let zero = vec![0f32; self.n_vertices as usize * 3];
        self.queue
            .write_buffer(&self.buf_external_forces, 0, bytemuck::cast_slice(&zero));
    }

    /// Run one physics substep on the GPU. Mirrors `PhysicsSolver::step`
    /// (with PCG-DPD rather than `StdRng`-DPD).
    pub fn step(&mut self) -> Result<()> {
        // Compute σ for current temperature.
        let kbt_J = 1.38e-23_f32 * self.config.temperature_K;
        let kbt = kbt_J * 1e12_f32;
        let gamma = 4.5_f32; // matches DPDParameters::default
        let sigma = (2.0_f32 * gamma * kbt).sqrt();

        let uniforms = StepUniforms {
            n_vertices: self.n_vertices,
            n_elements: self.n_elements,
            half_dt: self.config.dt_sec * 0.5,
            dt: self.config.dt_sec,
            inv_mass: 1.0 / self.config.vertex_mass_pg,
            accel_scale: 10.0,
            max_velocity: self.config.max_velocity_um_per_sec,
            max_displacement: self.config.max_displacement_um,
            k_spring: self.skalak_material.shear_modulus_uN_per_m * 1.732,
            area_modulus: self.skalak_material.area_modulus_uN_per_m,
            area_force_scale: 0.001,
            membrane_damping: self.config.membrane_damping,
            gamma,
            sigma,
            visc_scale: 0.001,
            rand_scale: 0.0001,
            step_count: self.step_count,
            seed: self.config.seed,
            _pad0: 0.0,
            _pad1: 0.0,
        };
        self.queue
            .write_buffer(&self.buf_uniforms, 0, bytemuck::bytes_of(&uniforms));

        let mut encoder = self.device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("backend step encoder"),
        });
        {
            let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("backend step pass"),
                timestamp_writes: None,
            });
            pass.set_bind_group(0, &self.bg, &[]);
            let groups_v = ((self.n_vertices) + 63) / 64;
            let groups_e = ((self.n_elements) + 63) / 64;

            // 1. Half-step velocity using current forces.
            pass.set_pipeline(&self.pipe_half_step);
            pass.dispatch_workgroups(groups_v.max(1), 1, 1);

            // 2. Update positions.
            pass.set_pipeline(&self.pipe_step_position);
            pass.dispatch_workgroups(groups_v.max(1), 1, 1);

            // 3. Skalak per-element compute.
            pass.set_pipeline(&self.pipe_skalak_per_element);
            pass.dispatch_workgroups(groups_e.max(1), 1, 1);

            // 4. Skalak aggregate + WLC baseline init.
            pass.set_pipeline(&self.pipe_skalak_init);
            pass.dispatch_workgroups(groups_v.max(1), 1, 1);

            // 5. DPD add.
            pass.set_pipeline(&self.pipe_dpd_add);
            pass.dispatch_workgroups(groups_v.max(1), 1, 1);

            // 6. Membrane damping + final half-step.
            pass.set_pipeline(&self.pipe_damping_half);
            pass.dispatch_workgroups(groups_v.max(1), 1, 1);
        }
        self.queue.submit(std::iter::once(encoder.finish()));
        self.step_count += 1;
        Ok(())
    }

    /// Download current positions and velocities to the host.
    pub fn read_state(&self) -> Result<(Vec<Vec3>, Vec<Vec3>)> {
        let bytes = self.pos_bytes;
        let pos_staging = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("backend pos staging"),
            size: bytes,
            usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
            mapped_at_creation: false,
        });
        let vel_staging = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("backend vel staging"),
            size: bytes,
            usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
            mapped_at_creation: false,
        });
        let mut encoder = self.device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("backend read encoder"),
        });
        encoder.copy_buffer_to_buffer(&self.buf_positions, 0, &pos_staging, 0, bytes);
        encoder.copy_buffer_to_buffer(&self.buf_velocities, 0, &vel_staging, 0, bytes);
        self.queue.submit(std::iter::once(encoder.finish()));

        let pos_slice = pos_staging.slice(..);
        let (tx_pos, rx_pos) = mpsc::channel();
        pos_slice.map_async(wgpu::MapMode::Read, move |r| { let _ = tx_pos.send(r); });
        let vel_slice = vel_staging.slice(..);
        let (tx_vel, rx_vel) = mpsc::channel();
        vel_slice.map_async(wgpu::MapMode::Read, move |r| { let _ = tx_vel.send(r); });

        self.device.poll(wgpu::Maintain::Wait);
        rx_pos.recv()
            .context("backend: pos map_async channel closed")?
            .map_err(|e| anyhow::anyhow!("backend: pos map_async failed: {e:?}"))?;
        rx_vel.recv()
            .context("backend: vel map_async channel closed")?
            .map_err(|e| anyhow::anyhow!("backend: vel map_async failed: {e:?}"))?;

        let mapped_pos = pos_slice.get_mapped_range();
        let pos_out: &[f32] = bytemuck::cast_slice(&mapped_pos);
        let mapped_vel = vel_slice.get_mapped_range();
        let vel_out: &[f32] = bytemuck::cast_slice(&mapped_vel);

        let n = self.n_vertices as usize;
        let mut positions = Vec::with_capacity(n);
        let mut velocities = Vec::with_capacity(n);
        for i in 0..n {
            positions.push(Vec3::new(
                pos_out[i * 3 + 0],
                pos_out[i * 3 + 1],
                pos_out[i * 3 + 2],
            ));
            velocities.push(Vec3::new(
                vel_out[i * 3 + 0],
                vel_out[i * 3 + 1],
                vel_out[i * 3 + 2],
            ));
        }

        drop(mapped_pos);
        drop(mapped_vel);
        pos_staging.unmap();
        vel_staging.unmap();

        Ok((positions, velocities))
    }

    /// Current step count, useful for diagnostics.
    pub fn step_count(&self) -> u32 {
        self.step_count
    }
}
