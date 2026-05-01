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

use crate::geometry::Mesh;
use crate::physics::{SkalakMaterial, SkalakSolver};

use super::ComputeContext;

const KERNEL_WGSL: &str = include_str!("../../shaders/compute/skalak.wgsl");

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
