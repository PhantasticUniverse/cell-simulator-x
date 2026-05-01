//! GPU biochemistry kernel host wrapper (Phase 11.2.A).
//!
//! Uploads a flat array of N cells × 17 glycolytic species (f32, cell-major)
//! to a storage buffer, dispatches `rk4_step` for `n_steps` RK4 ms-steps
//! per cell, and reads results back. The kernel matches the CPU
//! [`crate::biochemistry::MetabolismSolver`] structure: 11 glycolysis
//! enzymes + glucose transport + lactate export + external ATP demand.
//!
//! ## Scope (11.2.A)
//!
//! Glycolysis backbone only. The 2,3-BPG shunt, PPP, glutathione cycle,
//! Piezo1, ion homeostasis, hemoglobin, and pH buffer are deferred to
//! 11.2.B/C/D in subsequent sessions. The full 38-species port is a
//! mechanical extension of this kernel — every additional pathway adds
//! ~5–8 species and 2–7 enzyme rate functions matching this template.
//!
//! ## Numerical strategy
//!
//! CPU: f64 throughout (existing `MetabolismSolver`).
//! GPU: f32 throughout. Conversion happens at upload/download.
//! Glycolytic species span 0.001–5 mM, well within f32 precision.
//! Match between CPU and GPU is checked at < 1% relative error per species
//! after a 1-second simulation in [`tests/biochem_gpu_parity.rs`].

use std::sync::mpsc;

use anyhow::{Context as _, Result};
use bytemuck::{Pod, Zeroable};
use wgpu::util::DeviceExt;

use crate::biochemistry::{ExtendedMetaboliteIndices, MetabolitePool};

use super::ComputeContext;

const N_SPECIES: usize = 18;

/// Source of the glycolysis kernel.
const KERNEL_WGSL: &str = include_str!("../../shaders/compute/biochem_glycolysis.wgsl");

/// World-level uniform block (must match WGSL `struct U`).
#[repr(C)]
#[derive(Copy, Clone, Debug, Pod, Zeroable)]
struct GlycolysisUniforms {
    n_cells: u32,
    n_steps: u32,
    dt_sec: f32,
    external_glucose_mM: f32,
    atp_consumption_mM_per_sec: f32,
    enable_glucose_transport: u32,
    enable_lactate_export: u32,
    max_change_mM: f32,
    min_concentration_mM: f32,
    _pad0: f32,
    _pad1: f32,
    _pad2: f32,
}

/// Configuration for one glycolysis batch dispatch.
#[derive(Clone, Debug)]
pub struct GlycolysisBatchConfig {
    /// Integration timestep in seconds (matches CPU integrator default 1e-3).
    pub dt_sec: f64,
    /// Number of RK4 substeps per dispatch.
    pub n_steps: u32,
    /// External plasma glucose concentration in mM.
    pub external_glucose_mM: f64,
    /// External ATP consumption rate in mM/s (membrane pumps, etc.).
    pub atp_consumption_mM_per_sec: f64,
    /// Whether GLUT1-style glucose transport is enabled.
    pub enable_glucose_transport: bool,
    /// Whether MCT1-style lactate export is enabled.
    pub enable_lactate_export: bool,
    /// Per-step max-change clamp (mM); matches CPU `IntegratorConfig`.
    pub max_change_mM: f64,
    /// Minimum concentration floor (mM); matches CPU `IntegratorConfig`.
    pub min_concentration_mM: f64,
}

impl Default for GlycolysisBatchConfig {
    fn default() -> Self {
        Self {
            dt_sec: 0.001,
            n_steps: 1000, // 1 s of simulated time per dispatch.
            external_glucose_mM: 5.0,
            atp_consumption_mM_per_sec: 0.001,
            enable_glucose_transport: true,
            enable_lactate_export: true,
            max_change_mM: 0.1,
            min_concentration_mM: 1e-9,
        }
    }
}

/// Run a multi-cell glycolysis batch on the GPU.
///
/// `pools` is N parallel cells, each with their own metabolite pool. The
/// pools are mutated in place. Only the 17 glycolysis species are
/// touched; other indices (PPP, redox, ions) pass through unchanged.
///
/// One GPU dispatch advances every cell by `config.n_steps * config.dt_sec`
/// of simulated time. Compute pipeline + bind group layout are recreated
/// per call — for production hot-loop use, cache them in a longer-lived
/// `GlycolysisBackend` struct (added in 11.2.B once the kernel is
/// frozen).
pub fn run_glycolysis_batch(
    ctx: &ComputeContext,
    pools: &mut [MetabolitePool],
    config: &GlycolysisBatchConfig,
) -> Result<()> {
    let n_cells = pools.len();
    if n_cells == 0 {
        return Ok(());
    }

    let indices = ExtendedMetaboliteIndices::default();
    let g = &indices.glycolysis;

    // === Upload: flatten N pools × 18 f64 into one Vec<f32> ============
    let mut state_f32 = vec![0f32; n_cells * N_SPECIES];
    for (cell_idx, pool) in pools.iter().enumerate() {
        let base = cell_idx * N_SPECIES;
        let s = &mut state_f32[base..base + N_SPECIES];
        s[0]  = pool.get(g.glucose) as f32;
        s[1]  = pool.get(g.glucose_6_phosphate) as f32;
        s[2]  = pool.get(g.fructose_6_phosphate) as f32;
        s[3]  = pool.get(g.fructose_1_6_bisphosphate) as f32;
        s[4]  = pool.get(g.dihydroxyacetone_phosphate) as f32;
        s[5]  = pool.get(g.glyceraldehyde_3_phosphate) as f32;
        s[6]  = pool.get(g.bisphosphoglycerate_1_3) as f32;
        s[7]  = pool.get(g.phosphoglycerate_3) as f32;
        s[8]  = pool.get(g.phosphoglycerate_2) as f32;
        s[9]  = pool.get(g.phosphoenolpyruvate) as f32;
        s[10] = pool.get(g.pyruvate) as f32;
        s[11] = pool.get(g.lactate) as f32;
        s[12] = pool.get(g.atp) as f32;
        s[13] = pool.get(g.adp) as f32;
        s[14] = pool.get(g.nad) as f32;
        s[15] = pool.get(g.nadh) as f32;
        s[16] = pool.get(g.pi) as f32;
        s[17] = pool.get(indices.bisphosphoglycerate_2_3) as f32;
    }

    let bytes = (state_f32.len() * std::mem::size_of::<f32>()) as u64;

    // Storage buffer holds state; we write to it from the host then read
    // back via a staging buffer.
    let buf_state = ctx.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("glycolysis state"),
        contents: bytemuck::cast_slice(&state_f32),
        usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
    });

    // Uniforms.
    let uniforms = GlycolysisUniforms {
        n_cells: n_cells as u32,
        n_steps: config.n_steps,
        dt_sec: config.dt_sec as f32,
        external_glucose_mM: config.external_glucose_mM as f32,
        atp_consumption_mM_per_sec: config.atp_consumption_mM_per_sec as f32,
        enable_glucose_transport: u32::from(config.enable_glucose_transport),
        enable_lactate_export: u32::from(config.enable_lactate_export),
        max_change_mM: config.max_change_mM as f32,
        min_concentration_mM: config.min_concentration_mM as f32,
        _pad0: 0.0,
        _pad1: 0.0,
        _pad2: 0.0,
    };
    let buf_u = ctx.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("glycolysis uniforms"),
        contents: bytemuck::bytes_of(&uniforms),
        usage: wgpu::BufferUsages::UNIFORM,
    });

    let buf_staging = ctx.device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("glycolysis staging"),
        size: bytes,
        usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
        mapped_at_creation: false,
    });

    let bgl = ctx.device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
        label: Some("glycolysis bgl"),
        entries: &[
            wgpu::BindGroupLayoutEntry {
                binding: 0,
                visibility: wgpu::ShaderStages::COMPUTE,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Storage { read_only: false },
                    has_dynamic_offset: false,
                    min_binding_size: None,
                },
                count: None,
            },
            wgpu::BindGroupLayoutEntry {
                binding: 1,
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
        label: Some("glycolysis bg"),
        layout: &bgl,
        entries: &[
            wgpu::BindGroupEntry { binding: 0, resource: buf_state.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 1, resource: buf_u.as_entire_binding() },
        ],
    });

    let shader = ctx.device.create_shader_module(wgpu::ShaderModuleDescriptor {
        label: Some("glycolysis shader"),
        source: wgpu::ShaderSource::Wgsl(KERNEL_WGSL.into()),
    });

    let pipeline_layout = ctx.device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
        label: Some("glycolysis pl"),
        bind_group_layouts: &[&bgl],
        push_constant_ranges: &[],
    });

    let pipeline = ctx.device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
        label: Some("glycolysis pipeline"),
        layout: Some(&pipeline_layout),
        module: &shader,
        entry_point: "rk4_step",
        compilation_options: Default::default(),
    });

    let mut encoder = ctx.device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
        label: Some("glycolysis encoder"),
    });

    {
        let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
            label: Some("glycolysis pass"),
            timestamp_writes: None,
        });
        pass.set_pipeline(&pipeline);
        pass.set_bind_group(0, &bg, &[]);
        let groups = ((n_cells as u32) + 63) / 64;
        pass.dispatch_workgroups(groups, 1, 1);
    }

    encoder.copy_buffer_to_buffer(&buf_state, 0, &buf_staging, 0, bytes);
    ctx.queue.submit(std::iter::once(encoder.finish()));

    // === Download ======================================================
    let slice = buf_staging.slice(..);
    let (tx, rx) = mpsc::channel();
    slice.map_async(wgpu::MapMode::Read, move |result| {
        let _ = tx.send(result);
    });
    ctx.device.poll(wgpu::Maintain::Wait);
    rx.recv()
        .context("glycolysis: map_async channel closed")?
        .map_err(|e| anyhow::anyhow!("glycolysis: map_async failed: {e:?}"))?;

    let mapped = slice.get_mapped_range();
    let result_f32: &[f32] = bytemuck::cast_slice(&mapped);

    // === Distribute back to per-cell pools =============================
    for (cell_idx, pool) in pools.iter_mut().enumerate() {
        let base = cell_idx * N_SPECIES;
        let s = &result_f32[base..base + N_SPECIES];
        pool.set(g.glucose, s[0] as f64);
        pool.set(g.glucose_6_phosphate, s[1] as f64);
        pool.set(g.fructose_6_phosphate, s[2] as f64);
        pool.set(g.fructose_1_6_bisphosphate, s[3] as f64);
        pool.set(g.dihydroxyacetone_phosphate, s[4] as f64);
        pool.set(g.glyceraldehyde_3_phosphate, s[5] as f64);
        pool.set(g.bisphosphoglycerate_1_3, s[6] as f64);
        pool.set(g.phosphoglycerate_3, s[7] as f64);
        pool.set(g.phosphoglycerate_2, s[8] as f64);
        pool.set(g.phosphoenolpyruvate, s[9] as f64);
        pool.set(g.pyruvate, s[10] as f64);
        pool.set(g.lactate, s[11] as f64);
        pool.set(g.atp, s[12] as f64);
        pool.set(g.adp, s[13] as f64);
        pool.set(g.nad, s[14] as f64);
        pool.set(g.nadh, s[15] as f64);
        pool.set(g.pi, s[16] as f64);
        pool.set(indices.bisphosphoglycerate_2_3, s[17] as f64);
    }

    drop(mapped);
    buf_staging.unmap();

    Ok(())
}
