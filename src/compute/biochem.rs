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

use crate::biochemistry::{
    AdairConstants, ExtendedMetaboliteIndices, FullyIntegratedIndices, MetabolitePool,
};

use super::ComputeContext;

const N_SPECIES: usize = 18;

/// Number of species in the full 38-species kernel (Phase 11.2.C).
const N_SPECIES_FULL: usize = 38;

/// Source of the glycolysis kernel (18 species — Phases 11.2.A/B).
const KERNEL_WGSL: &str = include_str!("../../shaders/compute/biochem_glycolysis.wgsl");

/// Source of the full biochem kernel (38 species — Phase 11.2.C).
const KERNEL_FULL_WGSL: &str = include_str!("../../shaders/compute/biochem_full.wgsl");

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

// =============================================================================
// Phase 11.2.C/D — Full 38-species biochem kernel
// =============================================================================

/// Uniforms for the full biochem kernel (must match WGSL `struct U`).
#[repr(C)]
#[derive(Copy, Clone, Debug, Pod, Zeroable)]
struct FullBiochemUniforms {
    n_cells: u32,
    n_steps: u32,
    dt_sec: f32,
    external_glucose_mM: f32,
    atp_consumption_mM_per_sec: f32,
    enable_glucose_transport: u32,
    enable_lactate_export: u32,
    max_change_mM: f32,
    min_concentration_mM: f32,
    membrane_tension_pN_per_nm: f32,
    enable_piezo1: u32,
    enable_ion_homeostasis: u32,
    h2o2_production_rate_mM_per_sec: f32,
    na_external_mM: f32,
    k_external_mM: f32,
    g_na_per_sec: f32,
    g_k_per_sec: f32,
    // Phase 11.2.D additions:
    po2_mmHg: f32,
    pco2_mmHg: f32,
    temperature_K: f32,
    total_hb_mM: f32,
    base_p50_mmHg: f32,
    bohr_coefficient: f32,
    dpg_effect_mmHg_per_mM: f32,
    dpg_reference_mM: f32,
    delta_h_kcal_per_mol: f32,
    temperature_reference_K: f32,
    co2_effect_per_mmHg: f32,
    pco2_reference_mmHg: f32,
    k_on_per_mM_per_sec: f32,
    k_off_per_sec: f32,
    buffer_capacity_slykes: f32,
    ph_reference: f32,
    baseline_lactate_mM: f32,
    min_ph: f32,
    max_ph: f32,
    enable_hb: u32,
    adair_k1: f32,
    adair_k2: f32,
    adair_k3: f32,
    adair_k4: f32,
    _pad0: f32,
    _pad1: f32,
    _pad2: f32,
}

/// Configuration for one full biochem batch dispatch.
///
/// Mirrors `FullyIntegratedConfig` (CPU side) for the parameters consumed by
/// the GPU kernel. The three inline homeostasis hacks (basal NADPH
/// consumption, basal GSH oxidation, ATP regen) remain CPU-side, deferred to
/// 11.2.E.
#[derive(Clone, Debug)]
pub struct FullBiochemBatchConfig {
    /// Integration timestep in seconds.
    pub dt_sec: f64,
    /// Number of RK4 substeps per dispatch.
    pub n_steps: u32,
    /// External plasma glucose concentration in mM.
    pub external_glucose_mM: f64,
    /// External ATP consumption rate in mM/s (non-pump demand).
    pub atp_consumption_mM_per_sec: f64,
    /// Enable GLUT1-style glucose transport.
    pub enable_glucose_transport: bool,
    /// Enable MCT1-style lactate export.
    pub enable_lactate_export: bool,
    /// Per-step max-change clamp in mM (mirrors CPU `IntegratorConfig`).
    pub max_change_mM: f64,
    /// Minimum concentration floor in mM.
    pub min_concentration_mM: f64,
    /// Membrane tension fed into the Piezo1 Hill kinetics (pN/nm).
    pub membrane_tension_pN_per_nm: f64,
    /// Enable Piezo1 + PMCA + ATP-release block.
    pub enable_piezo1: bool,
    /// Enable Na⁺/K⁺-ATPase + Na⁺/K⁺ leaks.
    pub enable_ion_homeostasis: bool,
    /// Basal H₂O₂ production rate from Hb autoxidation (mM/s).
    pub h2o2_production_rate_mM_per_sec: f64,
    /// External Na⁺ concentration (plasma, mM).
    pub na_external_mM: f64,
    /// External K⁺ concentration (plasma, mM).
    pub k_external_mM: f64,
    /// Na⁺ leak conductance (mM/s per mM gradient).
    pub g_na_per_sec: f64,
    /// K⁺ leak conductance.
    pub g_k_per_sec: f64,

    // === Phase 11.2.D: Hb + pH ===

    /// Enable post-RK4 hemoglobin Adair update + pH buffer. When false
    /// (default), the kernel runs identically to 11.2.C and the `hb_state`
    /// buffer is untouched (caller may pass a dummy slice).
    pub enable_hb: bool,
    /// Constant world-level oxygen partial pressure (mmHg).
    pub po2_mmHg: f64,
    /// Constant world-level CO₂ partial pressure (mmHg).
    pub pco2_mmHg: f64,
    /// Temperature in Kelvin (used for van't Hoff).
    pub temperature_K: f64,
    /// Total hemoglobin tetramer concentration (mM).
    pub total_hb_mM: f64,
    /// Base P50 at standard conditions (mmHg).
    pub base_p50_mmHg: f64,
    /// Bohr coefficient ΔlogP50/ΔpH.
    pub bohr_coefficient: f64,
    /// 2,3-DPG effect on P50 (mmHg per mM).
    pub dpg_effect_mmHg_per_mM: f64,
    /// Reference 2,3-DPG concentration (mM).
    pub dpg_reference_mM: f64,
    /// Enthalpy of oxygenation (kcal/mol O₂).
    pub delta_h_kcal_per_mol: f64,
    /// Reference temperature for van't Hoff (K).
    pub temperature_reference_K: f64,
    /// CO₂ effect on logP50 per mmHg.
    pub co2_effect_per_mmHg: f64,
    /// Reference pCO₂ (mmHg).
    pub pco2_reference_mmHg: f64,
    /// O₂ binding rate constant (mM⁻¹·s⁻¹).
    pub k_on_per_mM_per_sec: f64,
    /// O₂ dissociation rate constant (s⁻¹).
    pub k_off_per_sec: f64,
    /// Total RBC cytoplasmic buffer capacity (slykes = mM/pH unit).
    pub buffer_capacity_slykes: f64,
    /// Reference pH at baseline lactate.
    pub ph_reference: f64,
    /// Baseline lactate concentration (mM).
    pub baseline_lactate_mM: f64,
    /// Min pH clamp.
    pub min_ph: f64,
    /// Max pH clamp.
    pub max_ph: f64,
    /// Base Adair K₁ (mmHg⁻¹). Defaults to `AdairConstants::default()`.
    pub adair_k1_per_mmHg: f64,
    /// Base Adair K₂ (mmHg⁻¹).
    pub adair_k2_per_mmHg: f64,
    /// Base Adair K₃ (mmHg⁻¹).
    pub adair_k3_per_mmHg: f64,
    /// Base Adair K₄ (mmHg⁻¹).
    pub adair_k4_per_mmHg: f64,
}

impl Default for FullBiochemBatchConfig {
    fn default() -> Self {
        let adair = AdairConstants::default();
        Self {
            dt_sec: 0.001,
            n_steps: 1000,
            external_glucose_mM: 5.0,
            atp_consumption_mM_per_sec: 0.001,
            enable_glucose_transport: true,
            enable_lactate_export: true,
            max_change_mM: 0.1,
            min_concentration_mM: 1e-9,
            membrane_tension_pN_per_nm: 0.5,
            enable_piezo1: true,
            enable_ion_homeostasis: true,
            h2o2_production_rate_mM_per_sec: 0.005,
            na_external_mM: 140.0,
            k_external_mM: 5.0,
            g_na_per_sec: 0.00024,
            g_k_per_sec: 0.00015,
            // Hb defaults match `HemoglobinSolver::default()` /
            // `AllostericParameters::default()` / `PhBufferModel::default()`.
            enable_hb: false,
            po2_mmHg: 100.0,
            pco2_mmHg: 40.0,
            temperature_K: 310.15,
            total_hb_mM: 5.0,
            base_p50_mmHg: 26.8,
            bohr_coefficient: -0.48,
            dpg_effect_mmHg_per_mM: 2.4,
            dpg_reference_mM: 5.0,
            delta_h_kcal_per_mol: -14.5,
            temperature_reference_K: 310.15,
            co2_effect_per_mmHg: 0.02,
            pco2_reference_mmHg: 40.0,
            k_on_per_mM_per_sec: 25.0,
            k_off_per_sec: 30.0,
            buffer_capacity_slykes: 60.0,
            ph_reference: 7.2,
            baseline_lactate_mM: 1.5,
            min_ph: 6.8,
            max_ph: 7.6,
            adair_k1_per_mmHg: adair.k1_per_mmHg,
            adair_k2_per_mmHg: adair.k2_per_mmHg,
            adair_k3_per_mmHg: adair.k3_per_mmHg,
            adair_k4_per_mmHg: adair.k4_per_mmHg,
        }
    }
}

/// Run a multi-cell full biochem batch on the GPU.
///
/// Operates on all 38 species of `FullyIntegratedIndices`. Each pool must
/// have at least 38 species; the kernel touches indices 0..38 in place and
/// leaves later indices (if any) untouched. The `hb_state` buffer is bound
/// but unused unless `config.enable_hb` is true (in which case use
/// [`run_full_biochem_batch_with_hb`] instead).
///
/// One dispatch advances every cell by `n_steps * dt_sec` of simulated time.
/// The pipeline + bind group are recreated per call (will be cached in a
/// `BiochemBackend` once the kernel is frozen — likely 11.2.E).
pub fn run_full_biochem_batch(
    ctx: &ComputeContext,
    pools: &mut [MetabolitePool],
    config: &FullBiochemBatchConfig,
) -> Result<()> {
    run_full_biochem_internal(ctx, pools, None, config)
}

/// Run a multi-cell full biochem batch on the GPU with the post-RK4 Hb +
/// pH update enabled. `hb_saturations` must have length equal to `pools`;
/// each entry is the per-cell hemoglobin saturation (0..1), updated in
/// place. The caller is responsible for setting `config.enable_hb = true`
/// to actually drive the Hb update — this function will set it for safety
/// regardless.
pub fn run_full_biochem_batch_with_hb(
    ctx: &ComputeContext,
    pools: &mut [MetabolitePool],
    hb_saturations: &mut [f64],
    config: &FullBiochemBatchConfig,
) -> Result<()> {
    if pools.len() != hb_saturations.len() {
        anyhow::bail!(
            "hb_saturations length ({}) must match pools length ({})",
            hb_saturations.len(),
            pools.len()
        );
    }
    let mut config = config.clone();
    config.enable_hb = true;
    run_full_biochem_internal(ctx, pools, Some(hb_saturations), &config)
}

fn run_full_biochem_internal(
    ctx: &ComputeContext,
    pools: &mut [MetabolitePool],
    hb_saturations: Option<&mut [f64]>,
    config: &FullBiochemBatchConfig,
) -> Result<()> {
    let n_cells = pools.len();
    if n_cells == 0 {
        return Ok(());
    }

    let idx = FullyIntegratedIndices::new();
    let g = &idx.glycolysis;
    let r = &idx.redox;
    let ions = &idx.ions;

    // === Upload: flatten N pools × 38 f64 → one Vec<f32> ===
    let mut state_f32 = vec![0f32; n_cells * N_SPECIES_FULL];
    for (cell_idx, pool) in pools.iter().enumerate() {
        let base = cell_idx * N_SPECIES_FULL;
        let s = &mut state_f32[base..base + N_SPECIES_FULL];
        // Glycolysis (0–16) + 2,3-BPG (17).
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
        s[17] = pool.get(idx.bisphosphoglycerate_2_3) as f32;
        // PPP (18–24).
        s[18] = pool.get(r.phosphogluconolactone_6) as f32;
        s[19] = pool.get(r.phosphogluconate_6) as f32;
        s[20] = pool.get(r.ribulose_5_phosphate) as f32;
        s[21] = pool.get(r.ribose_5_phosphate) as f32;
        s[22] = pool.get(r.xylulose_5_phosphate) as f32;
        s[23] = pool.get(r.sedoheptulose_7_phosphate) as f32;
        s[24] = pool.get(r.erythrose_4_phosphate) as f32;
        // Redox carriers (25–26).
        s[25] = pool.get(r.nadph) as f32;
        s[26] = pool.get(r.nadp_plus) as f32;
        // Glutathione (27–29).
        s[27] = pool.get(r.gsh) as f32;
        s[28] = pool.get(r.gssg) as f32;
        s[29] = pool.get(r.h2o2) as f32;
        // Ca²⁺ (30) — stored in µM units.
        s[30] = pool.get(r.ca2_plus_cytosolic) as f32;
        // GSH precursors (31–34).
        s[31] = pool.get(r.glutamate) as f32;
        s[32] = pool.get(r.cysteine) as f32;
        s[33] = pool.get(r.glycine) as f32;
        s[34] = pool.get(r.gamma_glu_cys) as f32;
        // Ions (35–37).
        s[35] = pool.get(ions.na_plus_cytosolic) as f32;
        s[36] = pool.get(ions.k_plus_cytosolic) as f32;
        s[37] = pool.get(ions.cl_minus_cytosolic) as f32;
    }

    let bytes = (state_f32.len() * std::mem::size_of::<f32>()) as u64;

    let buf_state = ctx.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("full biochem state"),
        contents: bytemuck::cast_slice(&state_f32),
        usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
    });

    // Hb-state buffer: always present (even when disabled) so the bind
    // group layout is stable across both call paths.
    let mut hb_init: Vec<f32> = vec![0.0; n_cells];
    if let Some(ref hbs) = hb_saturations {
        for (i, &s) in hbs.iter().enumerate() {
            hb_init[i] = s as f32;
        }
    }
    let buf_hb = ctx.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("full biochem hb_state"),
        contents: bytemuck::cast_slice(&hb_init),
        usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
    });
    let hb_bytes = (hb_init.len() * std::mem::size_of::<f32>()) as u64;

    let uniforms = FullBiochemUniforms {
        n_cells: n_cells as u32,
        n_steps: config.n_steps,
        dt_sec: config.dt_sec as f32,
        external_glucose_mM: config.external_glucose_mM as f32,
        atp_consumption_mM_per_sec: config.atp_consumption_mM_per_sec as f32,
        enable_glucose_transport: u32::from(config.enable_glucose_transport),
        enable_lactate_export: u32::from(config.enable_lactate_export),
        max_change_mM: config.max_change_mM as f32,
        min_concentration_mM: config.min_concentration_mM as f32,
        membrane_tension_pN_per_nm: config.membrane_tension_pN_per_nm as f32,
        enable_piezo1: u32::from(config.enable_piezo1),
        enable_ion_homeostasis: u32::from(config.enable_ion_homeostasis),
        h2o2_production_rate_mM_per_sec: config.h2o2_production_rate_mM_per_sec as f32,
        na_external_mM: config.na_external_mM as f32,
        k_external_mM: config.k_external_mM as f32,
        g_na_per_sec: config.g_na_per_sec as f32,
        g_k_per_sec: config.g_k_per_sec as f32,
        po2_mmHg: config.po2_mmHg as f32,
        pco2_mmHg: config.pco2_mmHg as f32,
        temperature_K: config.temperature_K as f32,
        total_hb_mM: config.total_hb_mM as f32,
        base_p50_mmHg: config.base_p50_mmHg as f32,
        bohr_coefficient: config.bohr_coefficient as f32,
        dpg_effect_mmHg_per_mM: config.dpg_effect_mmHg_per_mM as f32,
        dpg_reference_mM: config.dpg_reference_mM as f32,
        delta_h_kcal_per_mol: config.delta_h_kcal_per_mol as f32,
        temperature_reference_K: config.temperature_reference_K as f32,
        co2_effect_per_mmHg: config.co2_effect_per_mmHg as f32,
        pco2_reference_mmHg: config.pco2_reference_mmHg as f32,
        k_on_per_mM_per_sec: config.k_on_per_mM_per_sec as f32,
        k_off_per_sec: config.k_off_per_sec as f32,
        buffer_capacity_slykes: config.buffer_capacity_slykes as f32,
        ph_reference: config.ph_reference as f32,
        baseline_lactate_mM: config.baseline_lactate_mM as f32,
        min_ph: config.min_ph as f32,
        max_ph: config.max_ph as f32,
        enable_hb: u32::from(config.enable_hb),
        adair_k1: config.adair_k1_per_mmHg as f32,
        adair_k2: config.adair_k2_per_mmHg as f32,
        adair_k3: config.adair_k3_per_mmHg as f32,
        adair_k4: config.adair_k4_per_mmHg as f32,
        _pad0: 0.0,
        _pad1: 0.0,
        _pad2: 0.0,
    };
    let buf_u = ctx.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("full biochem uniforms"),
        contents: bytemuck::bytes_of(&uniforms),
        usage: wgpu::BufferUsages::UNIFORM,
    });

    let buf_staging = ctx.device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("full biochem staging"),
        size: bytes,
        usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
        mapped_at_creation: false,
    });

    let buf_hb_staging = ctx.device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("full biochem hb staging"),
        size: hb_bytes,
        usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
        mapped_at_creation: false,
    });

    let bgl = ctx.device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
        label: Some("full biochem bgl"),
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
            wgpu::BindGroupLayoutEntry {
                binding: 2,
                visibility: wgpu::ShaderStages::COMPUTE,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Storage { read_only: false },
                    has_dynamic_offset: false,
                    min_binding_size: None,
                },
                count: None,
            },
        ],
    });

    let bg = ctx.device.create_bind_group(&wgpu::BindGroupDescriptor {
        label: Some("full biochem bg"),
        layout: &bgl,
        entries: &[
            wgpu::BindGroupEntry { binding: 0, resource: buf_state.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 1, resource: buf_u.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 2, resource: buf_hb.as_entire_binding() },
        ],
    });

    let shader = ctx.device.create_shader_module(wgpu::ShaderModuleDescriptor {
        label: Some("full biochem shader"),
        source: wgpu::ShaderSource::Wgsl(KERNEL_FULL_WGSL.into()),
    });

    let pipeline_layout = ctx.device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
        label: Some("full biochem pl"),
        bind_group_layouts: &[&bgl],
        push_constant_ranges: &[],
    });

    let pipeline = ctx.device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
        label: Some("full biochem pipeline"),
        layout: Some(&pipeline_layout),
        module: &shader,
        entry_point: "rk4_step",
        compilation_options: Default::default(),
    });

    let mut encoder = ctx.device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
        label: Some("full biochem encoder"),
    });

    {
        let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
            label: Some("full biochem pass"),
            timestamp_writes: None,
        });
        pass.set_pipeline(&pipeline);
        pass.set_bind_group(0, &bg, &[]);
        let groups = ((n_cells as u32) + 63) / 64;
        pass.dispatch_workgroups(groups, 1, 1);
    }

    encoder.copy_buffer_to_buffer(&buf_state, 0, &buf_staging, 0, bytes);
    encoder.copy_buffer_to_buffer(&buf_hb, 0, &buf_hb_staging, 0, hb_bytes);
    ctx.queue.submit(std::iter::once(encoder.finish()));

    // === Download state ===
    let slice = buf_staging.slice(..);
    let (tx, rx) = mpsc::channel();
    slice.map_async(wgpu::MapMode::Read, move |result| {
        let _ = tx.send(result);
    });

    // Map Hb buffer in parallel.
    let hb_slice = buf_hb_staging.slice(..);
    let (tx_hb, rx_hb) = mpsc::channel();
    hb_slice.map_async(wgpu::MapMode::Read, move |result| {
        let _ = tx_hb.send(result);
    });

    ctx.device.poll(wgpu::Maintain::Wait);
    rx.recv()
        .context("full biochem: map_async channel closed")?
        .map_err(|e| anyhow::anyhow!("full biochem: map_async failed: {e:?}"))?;
    rx_hb.recv()
        .context("full biochem: hb map_async channel closed")?
        .map_err(|e| anyhow::anyhow!("full biochem: hb map_async failed: {e:?}"))?;

    let mapped = slice.get_mapped_range();
    let result_f32: &[f32] = bytemuck::cast_slice(&mapped);
    let mapped_hb = hb_slice.get_mapped_range();
    let hb_f32: &[f32] = bytemuck::cast_slice(&mapped_hb);

    if let Some(hbs) = hb_saturations {
        for (i, slot) in hbs.iter_mut().enumerate() {
            *slot = hb_f32[i] as f64;
        }
    }

    // === Distribute back ===
    for (cell_idx, pool) in pools.iter_mut().enumerate() {
        let base = cell_idx * N_SPECIES_FULL;
        let s = &result_f32[base..base + N_SPECIES_FULL];
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
        pool.set(idx.bisphosphoglycerate_2_3, s[17] as f64);
        pool.set(r.phosphogluconolactone_6, s[18] as f64);
        pool.set(r.phosphogluconate_6, s[19] as f64);
        pool.set(r.ribulose_5_phosphate, s[20] as f64);
        pool.set(r.ribose_5_phosphate, s[21] as f64);
        pool.set(r.xylulose_5_phosphate, s[22] as f64);
        pool.set(r.sedoheptulose_7_phosphate, s[23] as f64);
        pool.set(r.erythrose_4_phosphate, s[24] as f64);
        pool.set(r.nadph, s[25] as f64);
        pool.set(r.nadp_plus, s[26] as f64);
        pool.set(r.gsh, s[27] as f64);
        pool.set(r.gssg, s[28] as f64);
        pool.set(r.h2o2, s[29] as f64);
        pool.set(r.ca2_plus_cytosolic, s[30] as f64);
        pool.set(r.glutamate, s[31] as f64);
        pool.set(r.cysteine, s[32] as f64);
        pool.set(r.glycine, s[33] as f64);
        pool.set(r.gamma_glu_cys, s[34] as f64);
        pool.set(ions.na_plus_cytosolic, s[35] as f64);
        pool.set(ions.k_plus_cytosolic, s[36] as f64);
        pool.set(ions.cl_minus_cytosolic, s[37] as f64);
    }

    drop(mapped);
    drop(mapped_hb);
    buf_staging.unmap();
    buf_hb_staging.unmap();

    Ok(())
}
