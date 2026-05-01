//! GPU compute scaffolding (Phase 11.0).
//!
//! `ComputeContext` owns a wgpu device + queue dedicated to compute work.
//! It is intentionally separate from [`crate::render::RenderState`]: when
//! the simulator runs headless (CLI diagnostics, batch jobs, tests), there
//! is no window and no render surface, but compute still needs to run.
//!
//! ## Sub-phase scope (11.0)
//!
//! Phase 11.0 is the compute scaffolding deliverable. It establishes:
//!
//! - A pattern for headless wgpu device creation (no surface required).
//! - The `pollster::block_on` boundary between async wgpu APIs and
//!   synchronous CLI / test contexts.
//! - A reusable bind-group / pipeline-layout / dispatch shape that later
//!   sub-phases (11.2 biochem ODE, 11.3 mechanics) will follow.
//! - A first compute kernel (`vec_add`) used as a CPU↔GPU correctness
//!   sentinel in [`diagnostics`].
//!
//! ## Out of scope for 11.0
//!
//! - Sharing a `Device`/`Queue` with `RenderState` (deferred to 11.4).
//! - Compute kernels for the actual simulation (that's 11.2 and 11.3).
//! - Timestamp queries / `wgpu-profiler` integration (deferred to 11.5).
//! - Any wgpu version upgrade (the pre-flight on 2026-05-01 found that
//!   wgpu 0.20 → 29.x drags egui 0.28 → 0.34 with five distinct API
//!   breakages; the user opted to accept HUD freeze if needed and we
//!   chose to stay on 0.20 since its compute support is already complete).

use std::sync::Arc;

use anyhow::{Context as _, Result};

pub mod biochem;
pub mod diagnostics;
pub mod physics;

pub use biochem::{
    run_full_biochem_batch, run_full_biochem_batch_with_hb, run_glycolysis_batch,
    FullBiochemBatchConfig, GlycolysisBatchConfig,
};
pub use diagnostics::{run_diagnose_gpu, vec_add_gpu};
pub use physics::{
    run_dpd_membrane_forces, run_skalak_forces, run_verlet_step, run_wlc_forces, EdgeGpuPub,
    ElementGpuPub, SkalakBackendData, VerletConfig, WlcBackendData,
};

/// Headless wgpu compute context.
///
/// Owns its own `Device` and `Queue` and does not participate in any
/// surface. Multiple [`ComputeContext`]s can coexist with a
/// [`crate::render::RenderState`] but they will not share buffers
/// (that's the 11.4 deliverable).
pub struct ComputeContext {
    /// The adapter used to create the device. Held so the device's
    /// lifetime is tied to it and so callers can interrogate adapter
    /// info for diagnostics.
    pub adapter: wgpu::Adapter,
    /// The compute device. `Arc` so it can be borrowed by kernel
    /// helpers without forcing them to take a `&self` lifetime.
    pub device: Arc<wgpu::Device>,
    /// The command queue.
    pub queue: Arc<wgpu::Queue>,
}

impl ComputeContext {
    /// Create a headless compute context.
    ///
    /// Requests a high-performance adapter with `compatible_surface: None`
    /// so no window is required. Suitable for CLI diagnostics, tests,
    /// and batch jobs.
    pub async fn new_headless() -> Result<Self> {
        let instance = wgpu::Instance::new(wgpu::InstanceDescriptor {
            backends: wgpu::Backends::all(),
            ..Default::default()
        });

        let adapter = instance
            .request_adapter(&wgpu::RequestAdapterOptions {
                power_preference: wgpu::PowerPreference::HighPerformance,
                compatible_surface: None,
                force_fallback_adapter: false,
            })
            .await
            .context("failed to find a suitable GPU adapter for compute")?;

        log::info!("Compute adapter: {:?}", adapter.get_info());

        let (device, queue) = adapter
            .request_device(
                &wgpu::DeviceDescriptor {
                    label: Some("Cell Simulator X Compute Device"),
                    required_features: wgpu::Features::empty(),
                    required_limits: wgpu::Limits::default(),
                },
                None,
            )
            .await
            .context("failed to request compute device")?;

        Ok(Self {
            adapter,
            device: Arc::new(device),
            queue: Arc::new(queue),
        })
    }

    /// Synchronous wrapper around [`new_headless`].
    pub fn new_headless_blocking() -> Result<Self> {
        pollster::block_on(Self::new_headless())
    }

    /// Pretty-print adapter info.
    pub fn adapter_summary(&self) -> String {
        let info = self.adapter.get_info();
        format!(
            "{} ({:?} on {:?})",
            info.name, info.device_type, info.backend
        )
    }
}
