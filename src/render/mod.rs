//! Rendering module using wgpu (WebGPU/Metal backend).
//!
//! Provides GPU-accelerated rendering of the RBC mesh and spectrin network.

mod camera;
pub mod hud;
mod pipeline;

pub use camera::Camera;
pub use hud::{ExportAction, HudColors, HudOverlay, HudState, HudTheme};
pub use pipeline::RenderState;
