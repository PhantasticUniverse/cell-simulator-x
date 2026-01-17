//! Rendering module using wgpu (WebGPU/Metal backend).
//!
//! Provides GPU-accelerated rendering of the RBC mesh and spectrin network.

mod camera;
mod pipeline;

pub use camera::Camera;
pub use pipeline::RenderState;
