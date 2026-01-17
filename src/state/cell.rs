//! Cell state data structures.
//!
//! Primary structures for representing the geometric and overall state of the RBC.

use crate::{
    config::Parameters,
    geometry::{Mesh, SpectrinNetwork},
};
use glam::Vec3;

/// Complete state of a single red blood cell
pub struct CellState {
    /// Geometric representation
    pub geometry: GeometryState,
    /// Position in world space (μm)
    pub position_um: Vec3,
    /// Velocity (μm/s)
    pub velocity_um_per_sec: Vec3,
    /// Angular velocity (rad/s)
    pub angular_velocity_rad_per_sec: Vec3,
    /// Cell volume (μm³)
    /// Reference: Normal human RBC volume ~90 fL = 90 μm³
    /// Source: Mohandas & Gallagher, Blood 2008
    pub volume_um3: f32,
    /// Surface area (μm²)
    /// Reference: Normal human RBC surface area ~135 μm²
    /// Source: Evans & Fung, Microvasc Res 1972
    pub surface_area_um2: f32,
}

impl CellState {
    /// Create a new cell state from parameters
    pub fn new(params: &Parameters) -> Self {
        let mesh = Mesh::generate_rbc(&params.geometry);
        let spectrin_network = SpectrinNetwork::generate(&mesh, &params.geometry);

        // Calculate volume and surface area from mesh
        let volume_um3 = mesh.calculate_volume();
        let surface_area_um2 = mesh.calculate_surface_area();

        Self {
            geometry: GeometryState {
                mesh,
                spectrin_network,
            },
            position_um: Vec3::ZERO,
            velocity_um_per_sec: Vec3::ZERO,
            angular_velocity_rad_per_sec: Vec3::ZERO,
            volume_um3,
            surface_area_um2,
        }
    }
}

/// Geometric state containing mesh and network data
pub struct GeometryState {
    /// Surface mesh of the RBC (biconcave disc)
    pub mesh: Mesh,
    /// Spectrin-actin cytoskeleton network
    pub spectrin_network: SpectrinNetwork,
}
