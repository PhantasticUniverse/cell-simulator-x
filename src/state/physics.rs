//! Physics state data structures.
//!
//! Represents the mechanical state of the red blood cell membrane
//! and cytoskeleton.

use glam::Vec3;

/// Physics/mechanics state of the cell
#[derive(Debug, Clone)]
pub struct PhysicsState {
    /// Membrane tension state
    pub membrane: MembraneState,
    /// Per-vertex forces (μN)
    pub vertex_forces_uN: Vec<Vec3>,
    /// Per-vertex velocities (μm/s)
    pub vertex_velocities_um_per_sec: Vec<Vec3>,
    /// Total elastic energy of the membrane (pJ)
    pub elastic_energy_pJ: f32,
}

impl PhysicsState {
    /// Create a new physics state for a mesh with given vertex count
    pub fn new(vertex_count: usize) -> Self {
        Self {
            membrane: MembraneState::default(),
            vertex_forces_uN: vec![Vec3::ZERO; vertex_count],
            vertex_velocities_um_per_sec: vec![Vec3::ZERO; vertex_count],
            elastic_energy_pJ: 0.0,
        }
    }
}

/// Membrane mechanical state
#[derive(Debug, Clone)]
pub struct MembraneState {
    /// Shear modulus (μN/m)
    /// Reference: 5.5 ± 1.8 μN/m
    /// Source: Evans & Waugh, Biophys J 1977; Henon et al., Biophys J 1999
    pub shear_modulus_uN_per_m: f32,
    /// Bending modulus (pN·μm = 10⁻¹⁸ J)
    /// Reference: 1.8 × 10⁻¹⁹ J ≈ 44 kT at 37°C
    /// Source: Evans, Biophys J 1983
    pub bending_modulus_pN_um: f32,
    /// Area expansion modulus (mN/m)
    /// Reference: ~400-500 mN/m
    /// Source: Evans et al., Biophys J 1976
    pub area_modulus_mN_per_m: f32,
    /// Current area strain (ΔA/A₀)
    pub area_strain: f32,
    /// Mean curvature (1/μm)
    pub mean_curvature_per_um: f32,
}

impl Default for MembraneState {
    fn default() -> Self {
        Self {
            shear_modulus_uN_per_m: 5.5,
            bending_modulus_pN_um: 0.18, // ~44 kT
            area_modulus_mN_per_m: 450.0,
            area_strain: 0.0,
            mean_curvature_per_um: 0.0,
        }
    }
}
