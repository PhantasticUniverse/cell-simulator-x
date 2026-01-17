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
    /// Per-vertex forces (μN) - computed each step
    pub vertex_forces_uN: Vec<Vec3>,
    /// Per-vertex external forces (μN) - persistent until cleared
    pub external_forces_uN: Vec<Vec3>,
    /// Per-vertex velocities (μm/s)
    pub vertex_velocities_um_per_sec: Vec<Vec3>,
    /// Total elastic energy of the membrane (pJ)
    pub elastic_energy_pJ: f32,
    /// Spectrin network elastic energy (pJ)
    pub spectrin_energy_pJ: f32,
    /// Total kinetic energy (pJ)
    pub kinetic_energy_pJ: f32,
    /// Simulation time in seconds
    pub simulation_time_sec: f64,
    /// Number of physics steps taken
    pub step_count: u64,
    /// Reference positions for strain calculation
    pub reference_positions: Vec<Vec3>,
}

impl PhysicsState {
    /// Create a new physics state for a mesh with given vertex count
    pub fn new(vertex_count: usize) -> Self {
        Self {
            membrane: MembraneState::default(),
            vertex_forces_uN: vec![Vec3::ZERO; vertex_count],
            external_forces_uN: vec![Vec3::ZERO; vertex_count],
            vertex_velocities_um_per_sec: vec![Vec3::ZERO; vertex_count],
            elastic_energy_pJ: 0.0,
            spectrin_energy_pJ: 0.0,
            kinetic_energy_pJ: 0.0,
            simulation_time_sec: 0.0,
            step_count: 0,
            reference_positions: Vec::new(),
        }
    }

    /// Set external force on a vertex (persists until cleared)
    pub fn set_external_force(&mut self, vertex_idx: usize, force: Vec3) {
        if vertex_idx < self.external_forces_uN.len() {
            self.external_forces_uN[vertex_idx] = force;
        }
    }

    /// Clear all external forces
    pub fn clear_external_forces(&mut self) {
        for f in self.external_forces_uN.iter_mut() {
            *f = Vec3::ZERO;
        }
    }

    /// Initialize reference positions from current mesh
    pub fn init_reference_positions(&mut self, positions: &[Vec3]) {
        self.reference_positions = positions.to_vec();
    }

    /// Get total mechanical energy (elastic + kinetic)
    pub fn total_energy_pJ(&self) -> f32 {
        self.elastic_energy_pJ + self.spectrin_energy_pJ + self.kinetic_energy_pJ
    }

    /// Reset all velocities to zero
    pub fn reset_velocities(&mut self) {
        for v in self.vertex_velocities_um_per_sec.iter_mut() {
            *v = Vec3::ZERO;
        }
    }

    /// Get maximum velocity magnitude
    pub fn max_velocity_um_per_sec(&self) -> f32 {
        self.vertex_velocities_um_per_sec
            .iter()
            .map(|v| v.length())
            .fold(0.0f32, f32::max)
    }

    /// Get maximum force magnitude
    pub fn max_force_uN(&self) -> f32 {
        self.vertex_forces_uN
            .iter()
            .map(|f| f.length())
            .fold(0.0f32, f32::max)
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
