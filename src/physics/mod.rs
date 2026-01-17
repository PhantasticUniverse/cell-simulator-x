//! Physics simulation module for RBC mechanics.
//!
//! This module implements:
//! - DPD (Dissipative Particle Dynamics) fluid solver
//! - Skalak membrane strain energy model
//! - WLC (Worm-Like Chain) spectrin elasticity
//! - Velocity-Verlet time integration
//!
//! References:
//! - DPD: Groot & Warren, J Chem Phys 1997
//! - Skalak: Skalak et al., Biophys J 1973
//! - WLC: Marko & Siggia, Macromolecules 1995

pub mod dpd;
pub mod integrator;
pub mod membrane;
pub mod wlc;

pub use dpd::{DPDParameters, DPDParticle, DPDSolver};
pub use integrator::{IntegratorState, VelocityVerlet};
pub use membrane::{MembraneElement, SkalakMaterial, SkalakSolver};
pub use wlc::{WLCParameters, WLCSolver};

use glam::Vec3;

use crate::geometry::{Mesh, SpectrinNetwork};
use crate::state::PhysicsState;

/// Configuration parameters for the physics simulation
#[derive(Debug, Clone)]
pub struct PhysicsConfig {
    /// Timestep in seconds
    pub dt_sec: f32,
    /// Temperature in Kelvin
    pub temperature_K: f32,
    /// Enable DPD thermal fluctuations
    pub enable_thermal_noise: bool,
    /// Damping coefficient for membrane vertices
    pub membrane_damping: f32,
}

impl Default for PhysicsConfig {
    fn default() -> Self {
        Self {
            dt_sec: 1e-6,      // 1 microsecond timestep
            temperature_K: 310.0, // 37°C body temperature
            enable_thermal_noise: true,
            membrane_damping: 0.1,
        }
    }
}

/// Main physics solver combining all force calculations
pub struct PhysicsSolver {
    /// WLC solver for spectrin elasticity
    pub wlc_solver: WLCSolver,
    /// Skalak solver for membrane mechanics
    pub skalak_solver: SkalakSolver,
    /// DPD solver for fluid dynamics
    pub dpd_solver: DPDSolver,
    /// Time integrator
    pub integrator: VelocityVerlet,
    /// Configuration
    pub config: PhysicsConfig,
    /// Simulation time in seconds
    pub time_sec: f32,
}

impl PhysicsSolver {
    /// Create a new physics solver
    pub fn new(mesh: &Mesh, config: PhysicsConfig) -> Self {
        let wlc_solver = WLCSolver::new(WLCParameters::default());
        let skalak_solver = SkalakSolver::new(mesh, SkalakMaterial::default());
        let dpd_solver = DPDSolver::new(DPDParameters::default());
        let integrator = VelocityVerlet::new();

        Self {
            wlc_solver,
            skalak_solver,
            dpd_solver,
            integrator,
            config,
            time_sec: 0.0,
        }
    }

    /// Perform one physics timestep
    ///
    /// Updates vertex positions and velocities based on:
    /// 1. WLC spectrin forces
    /// 2. Skalak membrane forces
    /// 3. DPD dissipative and random forces
    pub fn step(
        &mut self,
        mesh: &mut Mesh,
        spectrin: &SpectrinNetwork,
        physics_state: &mut PhysicsState,
    ) {
        let dt = self.config.dt_sec;
        let n_vertices = mesh.vertices.len();

        // Ensure force and velocity vectors are sized correctly
        if physics_state.vertex_forces_uN.len() != n_vertices {
            physics_state.vertex_forces_uN = vec![Vec3::ZERO; n_vertices];
        }
        if physics_state.external_forces_uN.len() != n_vertices {
            physics_state.external_forces_uN = vec![Vec3::ZERO; n_vertices];
        }
        if physics_state.vertex_velocities_um_per_sec.len() != n_vertices {
            physics_state.vertex_velocities_um_per_sec = vec![Vec3::ZERO; n_vertices];
        }

        // For first step, initialize forces from external forces
        // (subsequent steps use forces computed at end of previous step)
        if physics_state.step_count == 0 {
            for i in 0..n_vertices {
                physics_state.vertex_forces_uN[i] = physics_state.external_forces_uN[i];
            }
        }

        // Get current positions
        let positions: Vec<Vec3> = mesh.vertices.iter().map(|v| v.position_vec3()).collect();

        // === Velocity-Verlet: First half-step for velocities ===
        self.integrator.half_step_velocity(
            &mut physics_state.vertex_velocities_um_per_sec,
            &physics_state.vertex_forces_uN,
            dt,
        );

        // === Update positions ===
        let new_positions = self.integrator.step_position(
            &positions,
            &physics_state.vertex_velocities_um_per_sec,
            dt,
        );

        // Update mesh vertices with new positions
        for (i, pos) in new_positions.iter().enumerate() {
            mesh.vertices[i].position = pos.to_array();
        }

        // Recompute normals after position update
        self.recompute_normals(mesh);

        // === Compute forces at new positions ===
        // Start with external forces (persistent)
        for (i, force) in physics_state.vertex_forces_uN.iter_mut().enumerate() {
            *force = physics_state.external_forces_uN.get(i).copied().unwrap_or(Vec3::ZERO);
        }

        // 1. WLC spectrin forces
        let spectrin_forces = self.wlc_solver.compute_forces(spectrin, &new_positions);
        for (node_idx, force) in spectrin_forces {
            // Map spectrin node to nearest mesh vertex
            if let Some(node) = spectrin.nodes.get(node_idx) {
                let vertex_idx = node.mesh_vertex_idx as usize;
                if vertex_idx < n_vertices {
                    physics_state.vertex_forces_uN[vertex_idx] += force;
                }
            }
        }

        // 2. Skalak membrane forces
        let (membrane_forces, elastic_energy) = self.skalak_solver.compute_forces(mesh);
        for (i, force) in membrane_forces.into_iter().enumerate() {
            if i < n_vertices {
                physics_state.vertex_forces_uN[i] += force;
            }
        }
        physics_state.elastic_energy_pJ = elastic_energy;

        // 3. DPD dissipative and random forces
        if self.config.enable_thermal_noise {
            let dpd_forces = self.dpd_solver.compute_membrane_forces(
                &new_positions,
                &physics_state.vertex_velocities_um_per_sec,
                self.config.temperature_K,
            );
            for (i, force) in dpd_forces.into_iter().enumerate() {
                if i < n_vertices {
                    physics_state.vertex_forces_uN[i] += force;
                }
            }
        }

        // 4. Apply damping
        for (i, vel) in physics_state.vertex_velocities_um_per_sec.iter().enumerate() {
            physics_state.vertex_forces_uN[i] -= *vel * self.config.membrane_damping;
        }

        // === Velocity-Verlet: Second half-step for velocities ===
        self.integrator.half_step_velocity(
            &mut physics_state.vertex_velocities_um_per_sec,
            &physics_state.vertex_forces_uN,
            dt,
        );

        // Update simulation time
        self.time_sec += dt;
        physics_state.simulation_time_sec += dt as f64;
        physics_state.step_count += 1;
    }

    /// Recompute vertex normals after position update
    fn recompute_normals(&self, mesh: &mut Mesh) {
        let n_vertices = mesh.vertices.len();
        let mut normals = vec![Vec3::ZERO; n_vertices];

        // Accumulate face normals for each vertex
        for chunk in mesh.indices.chunks(3) {
            let i0 = chunk[0] as usize;
            let i1 = chunk[1] as usize;
            let i2 = chunk[2] as usize;

            let v0 = mesh.vertices[i0].position_vec3();
            let v1 = mesh.vertices[i1].position_vec3();
            let v2 = mesh.vertices[i2].position_vec3();

            let edge1 = v1 - v0;
            let edge2 = v2 - v0;
            let face_normal = edge1.cross(edge2);

            normals[i0] += face_normal;
            normals[i1] += face_normal;
            normals[i2] += face_normal;
        }

        // Normalize and update
        for (i, normal) in normals.into_iter().enumerate() {
            let normalized = if normal.length_squared() > 0.0 {
                normal.normalize()
            } else {
                Vec3::Z
            };
            mesh.vertices[i].normal = normalized.to_array();
        }
    }

    /// Apply an external force to a specific vertex (persists until cleared)
    pub fn apply_external_force(&self, physics_state: &mut PhysicsState, vertex_idx: usize, force_uN: Vec3) {
        physics_state.set_external_force(vertex_idx, force_uN);
    }

    /// Clear all external forces
    pub fn clear_external_forces(&self, physics_state: &mut PhysicsState) {
        physics_state.clear_external_forces();
    }

    /// Get total kinetic energy of the system
    pub fn kinetic_energy_pJ(&self, physics_state: &PhysicsState) -> f32 {
        // Assume mass per vertex is ~1 pg (10^-15 kg) for simplicity
        // KE = 0.5 * m * v^2
        let mass_pg = 1.0; // picograms
        physics_state
            .vertex_velocities_um_per_sec
            .iter()
            .map(|v| 0.5 * mass_pg * v.length_squared())
            .sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::GeometryParameters;

    #[test]
    fn test_physics_solver_creation() {
        let params = GeometryParameters {
            cell_radius_um: 3.91,
            fung_tong_c0_um: 0.81,
            fung_tong_c2_um: 7.83,
            fung_tong_c4_um: -4.39,
            mesh_resolution: 10, // Small for testing
            spectrin_target_count: 100,
        };
        let mesh = Mesh::generate_rbc(&params);
        let config = PhysicsConfig::default();

        let solver = PhysicsSolver::new(&mesh, config);
        assert!(solver.time_sec == 0.0);
    }
}
