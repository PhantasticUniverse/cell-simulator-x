//! Membrane tension computation from physics state.
//!
//! Computes global membrane tension (pN/nm) from the Skalak strain invariants
//! computed by the physics solver. This tension value drives Piezo1
//! mechanotransduction in the biochemistry module.
//!
//! ## Algorithm
//! Tension is computed from the area-weighted average of strain invariants:
//! - T = Gs × (|I₁| + |I₂|) / 2
//!
//! Where:
//! - Gs = 5.5 μN/m = 5.5 pN/nm (shear modulus)
//! - I₁ = α² + β² - 2 (first strain invariant)
//! - I₂ = α²β² - 1 (second strain invariant, area change)
//!
//! ## Temporal Averaging
//! A rolling average over recent values provides stability:
//! - Prevents rapid tension oscillations from coupling instability
//! - Smooths out numerical noise from physics integration
//!
//! ## References
//! - Skalak et al., Biophys J 1973
//! - Evans & Waugh, Biophys J 1977

use std::collections::VecDeque;

use crate::geometry::Mesh;
use crate::physics::SkalakSolver;

/// Computes membrane tension from physics state for Piezo1 activation.
pub struct TensionComputer {
    /// Number of samples for temporal averaging
    averaging_window: usize,
    /// History of tension values for averaging
    tension_history: VecDeque<f64>,
    /// Shear modulus in pN/nm (converted from μN/m)
    shear_modulus_pN_per_nm: f64,
}

impl TensionComputer {
    /// Create a new tension computer.
    ///
    /// # Arguments
    /// * `averaging_window` - Number of samples for temporal averaging (default: 10)
    pub fn new(averaging_window: usize) -> Self {
        Self {
            averaging_window,
            tension_history: VecDeque::with_capacity(averaging_window),
            shear_modulus_pN_per_nm: 5.5, // 5.5 μN/m = 5.5 pN/nm (same numerical value)
        }
    }

    /// Compute global membrane tension from strain invariants.
    ///
    /// Returns tension in pN/nm, suitable for Piezo1 activation.
    ///
    /// # Arguments
    /// * `skalak_solver` - The membrane mechanics solver with element data
    /// * `mesh` - Current mesh with deformed positions
    ///
    /// # Algorithm
    /// 1. Compute strain invariants (I₁, I₂) for each membrane element
    /// 2. Area-weight the contributions
    /// 3. Apply shear modulus to get tension
    /// 4. Temporally average for stability
    pub fn compute_global_tension_pN_per_nm(
        &mut self,
        skalak_solver: &SkalakSolver,
        mesh: &Mesh,
    ) -> f64 {
        let mut total_weighted_strain = 0.0;
        let mut total_area = 0.0;

        // Compute area-weighted average of strain invariants
        for element in &skalak_solver.elements {
            let i0 = element.vertex_indices[0] as usize;
            let i1 = element.vertex_indices[1] as usize;
            let i2 = element.vertex_indices[2] as usize;

            let positions = [
                mesh.vertices[i0].position_vec3(),
                mesh.vertices[i1].position_vec3(),
                mesh.vertices[i2].position_vec3(),
            ];

            // Get strain invariants
            let (i1, i2) = element.strain_invariants(positions);

            // Current area for weighting
            let area = element.current_area(positions);

            // Use absolute values of strain invariants (tension can come from both)
            // Weighted by current element area
            let strain_magnitude = (i1.abs() + i2.abs()) as f64;
            total_weighted_strain += strain_magnitude * area as f64;
            total_area += area as f64;
        }

        // Compute average strain magnitude
        let avg_strain = if total_area > 1e-12 {
            total_weighted_strain / total_area
        } else {
            0.0
        };

        // Convert to tension: T = Gs × strain / 2
        // The factor of 2 accounts for averaging I₁ + I₂
        let instantaneous_tension = self.shear_modulus_pN_per_nm * avg_strain / 2.0;

        // Add to history for temporal averaging
        self.tension_history.push_back(instantaneous_tension);
        if self.tension_history.len() > self.averaging_window {
            self.tension_history.pop_front();
        }

        // Return averaged tension
        self.averaged_tension()
    }

    /// Get the temporally averaged tension.
    pub fn averaged_tension(&self) -> f64 {
        if self.tension_history.is_empty() {
            return 0.0;
        }
        self.tension_history.iter().sum::<f64>() / self.tension_history.len() as f64
    }

    /// Get the most recent instantaneous tension (before averaging).
    pub fn instantaneous_tension(&self) -> f64 {
        self.tension_history.back().copied().unwrap_or(0.0)
    }

    /// Clear the tension history (e.g., when resetting simulation).
    pub fn reset(&mut self) {
        self.tension_history.clear();
    }

    /// Set the shear modulus (useful for testing or disease models).
    pub fn set_shear_modulus(&mut self, gs_pN_per_nm: f64) {
        self.shear_modulus_pN_per_nm = gs_pN_per_nm;
    }
}

impl Default for TensionComputer {
    fn default() -> Self {
        Self::new(10) // 10-sample averaging window
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::GeometryParameters;
    use crate::physics::SkalakMaterial;

    fn create_test_mesh() -> Mesh {
        let params = GeometryParameters {
            cell_radius_um: 3.91,
            fung_tong_c0_um: 0.81,
            fung_tong_c2_um: 7.83,
            fung_tong_c4_um: -4.39,
            mesh_resolution: 8, // Small for testing
            spectrin_target_count: 50,
        };
        Mesh::generate_rbc(&params)
    }

    #[test]
    fn test_tension_at_rest() {
        let mesh = create_test_mesh();
        let skalak = SkalakSolver::new(&mesh, SkalakMaterial::default());
        let mut computer = TensionComputer::default();

        // At rest, strain invariants should be ~0, so tension ~0
        let tension = computer.compute_global_tension_pN_per_nm(&skalak, &mesh);

        // Tension should be very small at equilibrium
        assert!(
            tension < 0.5,
            "Tension at rest should be near zero, got: {:.4} pN/nm",
            tension
        );
    }

    #[test]
    fn test_averaging() {
        let mut computer = TensionComputer::new(5);

        // Manually add history values
        for i in 1..=5 {
            computer.tension_history.push_back(i as f64);
        }

        // Average of 1,2,3,4,5 = 3.0
        assert!((computer.averaged_tension() - 3.0).abs() < 0.01);
    }

    #[test]
    fn test_reset() {
        let mut computer = TensionComputer::default();
        computer.tension_history.push_back(1.0);
        computer.tension_history.push_back(2.0);

        computer.reset();
        assert!(computer.tension_history.is_empty());
        assert!(computer.averaged_tension().abs() < 0.01);
    }
}
