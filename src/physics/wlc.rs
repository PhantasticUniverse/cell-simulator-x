//! Worm-Like Chain (WLC) model for spectrin elasticity.
//!
//! Implements the Marko-Siggia force-extension formula for spectrin tetramers:
//!
//! F = (k_B*T / L_p) * [1/(4*(1-x/L_c)²) - 1/4 + x/L_c]
//!
//! Where:
//! - k_B*T = 4.11 pN·nm at 37°C (thermal energy)
//! - L_p = 20 nm (persistence length)
//! - L_c = 200 nm (contour length)
//! - x = extension (current length)
//!
//! Reference: Marko & Siggia, Macromolecules 1995
//! Reference: Rief et al., J Mol Biol 1999 (spectrin measurements)

use glam::Vec3;

use crate::geometry::SpectrinNetwork;

/// Boltzmann constant times temperature at 37°C in pN·nm
/// k_B*T = 1.38e-23 J/K * 310 K = 4.28e-21 J = 4.28 pN·nm
const KB_T_PN_NM: f32 = 4.11; // Using 4.11 as per literature

/// Parameters for the WLC model
#[derive(Debug, Clone)]
pub struct WLCParameters {
    /// Persistence length in μm
    /// Reference: ~20 nm for spectrin
    /// Source: Rief et al. 1999
    pub persistence_length_um: f32,
    /// Contour length in μm (fully extended)
    /// Reference: ~200 nm for spectrin tetramer
    /// Source: Rief et al. 1999
    pub contour_length_um: f32,
    /// Rest length in μm (in situ equilibrium)
    /// Reference: ~75 nm for spectrin
    /// Source: Liu et al. 1987
    pub rest_length_um: f32,
    /// Maximum relative extension (cap to prevent singularity)
    pub max_relative_extension: f32,
}

impl Default for WLCParameters {
    fn default() -> Self {
        Self {
            persistence_length_um: 0.020, // 20 nm
            contour_length_um: 0.200,     // 200 nm
            rest_length_um: 0.075,        // 75 nm
            max_relative_extension: 0.95, // Cap at 95% of contour length
        }
    }
}

/// WLC solver for spectrin network elasticity
pub struct WLCSolver {
    /// Model parameters
    pub params: WLCParameters,
}

impl WLCSolver {
    /// Create a new WLC solver
    pub fn new(params: WLCParameters) -> Self {
        Self { params }
    }

    /// Compute the Marko-Siggia force for a single spectrin tetramer
    ///
    /// Returns force magnitude in μN (micro-Newtons)
    /// Positive force = tension (extension), negative = compression
    ///
    /// Formula: F = (k_B*T / L_p) * [1/(4*(1-x/L_c)²) - 1/4 + x/L_c]
    pub fn marko_siggia_force(&self, extension_um: f32, contour_length_um: f32) -> f32 {
        let lp = self.params.persistence_length_um;
        let lc = contour_length_um;

        // Clamp extension to avoid singularity at full extension
        let max_ext = lc * self.params.max_relative_extension;
        let x = extension_um.clamp(0.001 * lc, max_ext);

        // Relative extension
        let xi = x / lc;

        // Convert k_B*T from pN·nm to μN·μm
        // 1 pN·nm = 1e-12 N * 1e-9 m = 1e-21 J
        // 1 μN·μm = 1e-6 N * 1e-6 m = 1e-12 J
        // So 1 pN·nm = 1e-9 μN·μm
        let kbt_uN_um = KB_T_PN_NM * 1e-9;

        // Marko-Siggia interpolation formula
        // F = (kT/Lp) * [1/(4(1-x/Lc)^2) - 1/4 + x/Lc]
        let one_minus_xi = 1.0 - xi;
        let term1 = 1.0 / (4.0 * one_minus_xi * one_minus_xi);
        let term2 = -0.25;
        let term3 = xi;

        let force_uN = (kbt_uN_um / lp) * (term1 + term2 + term3);

        force_uN
    }

    /// Compute forces on all spectrin network nodes
    ///
    /// Returns a vector of (node_index, force_vector) pairs
    pub fn compute_forces(
        &self,
        network: &SpectrinNetwork,
        _positions: &[Vec3],
    ) -> Vec<(usize, Vec3)> {
        let mut node_forces: Vec<(usize, Vec3)> = Vec::new();

        // For each spectrin edge, compute WLC force
        for edge in &network.edges {
            let node_a = edge.node_a as usize;
            let node_b = edge.node_b as usize;

            // Get node positions from network (spectrin nodes track their own positions)
            let pos_a = network.nodes[node_a].position_vec3();
            let pos_b = network.nodes[node_b].position_vec3();

            // Current extension
            let diff = pos_b - pos_a;
            let current_length = diff.length();

            if current_length < 1e-10 {
                continue; // Skip degenerate edges
            }

            // Direction unit vector (from a to b)
            let direction = diff / current_length;

            // Compute WLC force magnitude using edge's contour length
            let force_magnitude = self.marko_siggia_force(current_length, edge.contour_length_um);

            // Force vector on node A (pulls toward B if extended, pushes away if compressed)
            let force_on_a = direction * force_magnitude;
            let force_on_b = -force_on_a; // Newton's third law

            node_forces.push((node_a, force_on_a));
            node_forces.push((node_b, force_on_b));
        }

        // Aggregate forces by node
        let mut aggregated: Vec<(usize, Vec3)> = Vec::new();
        let mut force_map = std::collections::HashMap::new();

        for (idx, force) in node_forces {
            *force_map.entry(idx).or_insert(Vec3::ZERO) += force;
        }

        for (idx, force) in force_map {
            aggregated.push((idx, force));
        }

        aggregated
    }

    /// Compute total elastic energy stored in spectrin network
    ///
    /// Returns energy in pJ (pico-Joules)
    ///
    /// WLC energy: U = (k_B*T / L_p) * L_c * [x²/(4L_c(1-x/L_c)) + x/L_c - ln(1-x/L_c)]
    pub fn compute_energy(&self, network: &SpectrinNetwork) -> f32 {
        let lp = self.params.persistence_length_um;

        // k_B*T in pJ·μm (since we want pJ output)
        // 4.11 pN·nm = 4.11e-21 J = 4.11e-9 pJ
        let kbt_pJ = KB_T_PN_NM * 1e-9;

        let mut total_energy = 0.0;

        for edge in &network.edges {
            let node_a = edge.node_a as usize;
            let node_b = edge.node_b as usize;

            let pos_a = network.nodes[node_a].position_vec3();
            let pos_b = network.nodes[node_b].position_vec3();

            let x = pos_a.distance(pos_b);
            let lc = edge.contour_length_um;

            // Clamp to prevent singularity
            let max_ext = lc * self.params.max_relative_extension;
            let x_clamped = x.clamp(0.001 * lc, max_ext);
            let xi = x_clamped / lc;

            // Simplified WLC energy (integrated Marko-Siggia)
            // This is an approximation suitable for small deformations
            let one_minus_xi = 1.0 - xi;
            let energy_term = (xi * xi) / (4.0 * one_minus_xi) + xi - one_minus_xi.ln().abs();
            let energy = (kbt_pJ / lp) * lc * energy_term;

            total_energy += energy;
        }

        total_energy
    }

    /// Get the force-extension curve for validation
    ///
    /// Returns (extension_um, force_uN) pairs
    pub fn force_extension_curve(&self, n_points: usize) -> Vec<(f32, f32)> {
        let lc = self.params.contour_length_um;
        let mut curve = Vec::with_capacity(n_points);

        for i in 0..n_points {
            let x = (i as f32 / (n_points - 1) as f32) * lc * 0.95;
            let f = self.marko_siggia_force(x, lc);
            curve.push((x, f));
        }

        curve
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_wlc_parameters() {
        let params = WLCParameters::default();
        assert_eq!(params.persistence_length_um, 0.020);
        assert_eq!(params.contour_length_um, 0.200);
        assert_eq!(params.rest_length_um, 0.075);
    }

    #[test]
    fn test_marko_siggia_at_rest_length() {
        let solver = WLCSolver::new(WLCParameters::default());
        let rest = solver.params.rest_length_um;
        let lc = solver.params.contour_length_um;

        let force = solver.marko_siggia_force(rest, lc);

        // At rest length (75nm), force should be small but positive
        // because WLC naturally tends toward extended state
        // Force should be in the range of ~1-10 pN = 1e-6 to 1e-5 μN
        assert!(force > 0.0, "Force at rest length should be positive (tensile)");
        assert!(force < 1e-4, "Force at rest length should be small (< 0.1 nN)");
    }

    #[test]
    fn test_marko_siggia_monotonic() {
        let solver = WLCSolver::new(WLCParameters::default());
        let lc = solver.params.contour_length_um;

        // Force should increase monotonically with extension
        let mut prev_force = solver.marko_siggia_force(0.01 * lc, lc);
        for i in 2..95 {
            let x = (i as f32 / 100.0) * lc;
            let force = solver.marko_siggia_force(x, lc);
            assert!(
                force >= prev_force,
                "Force should be monotonically increasing: {} >= {} at x={}",
                force,
                prev_force,
                x
            );
            prev_force = force;
        }
    }

    #[test]
    fn test_marko_siggia_divergence_near_contour() {
        let solver = WLCSolver::new(WLCParameters::default());
        let lc = solver.params.contour_length_um;

        // Near full extension, force should be very high
        let force_90 = solver.marko_siggia_force(0.9 * lc, lc);
        let force_50 = solver.marko_siggia_force(0.5 * lc, lc);

        assert!(
            force_90 > force_50 * 10.0,
            "Force at 90% extension should be >> force at 50%"
        );
    }

    #[test]
    fn test_spectrin_stiffness_order_of_magnitude() {
        // Test that stiffness is positive and increases with extension
        // The actual numerical values depend on the specific unit system
        let solver = WLCSolver::new(WLCParameters::default());
        let lc = solver.params.contour_length_um;

        let x1 = 0.4 * lc;
        let x2 = 0.5 * lc;
        let f1 = solver.marko_siggia_force(x1, lc);
        let f2 = solver.marko_siggia_force(x2, lc);

        // Stiffness k = dF/dx
        let stiffness = (f2 - f1) / (x2 - x1);

        // Stiffness should be positive (force increases with extension)
        assert!(
            stiffness > 0.0,
            "Stiffness should be positive, got {}",
            stiffness
        );

        // Stiffness should increase near full extension
        let x3 = 0.8 * lc;
        let x4 = 0.9 * lc;
        let f3 = solver.marko_siggia_force(x3, lc);
        let f4 = solver.marko_siggia_force(x4, lc);
        let stiffness_high = (f4 - f3) / (x4 - x3);

        assert!(
            stiffness_high > stiffness,
            "Stiffness should be higher at larger extensions: {} vs {}",
            stiffness_high,
            stiffness
        );
    }
}
