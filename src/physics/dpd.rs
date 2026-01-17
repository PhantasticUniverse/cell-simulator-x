//! Dissipative Particle Dynamics (DPD) solver.
//!
//! DPD is a mesoscale simulation method that includes:
//! - Conservative forces (soft repulsion)
//! - Dissipative forces (velocity-dependent friction)
//! - Random forces (thermal fluctuations)
//!
//! Force components:
//! F_C = a_ij * (1 - r/r_c) * r̂        (conservative)
//! F_D = -γ * w(r) * (r̂·v_ij) * r̂     (dissipative)
//! F_R = σ * w(r) * θ * r̂ / √dt       (random)
//!
//! Fluctuation-dissipation theorem: σ² = 2*γ*k_B*T
//!
//! Reference: Groot & Warren, J Chem Phys 1997
//! Reference: Español & Warren, Europhys Lett 1995

use glam::Vec3;
use rand::prelude::*;
use rand_distr::StandardNormal;

/// DPD particle for fluid dynamics
#[derive(Debug, Clone)]
pub struct DPDParticle {
    /// Position in μm
    pub position_um: Vec3,
    /// Velocity in μm/s
    pub velocity_um_per_sec: Vec3,
    /// Force in μN
    pub force_uN: Vec3,
    /// Mass in pg (picograms)
    pub mass_pg: f32,
}

impl DPDParticle {
    pub fn new(position: Vec3, mass: f32) -> Self {
        Self {
            position_um: position,
            velocity_um_per_sec: Vec3::ZERO,
            force_uN: Vec3::ZERO,
            mass_pg: mass,
        }
    }
}

/// DPD solver parameters
#[derive(Debug, Clone)]
pub struct DPDParameters {
    /// Conservative force coefficient (μN/μm)
    pub conservative_coeff_a: f32,
    /// Dissipation coefficient γ (μN·s/μm²)
    /// Reference: γ = 4.5 for typical DPD
    pub dissipation_gamma: f32,
    /// Random force coefficient σ
    /// Satisfies σ² = 2*γ*k_B*T
    pub random_sigma: f32,
    /// Cutoff radius in μm
    pub cutoff_radius_um: f32,
    /// Weight function exponent (typically 0.5 or 1.0)
    pub weight_exponent: f32,
}

impl Default for DPDParameters {
    fn default() -> Self {
        // Default values based on typical DPD simulations
        // γ = 4.5, at T = 310K
        let gamma: f32 = 4.5;
        let kbt: f32 = 4.11e-9; // pN·nm in μN·μm units (4.11 pN·nm = 4.11e-9 μN·μm)
        let sigma: f32 = (2.0 * gamma * kbt).sqrt();

        Self {
            conservative_coeff_a: 25.0, // Typical DPD value
            dissipation_gamma: gamma,
            random_sigma: sigma,
            cutoff_radius_um: 1.0,  // 1 μm cutoff
            weight_exponent: 0.5,   // w(r) = (1-r/r_c)^0.5
        }
    }
}

impl DPDParameters {
    /// Create parameters for specific temperature
    pub fn at_temperature(temperature_K: f32) -> Self {
        let mut params = Self::default();

        // k_B*T at given temperature
        // k_B = 1.38e-23 J/K
        let kbt_J = 1.38e-23_f32 * temperature_K;
        // Convert to μN·μm: 1 J = 1e6 μN * 1e6 μm = 1e12 μN·μm
        let kbt = kbt_J * 1e12_f32;

        // Update sigma for fluctuation-dissipation
        params.random_sigma = (2.0_f32 * params.dissipation_gamma * kbt).sqrt();

        params
    }
}

/// DPD solver
pub struct DPDSolver {
    /// Solver parameters
    pub params: DPDParameters,
    /// Random number generator
    rng: StdRng,
}

impl DPDSolver {
    /// Create a new DPD solver
    pub fn new(params: DPDParameters) -> Self {
        Self {
            params,
            rng: StdRng::from_entropy(),
        }
    }

    /// Compute DPD weight function w(r)
    /// w(r) = (1 - r/r_c)^s for r < r_c, 0 otherwise
    fn weight(&self, r: f32) -> f32 {
        let rc = self.params.cutoff_radius_um;
        if r >= rc || r < 1e-10 {
            0.0
        } else {
            let ratio = 1.0 - r / rc;
            ratio.powf(self.params.weight_exponent)
        }
    }

    /// Compute conservative force between two particles
    /// F_C = a * (1 - r/r_c) * r̂
    fn conservative_force(&self, r: f32, r_hat: Vec3) -> Vec3 {
        let rc = self.params.cutoff_radius_um;
        if r >= rc || r < 1e-10 {
            Vec3::ZERO
        } else {
            let a = self.params.conservative_coeff_a;
            r_hat * a * (1.0 - r / rc)
        }
    }

    /// Compute dissipative force between two particles
    /// F_D = -γ * w(r)² * (r̂·v_ij) * r̂
    fn dissipative_force(&self, r: f32, r_hat: Vec3, v_rel: Vec3) -> Vec3 {
        let w = self.weight(r);
        if w < 1e-10 {
            return Vec3::ZERO;
        }

        let gamma = self.params.dissipation_gamma;
        let rdotv = r_hat.dot(v_rel);

        -r_hat * gamma * w * w * rdotv
    }

    /// Compute random force between two particles
    /// F_R = σ * w(r) * θ * r̂ / √dt
    /// θ is a Gaussian random variable with mean 0 and variance 1
    fn random_force(&mut self, r: f32, r_hat: Vec3, dt: f32) -> Vec3 {
        let w = self.weight(r);
        if w < 1e-10 {
            return Vec3::ZERO;
        }

        let sigma = self.params.random_sigma;
        let theta: f32 = self.rng.sample(StandardNormal);

        // Scale by 1/sqrt(dt) for proper fluctuation-dissipation
        let dt_factor = if dt > 0.0 { 1.0 / dt.sqrt() } else { 0.0 };

        r_hat * sigma * w * theta * dt_factor
    }

    /// Compute all DPD forces between particles
    ///
    /// Returns forces per particle
    pub fn compute_forces(
        &mut self,
        particles: &[DPDParticle],
        dt: f32,
    ) -> Vec<Vec3> {
        let n = particles.len();
        let mut forces = vec![Vec3::ZERO; n];
        let rc = self.params.cutoff_radius_um;

        // Pairwise interactions (O(n²) - should use cell lists for large systems)
        for i in 0..n {
            for j in (i + 1)..n {
                let r_vec = particles[j].position_um - particles[i].position_um;
                let r = r_vec.length();

                if r >= rc || r < 1e-10 {
                    continue;
                }

                let r_hat = r_vec / r;
                let v_rel = particles[i].velocity_um_per_sec - particles[j].velocity_um_per_sec;

                // Compute force components
                let f_c = self.conservative_force(r, r_hat);
                let f_d = self.dissipative_force(r, r_hat, v_rel);
                let f_r = self.random_force(r, r_hat, dt);

                let total_force = f_c + f_d + f_r;

                // Newton's third law
                forces[i] += total_force;
                forces[j] -= total_force;
            }
        }

        forces
    }

    /// Compute DPD forces on membrane vertices (simplified version)
    ///
    /// This applies dissipative and random forces to membrane vertices
    /// to model thermal fluctuations and viscous damping.
    pub fn compute_membrane_forces(
        &mut self,
        positions: &[Vec3],
        velocities: &[Vec3],
        temperature_K: f32,
    ) -> Vec<Vec3> {
        let n = positions.len();
        let mut forces = vec![Vec3::ZERO; n];

        // Update sigma for current temperature
        let kbt_J = 1.38e-23 * temperature_K;
        let kbt = kbt_J * 1e12; // Convert to μN·μm
        let sigma = (2.0 * self.params.dissipation_gamma * kbt).sqrt();
        let gamma = self.params.dissipation_gamma;

        // Apply dissipative and random forces to each vertex
        // This is a simplified single-particle version for membrane dynamics
        for i in 0..n {
            let vel = velocities[i];

            // Dissipative force: -γ * v (simple viscous damping)
            let f_d = -vel * gamma * 0.001; // Scale factor for stability

            // Random force: σ * ξ where ξ is Gaussian white noise
            let xi = Vec3::new(
                self.rng.sample::<f32, _>(StandardNormal),
                self.rng.sample::<f32, _>(StandardNormal),
                self.rng.sample::<f32, _>(StandardNormal),
            );
            let f_r = xi * sigma * 0.0001; // Scale factor for stability

            forces[i] = f_d + f_r;
        }

        forces
    }

    /// Perform a full DPD timestep on particles
    pub fn step(&mut self, particles: &mut [DPDParticle], dt: f32) {
        // Compute forces
        let forces = self.compute_forces(particles, dt);

        // Update velocities and positions (simple Euler)
        for (i, particle) in particles.iter_mut().enumerate() {
            let accel = forces[i] / particle.mass_pg;
            particle.velocity_um_per_sec += accel * dt;
            particle.position_um += particle.velocity_um_per_sec * dt;
            particle.force_uN = forces[i];
        }
    }

    /// Get total kinetic energy of DPD particles
    pub fn kinetic_energy(&self, particles: &[DPDParticle]) -> f32 {
        particles
            .iter()
            .map(|p| 0.5 * p.mass_pg * p.velocity_um_per_sec.length_squared())
            .sum()
    }

    /// Get temperature from equipartition theorem
    /// T = (2/3) * KE / (N * k_B)
    pub fn temperature(&self, particles: &[DPDParticle]) -> f32 {
        let ke = self.kinetic_energy(particles);
        let n = particles.len() as f32;
        let kb = 1.38e-23 * 1e12; // k_B in pg·μm²/s²/K (converted units)

        if n > 0.0 && kb > 0.0 {
            (2.0 / 3.0) * ke / (n * kb)
        } else {
            0.0
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dpd_parameters_default() {
        let params = DPDParameters::default();
        assert_eq!(params.dissipation_gamma, 4.5);
        assert!(params.random_sigma > 0.0);
    }

    #[test]
    fn test_weight_function() {
        let solver = DPDSolver::new(DPDParameters::default());

        // At r = 0, w = 1
        let w0 = solver.weight(0.001);
        assert!(w0 > 0.99, "Weight at r≈0 should be ~1, got {}", w0);

        // At r = r_c, w = 0
        let w_rc = solver.weight(solver.params.cutoff_radius_um);
        assert!(w_rc < 0.01, "Weight at r=r_c should be ~0, got {}", w_rc);

        // Weight should decrease with r
        let w_half = solver.weight(solver.params.cutoff_radius_um / 2.0);
        assert!(w_half > 0.0 && w_half < 1.0);
    }

    #[test]
    fn test_conservative_force_repulsive() {
        let solver = DPDSolver::new(DPDParameters::default());

        let r = 0.5; // Within cutoff
        let r_hat = Vec3::X;

        let force = solver.conservative_force(r, r_hat);

        // Conservative force should be repulsive (positive in r_hat direction)
        assert!(force.x > 0.0, "Conservative force should be repulsive");
    }

    #[test]
    fn test_dissipative_force_direction() {
        let solver = DPDSolver::new(DPDParameters::default());

        let r = 0.5;
        let r_hat = Vec3::X;

        // Particles approaching (v_rel points opposite to r_hat)
        let v_rel_approach = Vec3::new(-1.0, 0.0, 0.0);
        let f_approach = solver.dissipative_force(r, r_hat, v_rel_approach);

        // Dissipative force should oppose relative motion
        // When approaching, r̂·v < 0, so F_D should point in +r̂ direction
        assert!(f_approach.x > 0.0, "Dissipative force should resist approach");

        // Particles separating
        let v_rel_separate = Vec3::new(1.0, 0.0, 0.0);
        let f_separate = solver.dissipative_force(r, r_hat, v_rel_separate);

        // When separating, r̂·v > 0, so F_D should point in -r̂ direction
        assert!(f_separate.x < 0.0, "Dissipative force should resist separation");
    }

    #[test]
    fn test_fluctuation_dissipation() {
        // Verify σ² = 2*γ*k_B*T
        let params = DPDParameters::default();

        let gamma = params.dissipation_gamma;
        let sigma = params.random_sigma;

        let kbt = 4.11e-9; // k_B*T at 37°C in μN·μm

        let expected_sigma_sq = 2.0 * gamma * kbt;
        let actual_sigma_sq = sigma * sigma;

        assert!(
            (actual_sigma_sq - expected_sigma_sq).abs() / expected_sigma_sq < 0.01,
            "Fluctuation-dissipation relation not satisfied: {} vs {}",
            actual_sigma_sq,
            expected_sigma_sq
        );
    }

    #[test]
    fn test_membrane_forces_generation() {
        let mut solver = DPDSolver::new(DPDParameters::default());

        let positions = vec![Vec3::ZERO, Vec3::X, Vec3::Y];
        let velocities = vec![Vec3::ZERO; 3];

        let forces = solver.compute_membrane_forces(&positions, &velocities, 310.0);

        assert_eq!(forces.len(), 3);
        // Forces should be non-zero due to random component
        // (could be zero by chance, but very unlikely for all three)
    }
}
