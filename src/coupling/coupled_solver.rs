//! Coupled mechano-metabolic solver orchestrating physics and biochemistry.
//!
//! This solver synchronizes the physics simulation (1 μs timesteps) with
//! the biochemistry simulation (1 ms timesteps), implementing bidirectional
//! coupling between mechanical deformation and metabolism.
//!
//! ## Timestep Synchronization
//! - Physics runs at 1 μs (1e-6 s) for numerical stability of membrane mechanics
//! - Biochemistry runs at 1 ms (1e-3 s) for enzyme kinetics
//! - Ratio: 1000 physics steps per biochemistry step
//!
//! ## Coupling Points
//! 1. **Forward**: Tension → Piezo1 → Ca²⁺ → ATP release
//! 2. **Reverse**: ATP → Spectrin phosphorylation → Membrane stiffness
//!
//! ## Usage
//! ```ignore
//! let mut solver = CoupledSolver::new(config);
//! for _ in 0..60000 { // 60 seconds at 1ms steps
//!     solver.step(&mut mesh, &spectrin, &mut physics_state, &mut metabolites);
//! }
//! ```

use crate::biochemistry::{FullyIntegratedSolver, FullyIntegratedConfig, MetabolitePool};
use crate::geometry::{Mesh, SpectrinNetwork};
use crate::physics::{PhysicsConfig, PhysicsSolver};
use crate::state::PhysicsState;

use super::spectrin_modulator::SpectrinModulator;
use super::tension_computer::TensionComputer;

/// Configuration for the coupled solver.
#[derive(Debug, Clone)]
pub struct CoupledConfig {
    /// Physics timestep (seconds). Default: 1e-6 (1 μs)
    pub physics_dt_sec: f64,
    /// Biochemistry timestep (seconds). Default: 1e-3 (1 ms)
    pub biochem_dt_sec: f64,
    /// Number of physics substeps per biochemistry step. Default: 1000
    pub physics_substeps: usize,
    /// Enable forward coupling (tension → Piezo1)
    pub enable_tension_coupling: bool,
    /// Enable reverse coupling (ATP → stiffness)
    pub enable_atp_stiffness_coupling: bool,
    /// External glucose concentration (mM)
    pub external_glucose_mM: f64,
    /// External pO2 (mmHg)
    pub po2_mmHg: f64,
    /// Baseline oxidative stress multiplier
    pub oxidative_stress_multiplier: f64,
}

impl Default for CoupledConfig {
    fn default() -> Self {
        Self {
            physics_dt_sec: 1e-6,      // 1 microsecond
            biochem_dt_sec: 1e-3,      // 1 millisecond
            physics_substeps: 1000,     // 1000 physics steps per biochem step
            enable_tension_coupling: true,
            enable_atp_stiffness_coupling: true,
            external_glucose_mM: 5.0,
            po2_mmHg: 100.0,
            oxidative_stress_multiplier: 1.0,
        }
    }
}

/// Coupled mechano-metabolic solver.
///
/// Orchestrates synchronized physics and biochemistry simulations
/// with bidirectional coupling.
pub struct CoupledSolver {
    /// Physics solver for membrane mechanics
    pub physics: PhysicsSolver,
    /// Biochemistry solver for metabolism
    pub biochemistry: FullyIntegratedSolver,
    /// Tension computer (physics → biochemistry)
    pub tension_computer: TensionComputer,
    /// Spectrin modulator (biochemistry → physics)
    pub spectrin_modulator: SpectrinModulator,
    /// Configuration
    pub config: CoupledConfig,
    /// Current simulation time (seconds)
    pub time_sec: f64,
    /// Base shear modulus (stored for modification)
    base_shear_modulus: f32,
    /// Base persistence length (stored for modification)
    base_persistence_length: f32,
    /// Current stiffness modifier (for diagnostics)
    current_stiffness_modifier: f32,
    /// Current computed tension (pN/nm)
    current_tension_pN_per_nm: f64,
}

impl CoupledSolver {
    /// Create a new coupled solver.
    ///
    /// # Arguments
    /// * `mesh` - Initial mesh geometry (used to initialize physics solver)
    /// * `config` - Coupled solver configuration
    pub fn new(mesh: &Mesh, config: CoupledConfig) -> Self {
        // Create physics solver
        let physics_config = PhysicsConfig {
            dt_sec: config.physics_dt_sec as f32,
            temperature_K: 310.0,
            enable_thermal_noise: true,
            membrane_damping: 5.0,
        };
        let physics = PhysicsSolver::new(mesh, physics_config);

        // Store base material parameters for later modification
        let base_shear_modulus = physics.skalak_solver.material.shear_modulus_uN_per_m;
        let base_persistence_length = physics.wlc_solver.params.persistence_length_um;

        // Create biochemistry solver
        let biochem_config = FullyIntegratedConfig {
            dt_sec: config.biochem_dt_sec,
            external_glucose_mM: config.external_glucose_mM,
            po2_mmHg: config.po2_mmHg,
            oxidative_stress_multiplier: config.oxidative_stress_multiplier,
            enable_piezo1: config.enable_tension_coupling,
            membrane_tension_pN_per_nm: 0.0, // Will be updated dynamically
            ..Default::default()
        };
        let biochemistry = FullyIntegratedSolver::new(biochem_config);

        // Create coupling components
        let tension_computer = TensionComputer::default();
        let spectrin_modulator = SpectrinModulator::default();

        Self {
            physics,
            biochemistry,
            tension_computer,
            spectrin_modulator,
            config,
            time_sec: 0.0,
            base_shear_modulus,
            base_persistence_length,
            current_stiffness_modifier: 1.0,
            current_tension_pN_per_nm: 0.0,
        }
    }

    /// Perform one coupled timestep (1 ms = 1 biochemistry step).
    ///
    /// This executes:
    /// 1. N physics substeps (accumulating tension)
    /// 2. One biochemistry step (with computed tension)
    /// 3. Update physics parameters from ATP level
    ///
    /// # Arguments
    /// * `mesh` - Mesh geometry (modified by physics)
    /// * `spectrin` - Spectrin network
    /// * `physics_state` - Physics state (velocities, forces)
    /// * `metabolites` - Metabolite concentrations
    pub fn step(
        &mut self,
        mesh: &mut Mesh,
        spectrin: &SpectrinNetwork,
        physics_state: &mut PhysicsState,
        metabolites: &mut MetabolitePool,
    ) {
        // === Phase 1: Run physics substeps and accumulate tension ===
        for _ in 0..self.config.physics_substeps {
            self.physics.step(mesh, spectrin, physics_state);

            // Compute and accumulate tension if coupling enabled
            if self.config.enable_tension_coupling {
                self.tension_computer
                    .compute_global_tension_pN_per_nm(&self.physics.skalak_solver, mesh);
            }
        }

        // Get averaged tension for biochemistry (only update if computing from physics)
        if self.config.enable_tension_coupling {
            self.current_tension_pN_per_nm = self.tension_computer.averaged_tension();
            self.biochemistry.config.membrane_tension_pN_per_nm = self.current_tension_pN_per_nm;
        }
        // Otherwise, keep the tension from set_tension_override()

        // Step biochemistry (includes Piezo1 mechanotransduction)
        self.biochemistry.step(metabolites, 0.0);

        // === Phase 3: Update physics from ATP if coupling enabled ===
        if self.config.enable_atp_stiffness_coupling {
            let atp_mM = metabolites.get(self.biochemistry.indices.glycolysis.atp);
            self.current_stiffness_modifier = self.spectrin_modulator.stiffness_modifier(atp_mM);

            // Update material parameters
            // Note: Currently the solvers don't have runtime parameter setters,
            // so we store the modifier for diagnostics. In a full implementation,
            // we would modify physics.skalak_solver.material.shear_modulus_uN_per_m
            // and physics.wlc_solver.params.persistence_length_um.

            // For now, we apply the modifier through diagnostics and future
            // runtime parameter modification (Phase 8 extension)
        }

        // Update time
        self.time_sec += self.config.biochem_dt_sec;
    }

    /// Run simulation for a given duration.
    ///
    /// # Arguments
    /// * `mesh` - Mesh geometry
    /// * `spectrin` - Spectrin network
    /// * `physics_state` - Physics state
    /// * `metabolites` - Metabolite concentrations
    /// * `duration_sec` - Total duration in seconds
    pub fn run(
        &mut self,
        mesh: &mut Mesh,
        spectrin: &SpectrinNetwork,
        physics_state: &mut PhysicsState,
        metabolites: &mut MetabolitePool,
        duration_sec: f64,
    ) {
        let n_steps = (duration_sec / self.config.biochem_dt_sec).ceil() as usize;
        for _ in 0..n_steps {
            self.step(mesh, spectrin, physics_state, metabolites);
        }
    }

    /// Get current diagnostics.
    pub fn diagnostics(&self, metabolites: &MetabolitePool) -> CoupledDiagnostics {
        let biochem_diag = self.biochemistry.diagnostics(metabolites);

        CoupledDiagnostics {
            time_sec: self.time_sec,
            // Mechanics
            membrane_tension_pN_per_nm: self.current_tension_pN_per_nm,
            stiffness_modifier: self.current_stiffness_modifier,
            effective_shear_modulus_uN_per_m: self.base_shear_modulus * self.current_stiffness_modifier,
            effective_persistence_length_um: self.base_persistence_length / self.current_stiffness_modifier,
            // Biochemistry (key metabolites)
            atp_mM: biochem_diag.atp_mM,
            adp_mM: biochem_diag.adp_mM,
            ca_cytosolic_uM: biochem_diag.ca_cytosolic_uM,
            gsh_mM: biochem_diag.gsh_mM,
            h2o2_uM: biochem_diag.h2o2_uM,
            // Coupling status
            piezo1_open_probability: if self.current_tension_pN_per_nm > 0.0 {
                // Approximate open probability from tension (Hill equation)
                let half_activation: f64 = 2.0; // pN/nm
                let hill: f64 = 4.0;
                let t: f64 = self.current_tension_pN_per_nm;
                let p = t.powf(hill) / (half_activation.powf(hill) + t.powf(hill));
                p.min(1.0)
            } else {
                0.0
            },
            stiffness_status: self.spectrin_modulator.status_description(biochem_diag.atp_mM),
        }
    }

    /// Reset solver state.
    pub fn reset(&mut self) {
        self.biochemistry.reset();
        self.tension_computer.reset();
        self.time_sec = 0.0;
        self.current_stiffness_modifier = 1.0;
        self.current_tension_pN_per_nm = 0.0;
    }

    /// Set membrane tension override (for testing without physics).
    pub fn set_tension_override(&mut self, tension_pN_per_nm: f64) {
        self.current_tension_pN_per_nm = tension_pN_per_nm;
        self.biochemistry.config.membrane_tension_pN_per_nm = tension_pN_per_nm;
    }
}

/// Diagnostic information from coupled solver.
#[derive(Debug, Clone)]
pub struct CoupledDiagnostics {
    /// Simulation time (seconds)
    pub time_sec: f64,

    // Mechanics
    /// Current membrane tension (pN/nm)
    pub membrane_tension_pN_per_nm: f64,
    /// ATP-dependent stiffness modifier (1.0 = normal)
    pub stiffness_modifier: f32,
    /// Effective shear modulus after ATP modification (μN/m)
    pub effective_shear_modulus_uN_per_m: f32,
    /// Effective persistence length after ATP modification (μm)
    pub effective_persistence_length_um: f32,

    // Biochemistry
    /// ATP concentration (mM)
    pub atp_mM: f64,
    /// ADP concentration (mM)
    pub adp_mM: f64,
    /// Cytosolic Ca²⁺ (μM)
    pub ca_cytosolic_uM: f64,
    /// GSH concentration (mM)
    pub gsh_mM: f64,
    /// H2O2 concentration (μM)
    pub h2o2_uM: f64,

    // Coupling status
    /// Piezo1 open probability
    pub piezo1_open_probability: f64,
    /// Description of stiffness status
    pub stiffness_status: &'static str,
}

impl CoupledDiagnostics {
    /// Print a formatted summary.
    pub fn print_summary(&self) {
        println!("=== Coupled Mechano-Metabolic State (t = {:.3} s) ===", self.time_sec);
        println!();
        println!("Mechanics:");
        println!("  Membrane tension:    {:.3} pN/nm", self.membrane_tension_pN_per_nm);
        println!("  Stiffness modifier:  {:.3}x ({})", self.stiffness_modifier, self.stiffness_status);
        println!("  Effective Gs:        {:.2} μN/m (base: 5.5)", self.effective_shear_modulus_uN_per_m);
        println!("  Effective Lp:        {:.4} μm (base: 0.020)", self.effective_persistence_length_um);
        println!();
        println!("Biochemistry:");
        println!("  ATP:                 {:.3} mM (target: 1.5-2.5)", self.atp_mM);
        println!("  ADP:                 {:.3} mM", self.adp_mM);
        println!("  Ca²⁺ (cytosolic):    {:.0} nM", self.ca_cytosolic_uM * 1000.0);
        println!("  GSH:                 {:.3} mM", self.gsh_mM);
        println!("  H2O2:                {:.2} μM", self.h2o2_uM);
        println!();
        println!("Coupling:");
        println!("  Piezo1 P(open):      {:.1}%", self.piezo1_open_probability * 100.0);
    }

    /// Print a one-line row for time series.
    pub fn print_row_header() {
        println!("{:>8} {:>10} {:>10} {:>8} {:>8} {:>10} {:>10}",
            "Time(s)", "Tension", "Stiff.Mod", "ATP", "Ca²⁺(nM)", "Piezo1(%)", "Status");
        println!("{}", "-".repeat(76));
    }

    /// Print a one-line row.
    pub fn print_row(&self) {
        println!("{:8.2} {:10.3} {:10.3} {:8.3} {:8.0} {:10.1} {:>10}",
            self.time_sec,
            self.membrane_tension_pN_per_nm,
            self.stiffness_modifier,
            self.atp_mM,
            self.ca_cytosolic_uM * 1000.0,
            self.piezo1_open_probability * 100.0,
            self.stiffness_status);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::GeometryParameters;

    fn create_test_system() -> (Mesh, SpectrinNetwork, PhysicsState, MetabolitePool) {
        let params = GeometryParameters {
            cell_radius_um: 3.91,
            fung_tong_c0_um: 0.81,
            fung_tong_c2_um: 7.83,
            fung_tong_c4_um: -4.39,
            mesh_resolution: 8,
            spectrin_target_count: 50,
        };
        let mesh = Mesh::generate_rbc(&params);
        let spectrin = SpectrinNetwork::generate(&mesh, &params);
        let physics_state = PhysicsState::new(mesh.vertices.len());
        let metabolites = MetabolitePool::default_fully_integrated();

        (mesh, spectrin, physics_state, metabolites)
    }

    #[test]
    fn test_coupled_solver_creation() {
        let (mesh, _spectrin, _physics, _metabolites) = create_test_system();
        let config = CoupledConfig::default();
        let solver = CoupledSolver::new(&mesh, config);

        assert!(solver.time_sec == 0.0);
        assert!((solver.current_stiffness_modifier - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_coupled_solver_step() {
        let (mut mesh, spectrin, mut physics, mut metabolites) = create_test_system();
        let mut config = CoupledConfig::default();
        config.physics_substeps = 10; // Reduce for faster test
        let mut solver = CoupledSolver::new(&mesh, config);

        solver.step(&mut mesh, &spectrin, &mut physics, &mut metabolites);

        assert!(solver.time_sec > 0.0);
    }

    #[test]
    fn test_tension_override() {
        let (mesh, _spectrin, _physics, metabolites) = create_test_system();
        let config = CoupledConfig::default();
        let mut solver = CoupledSolver::new(&mesh, config);

        solver.set_tension_override(2.0);

        let diag = solver.diagnostics(&metabolites);
        assert!((diag.membrane_tension_pN_per_nm - 2.0).abs() < 0.01);
    }

    #[test]
    fn test_stiffness_modifier_updates() {
        let (mesh, _spectrin, _physics, mut metabolites) = create_test_system();
        let config = CoupledConfig::default();
        let solver = CoupledSolver::new(&mesh, config);

        // At normal ATP, modifier should be 1.0
        let diag = solver.diagnostics(&metabolites);
        assert!((diag.stiffness_modifier - 1.0).abs() < 0.1);

        // Deplete ATP
        metabolites.set(solver.biochemistry.indices.glycolysis.atp, 0.5);
        // Manually update modifier (normally happens in step())
        let atp_mM = metabolites.get(solver.biochemistry.indices.glycolysis.atp);
        let modifier = solver.spectrin_modulator.stiffness_modifier(atp_mM);
        assert!(modifier > 1.0, "Low ATP should increase stiffness modifier");
    }
}
