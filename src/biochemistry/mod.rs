//! Biochemistry module for RBC metabolism simulation.
//!
//! This module implements metabolic pathways specific to human red blood cells:
//! - Glycolysis (sole ATP source in mature RBCs)
//! - Rapoport-Luebering shunt (2,3-BPG regulation)
//!
//! RBCs are unique in that they:
//! - Lack mitochondria (no oxidative phosphorylation)
//! - Rely entirely on anaerobic glycolysis for ATP
//! - Contain high levels of 2,3-BPG for hemoglobin regulation
//!
//! Key targets (from PRD):
//! - ATP concentration: 1.5-2.5 mM
//! - 2,3-DPG concentration: 4.5-5.5 mM
//! - Glucose consumption: 1.2-1.5 μmol/hr/mL RBCs
//!
//! References:
//! - Rapoport TA et al. Eur J Biochem. 1976;69:571-584
//! - Joshi A, Palsson BO. J Theor Biol. 1989;141:515-528
//! - Mulquiney PJ, Kuchel PW. Biochem J. 1999;342:567-580

pub mod enzyme;
pub mod glycolysis;
pub mod integrator;
pub mod rapoport_luebering;

pub use enzyme::{Enzyme, ReactionStoichiometry};
pub use glycolysis::{GlycolysisSolver, MetaboliteIndices};
pub use integrator::{IntegratorConfig, RK4Integrator};
pub use rapoport_luebering::{RapoportLueberingSolver, calculate_p50_from_dpg, estimate_dpg_from_ph};

/// Metabolite concentration pool
///
/// Stores all metabolite concentrations as a contiguous vector for efficient
/// ODE integration. Indices are defined by `MetaboliteIndices`.
#[derive(Debug, Clone)]
pub struct MetabolitePool {
    /// Concentrations in mM
    pub concentrations_mM: Vec<f64>,
}

impl MetabolitePool {
    /// Create a new metabolite pool with given size
    pub fn new(n_metabolites: usize) -> Self {
        Self {
            concentrations_mM: vec![0.0; n_metabolites],
        }
    }

    /// Create a pool with default physiological concentrations
    pub fn default_physiological() -> Self {
        let indices = ExtendedMetaboliteIndices::default();
        let mut pool = Self::new(indices.total_count());

        // Glycolytic intermediates (mM)
        // Reference: Rapoport 1976, Joshi-Palsson 1989
        pool.set(indices.glycolysis.glucose, 5.0);
        pool.set(indices.glycolysis.glucose_6_phosphate, 0.05);
        pool.set(indices.glycolysis.fructose_6_phosphate, 0.02);
        pool.set(indices.glycolysis.fructose_1_6_bisphosphate, 0.01);
        pool.set(indices.glycolysis.dihydroxyacetone_phosphate, 0.1);
        pool.set(indices.glycolysis.glyceraldehyde_3_phosphate, 0.005);
        pool.set(indices.glycolysis.bisphosphoglycerate_1_3, 0.001);
        pool.set(indices.glycolysis.phosphoglycerate_3, 0.1);
        pool.set(indices.glycolysis.phosphoglycerate_2, 0.02);
        pool.set(indices.glycolysis.phosphoenolpyruvate, 0.02);
        pool.set(indices.glycolysis.pyruvate, 0.1);
        pool.set(indices.glycolysis.lactate, 1.5);

        // Cofactors (mM)
        // Reference: Minakami & Yoshikawa 1966, Zerez et al. 1987
        pool.set(indices.glycolysis.atp, 2.0);
        pool.set(indices.glycolysis.adp, 0.25);
        pool.set(indices.glycolysis.nad, 0.07);
        pool.set(indices.glycolysis.nadh, 0.03);
        pool.set(indices.glycolysis.pi, 1.0);

        // 2,3-BPG (mM)
        // Reference: Benesch & Benesch 1969
        pool.set(indices.bisphosphoglycerate_2_3, 5.0);

        pool
    }

    /// Get concentration at index
    #[inline]
    pub fn get(&self, idx: usize) -> f64 {
        self.concentrations_mM.get(idx).copied().unwrap_or(0.0)
    }

    /// Set concentration at index
    #[inline]
    pub fn set(&mut self, idx: usize, value: f64) {
        if idx < self.concentrations_mM.len() {
            self.concentrations_mM[idx] = value.max(0.0);
        }
    }

    /// Get mutable reference to concentrations for integration
    pub fn as_mut_slice(&mut self) -> &mut [f64] {
        &mut self.concentrations_mM
    }

    /// Get reference to concentrations
    pub fn as_slice(&self) -> &[f64] {
        &self.concentrations_mM
    }

    /// Total number of metabolites
    pub fn len(&self) -> usize {
        self.concentrations_mM.len()
    }

    /// Check if empty
    pub fn is_empty(&self) -> bool {
        self.concentrations_mM.is_empty()
    }
}

/// Extended metabolite indices including 2,3-BPG
#[derive(Debug, Clone, Copy)]
pub struct ExtendedMetaboliteIndices {
    /// Glycolysis indices
    pub glycolysis: MetaboliteIndices,
    /// 2,3-BPG index (extends glycolysis)
    pub bisphosphoglycerate_2_3: usize,
}

impl Default for ExtendedMetaboliteIndices {
    fn default() -> Self {
        Self {
            glycolysis: MetaboliteIndices::default(),
            bisphosphoglycerate_2_3: 17,  // After glycolysis indices (0-16)
        }
    }
}

impl ExtendedMetaboliteIndices {
    /// Total number of metabolites
    pub fn total_count(&self) -> usize {
        18  // 17 glycolysis + 1 for 2,3-BPG
    }
}

/// Configuration for metabolism simulation
#[derive(Debug, Clone)]
pub struct MetabolismConfig {
    /// Timestep for ODE integration (seconds)
    pub dt_sec: f64,
    /// Temperature in Kelvin
    pub temperature_K: f64,
    /// Intracellular pH
    pub ph: f64,
    /// External glucose concentration (mM)
    pub external_glucose_mM: f64,
    /// Enable glucose transport equilibration
    pub enable_glucose_transport: bool,
}

impl Default for MetabolismConfig {
    fn default() -> Self {
        Self {
            dt_sec: 0.001,        // 1 ms timestep
            temperature_K: 310.0,  // 37°C
            ph: 7.2,
            external_glucose_mM: 5.0,
            enable_glucose_transport: true,
        }
    }
}

/// Main metabolism solver combining all pathways
///
/// Integrates:
/// - Glycolysis (11 reactions)
/// - Rapoport-Luebering shunt (2 reactions)
/// - ATP consumption (external demand)
pub struct MetabolismSolver {
    /// Glycolysis pathway solver
    pub glycolysis: GlycolysisSolver,
    /// Rapoport-Luebering shunt solver
    pub shunt: RapoportLueberingSolver,
    /// ODE integrator
    pub integrator: RK4Integrator,
    /// Metabolite indices
    pub indices: ExtendedMetaboliteIndices,
    /// Configuration
    pub config: MetabolismConfig,
    /// Current simulation time (seconds)
    pub time_sec: f64,
}

impl MetabolismSolver {
    /// Create a new metabolism solver with default parameters
    pub fn new(config: MetabolismConfig) -> Self {
        let indices = ExtendedMetaboliteIndices::default();
        let glycolysis = GlycolysisSolver::new();
        let shunt = RapoportLueberingSolver::new(&indices.glycolysis, indices.bisphosphoglycerate_2_3);

        let integrator = RK4Integrator::new(
            indices.total_count(),
            IntegratorConfig {
                dt_sec: config.dt_sec,
                max_change_mM: 0.1,
                min_concentration_mM: 1e-9,
            },
        );

        Self {
            glycolysis,
            shunt,
            integrator,
            indices,
            config,
            time_sec: 0.0,
        }
    }

    /// Perform one integration step
    ///
    /// # Arguments
    /// * `metabolites` - Metabolite pool to update
    /// * `atp_consumption_mM_per_sec` - External ATP demand (membrane pumps, etc.)
    pub fn step(&mut self, metabolites: &mut MetabolitePool, atp_consumption_mM_per_sec: f64) {
        let indices = self.indices;
        let config = self.config.clone();

        // Closure that computes all derivatives
        let derivatives = |state: &[f64], dydt: &mut [f64]| {
            // Reset derivatives
            for d in dydt.iter_mut() {
                *d = 0.0;
            }

            // Create temporary pool for rate calculations
            let temp_pool = MetabolitePool {
                concentrations_mM: state.to_vec(),
            };

            // Glycolysis reactions
            self.glycolysis.compute_derivatives(&temp_pool, dydt);

            // Rapoport-Luebering shunt
            self.shunt.compute_derivatives(&temp_pool, dydt);

            // External ATP consumption (membrane pumps, etc.)
            dydt[indices.glycolysis.atp] -= atp_consumption_mM_per_sec;
            dydt[indices.glycolysis.adp] += atp_consumption_mM_per_sec;

            // Glucose transport (equilibration with external)
            if config.enable_glucose_transport {
                let glc_internal = state[indices.glycolysis.glucose];
                let glc_external = config.external_glucose_mM;
                // Simple facilitated transport (GLUT1)
                let transport_rate = 0.5 * (glc_external - glc_internal);  // mM/s
                dydt[indices.glycolysis.glucose] += transport_rate;
            }

            // Lactate export (constant rate assumption)
            let lac = state[indices.glycolysis.lactate];
            if lac > 1.0 {
                let export_rate = 0.1 * (lac - 1.0);  // mM/s, above 1mM threshold
                dydt[indices.glycolysis.lactate] -= export_rate;
            }
        };

        // Integrate
        self.integrator.step(metabolites.as_mut_slice(), derivatives);
        self.time_sec = self.integrator.time_sec;
    }

    /// Run simulation for a given duration
    ///
    /// # Arguments
    /// * `metabolites` - Metabolite pool to update
    /// * `duration_sec` - Simulation duration in seconds
    /// * `atp_consumption_mM_per_sec` - External ATP demand
    pub fn run(
        &mut self,
        metabolites: &mut MetabolitePool,
        duration_sec: f64,
        atp_consumption_mM_per_sec: f64,
    ) {
        let n_steps = (duration_sec / self.config.dt_sec).ceil() as usize;
        for _ in 0..n_steps {
            self.step(metabolites, atp_consumption_mM_per_sec);
        }
    }

    /// Get diagnostic information about current state
    pub fn diagnostics(&self, metabolites: &MetabolitePool) -> MetabolismDiagnostics {
        let idx = &self.indices;

        // Calculate key metrics
        let atp = metabolites.get(idx.glycolysis.atp);
        let adp = metabolites.get(idx.glycolysis.adp);
        let dpg = metabolites.get(idx.bisphosphoglycerate_2_3);
        let nad = metabolites.get(idx.glycolysis.nad);
        let nadh = metabolites.get(idx.glycolysis.nadh);
        let glucose = metabolites.get(idx.glycolysis.glucose);
        let lactate = metabolites.get(idx.glycolysis.lactate);

        // ATP/ADP ratio
        let atp_adp_ratio = if adp > 1e-9 { atp / adp } else { f64::INFINITY };

        // NADH/NAD+ ratio
        let nadh_nad_ratio = if nad > 1e-9 { nadh / nad } else { f64::INFINITY };

        // Get reaction rates
        let glycolysis_rates = self.glycolysis.get_rates(metabolites);
        let shunt_rates = self.shunt.get_rates(metabolites);

        // Flux through HK (glucose consumption)
        let hk_rate = glycolysis_rates.iter()
            .find(|(name, _)| *name == "Hexokinase")
            .map(|(_, r)| *r)
            .unwrap_or(0.0);

        // Flux through PGK (main pathway)
        let pgk_rate = glycolysis_rates.iter()
            .find(|(name, _)| *name == "Phosphoglycerate Kinase")
            .map(|(_, r)| *r)
            .unwrap_or(0.0);

        // Shunt ratio
        let shunt_ratio = self.shunt.shunt_ratio(metabolites, pgk_rate);

        // Estimated P50
        let p50_mmHg = calculate_p50_from_dpg(dpg);

        MetabolismDiagnostics {
            time_sec: self.time_sec,
            atp_mM: atp,
            adp_mM: adp,
            dpg_2_3_mM: dpg,
            glucose_mM: glucose,
            lactate_mM: lactate,
            nad_mM: nad,
            nadh_mM: nadh,
            atp_adp_ratio,
            nadh_nad_ratio,
            glucose_consumption_mM_per_sec: hk_rate,
            shunt_ratio,
            p50_mmHg,
            glycolysis_rates,
            shunt_rates,
        }
    }

    /// Reset solver state
    pub fn reset(&mut self) {
        self.integrator.reset();
        self.time_sec = 0.0;
    }

    /// Check if metabolite concentrations are within physiological range
    pub fn validate_state(&self, metabolites: &MetabolitePool) -> Vec<String> {
        let mut warnings = Vec::new();
        let idx = &self.indices;

        // ATP validation (1.5-2.5 mM)
        let atp = metabolites.get(idx.glycolysis.atp);
        if atp < 1.0 {
            warnings.push(format!("ATP critically low: {:.3} mM (target: 1.5-2.5 mM)", atp));
        } else if atp < 1.5 {
            warnings.push(format!("ATP below target: {:.3} mM (target: 1.5-2.5 mM)", atp));
        } else if atp > 3.0 {
            warnings.push(format!("ATP above target: {:.3} mM (target: 1.5-2.5 mM)", atp));
        }

        // 2,3-DPG validation (4.5-5.5 mM)
        let dpg = metabolites.get(idx.bisphosphoglycerate_2_3);
        if dpg < 3.0 {
            warnings.push(format!("2,3-DPG low: {:.3} mM (target: 4.5-5.5 mM)", dpg));
        } else if dpg > 7.0 {
            warnings.push(format!("2,3-DPG high: {:.3} mM (target: 4.5-5.5 mM)", dpg));
        }

        // NAD/NADH ratio validation (typically 2-7)
        let nad = metabolites.get(idx.glycolysis.nad);
        let nadh = metabolites.get(idx.glycolysis.nadh);
        if nadh > 1e-9 {
            let ratio = nad / nadh;
            if ratio < 1.0 {
                warnings.push(format!("NAD/NADH ratio too low: {:.2} (should be >2)", ratio));
            }
        }

        warnings
    }
}

impl Default for MetabolismSolver {
    fn default() -> Self {
        Self::new(MetabolismConfig::default())
    }
}

/// Diagnostic information from metabolism solver
#[derive(Debug, Clone)]
pub struct MetabolismDiagnostics {
    /// Current simulation time (seconds)
    pub time_sec: f64,
    /// ATP concentration (mM)
    pub atp_mM: f64,
    /// ADP concentration (mM)
    pub adp_mM: f64,
    /// 2,3-DPG concentration (mM)
    pub dpg_2_3_mM: f64,
    /// Glucose concentration (mM)
    pub glucose_mM: f64,
    /// Lactate concentration (mM)
    pub lactate_mM: f64,
    /// NAD+ concentration (mM)
    pub nad_mM: f64,
    /// NADH concentration (mM)
    pub nadh_mM: f64,
    /// ATP/ADP ratio
    pub atp_adp_ratio: f64,
    /// NADH/NAD+ ratio
    pub nadh_nad_ratio: f64,
    /// Glucose consumption rate (mM/s)
    pub glucose_consumption_mM_per_sec: f64,
    /// Fraction of 1,3-BPG going through shunt
    pub shunt_ratio: f64,
    /// Estimated P50 (mmHg) based on 2,3-DPG
    pub p50_mmHg: f64,
    /// Individual glycolysis reaction rates
    pub glycolysis_rates: Vec<(&'static str, f64)>,
    /// Shunt reaction rates
    pub shunt_rates: Vec<(&'static str, f64)>,
}

impl MetabolismDiagnostics {
    /// Print a formatted summary
    pub fn print_summary(&self) {
        println!("=== Metabolism State (t = {:.3} s) ===", self.time_sec);
        println!();
        println!("Key Metabolites:");
        println!("  ATP:      {:.3} mM  (target: 1.5-2.5 mM)", self.atp_mM);
        println!("  ADP:      {:.3} mM", self.adp_mM);
        println!("  2,3-DPG:  {:.3} mM  (target: 4.5-5.5 mM)", self.dpg_2_3_mM);
        println!("  Glucose:  {:.3} mM", self.glucose_mM);
        println!("  Lactate:  {:.3} mM", self.lactate_mM);
        println!();
        println!("Ratios:");
        println!("  ATP/ADP:   {:.2}", self.atp_adp_ratio);
        println!("  NADH/NAD+: {:.3}", self.nadh_nad_ratio);
        println!();
        println!("Fluxes:");
        println!("  Glucose consumption: {:.4} mM/s ({:.2} μmol/hr/mL)",
            self.glucose_consumption_mM_per_sec,
            self.glucose_consumption_mM_per_sec * 3600.0 * 1000.0);  // Convert to μmol/hr/mL
        println!("  Shunt ratio: {:.1}%", self.shunt_ratio * 100.0);
        println!();
        println!("O2 Affinity:");
        println!("  Estimated P50: {:.1} mmHg (target: ~27 mmHg)", self.p50_mmHg);
    }

    /// Print detailed reaction rates
    pub fn print_rates(&self) {
        println!("\nGlycolysis Reaction Rates (mM/s):");
        for (name, rate) in &self.glycolysis_rates {
            println!("  {}: {:.6}", name, rate);
        }
        println!("\nShunt Reaction Rates (mM/s):");
        for (name, rate) in &self.shunt_rates {
            println!("  {}: {:.6}", name, rate);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_metabolite_pool_default() {
        let pool = MetabolitePool::default_physiological();

        let indices = ExtendedMetaboliteIndices::default();
        assert!((pool.get(indices.glycolysis.atp) - 2.0).abs() < 0.01);
        assert!((pool.get(indices.glycolysis.glucose) - 5.0).abs() < 0.01);
        assert!((pool.get(indices.bisphosphoglycerate_2_3) - 5.0).abs() < 0.01);
    }

    #[test]
    fn test_metabolism_solver_creation() {
        let solver = MetabolismSolver::new(MetabolismConfig::default());
        assert_eq!(solver.time_sec, 0.0);
    }

    #[test]
    fn test_metabolism_step() {
        let mut solver = MetabolismSolver::new(MetabolismConfig::default());
        let mut metabolites = MetabolitePool::default_physiological();

        let _initial_atp = metabolites.get(solver.indices.glycolysis.atp);

        // Run a few steps
        for _ in 0..10 {
            solver.step(&mut metabolites, 0.001);  // Small ATP consumption
        }

        // Time should advance
        assert!(solver.time_sec > 0.0);

        // ATP should still be positive
        let final_atp = metabolites.get(solver.indices.glycolysis.atp);
        assert!(final_atp > 0.0);
    }

    #[test]
    fn test_diagnostics() {
        let solver = MetabolismSolver::new(MetabolismConfig::default());
        let metabolites = MetabolitePool::default_physiological();

        let diag = solver.diagnostics(&metabolites);

        assert!(diag.atp_mM > 0.0);
        assert!(diag.dpg_2_3_mM > 0.0);
        assert!(diag.p50_mmHg > 20.0 && diag.p50_mmHg < 35.0);
    }

    #[test]
    fn test_steady_state_approach() {
        let mut solver = MetabolismSolver::new(MetabolismConfig::default());
        let mut metabolites = MetabolitePool::default_physiological();

        // Run for 10 seconds of simulation time
        solver.run(&mut metabolites, 10.0, 0.001);

        let diag = solver.diagnostics(&metabolites);

        // Should approach steady state with ATP in physiological range
        // Note: exact values depend on kinetic parameters
        assert!(
            diag.atp_mM > 0.5 && diag.atp_mM < 4.0,
            "ATP out of range after steady state: {} mM",
            diag.atp_mM
        );
    }

    #[test]
    fn test_validation() {
        let solver = MetabolismSolver::new(MetabolismConfig::default());
        let metabolites = MetabolitePool::default_physiological();

        let warnings = solver.validate_state(&metabolites);
        // Default physiological state should have no warnings
        assert!(
            warnings.is_empty(),
            "Unexpected warnings for physiological state: {:?}",
            warnings
        );
    }
}
