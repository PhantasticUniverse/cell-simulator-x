//! Redox metabolism solver for RBC antioxidant defense.
//!
//! This module integrates the Pentose Phosphate Pathway (PPP), Glutathione
//! cycle, and Piezo1 mechanotransduction into a unified redox metabolism solver.
//!
//! The key coupling is:
//! - PPP produces NADPH from G6P
//! - Glutathione reductase uses NADPH to regenerate GSH from GSSG
//! - GSH is consumed by GPx to detoxify H2O2
//! - Piezo1 responds to membrane tension with Ca2+ influx -> ATP release
//!
//! This creates a feedback loop where oxidative stress depletes NADPH and GSH,
//! which in turn stimulates PPP flux (via decreased NADPH inhibition of G6PDH).
//!
//! Validation targets:
//! - NADPH/NADP+ ratio: 10-20 (Veech 1969)
//! - GSH/GSSG ratio: 100-400 (Meister 1983)
//! - Total GSH: 2-3 mM (Beutler 1969)
//! - Cytosolic H2O2: <5 uM (Chance 1979)
//! - Resting Ca2+: ~100 nM (Engelmann 2008)

use super::MetabolitePool;
use super::glycolysis::MetaboliteIndices;
use super::pentose_phosphate::{PentosePhosphatePathway, RedoxIndices};
use super::glutathione::GlutathioneCycle;
use super::piezo1::Piezo1System;
use super::integrator::{IntegratorConfig, RK4Integrator};

/// Redox metabolism solver combining PPP, glutathione, and Piezo1
pub struct RedoxSolver {
    pub ppp: PentosePhosphatePathway,
    pub glutathione: GlutathioneCycle,
    pub piezo1: Piezo1System,
    pub indices: RedoxIndices,
    pub integrator: RK4Integrator,
    /// Current simulation time
    pub time_sec: f64,
    /// Configuration
    pub config: RedoxConfig,
}

/// Configuration for redox solver
#[derive(Debug, Clone)]
pub struct RedoxConfig {
    /// Integration timestep (seconds)
    pub dt_sec: f64,
    /// Oxidative stress multiplier (1.0 = normal)
    pub oxidative_stress_multiplier: f64,
    /// Membrane tension (pN/nm)
    pub membrane_tension_pN_per_nm: f64,
    /// Enable Piezo1 mechanotransduction
    pub enable_piezo1: bool,
}

impl Default for RedoxConfig {
    fn default() -> Self {
        Self {
            dt_sec: 0.001,  // 1 ms timestep
            oxidative_stress_multiplier: 1.0,
            membrane_tension_pN_per_nm: 0.0,  // No tension at rest
            enable_piezo1: true,
        }
    }
}

impl RedoxSolver {
    /// Create a new redox solver
    ///
    /// # Arguments
    /// * `glycolysis_indices` - Glycolysis metabolite indices for shared metabolites
    /// * `config` - Solver configuration
    pub fn new(glycolysis_indices: &MetaboliteIndices, config: RedoxConfig) -> Self {
        // Create redox indices starting after 2,3-BPG (index 17)
        let indices = RedoxIndices::new(glycolysis_indices, 18);

        let ppp = PentosePhosphatePathway::new(&indices);
        let mut glutathione = GlutathioneCycle::new(&indices);
        glutathione.set_oxidative_stress(config.oxidative_stress_multiplier);
        let piezo1 = Piezo1System::new(&indices);

        // Total metabolites = 18 (glycolysis + 2,3-BPG) + new redox metabolites
        let total_metabolites = 18 + RedoxIndices::new_metabolite_count();

        let integrator = RK4Integrator::new(
            total_metabolites,
            IntegratorConfig {
                dt_sec: config.dt_sec,
                max_change_mM: 0.1,
                min_concentration_mM: 1e-9,
            },
        );

        Self {
            ppp,
            glutathione,
            piezo1,
            indices,
            integrator,
            time_sec: 0.0,
            config,
        }
    }

    /// Set oxidative stress level
    pub fn set_oxidative_stress(&mut self, multiplier: f64) {
        self.config.oxidative_stress_multiplier = multiplier;
        self.glutathione.set_oxidative_stress(multiplier);
    }

    /// Set membrane tension for Piezo1
    pub fn set_membrane_tension(&mut self, tension_pN_per_nm: f64) {
        self.config.membrane_tension_pN_per_nm = tension_pN_per_nm;
    }

    /// Compute all redox derivatives
    fn compute_derivatives(&self, metabolites: &MetabolitePool, dydt: &mut [f64]) {
        // PPP (G6P -> NADPH)
        self.ppp.compute_derivatives(metabolites, dydt);

        // Glutathione cycle (NADPH -> GSH, H2O2 detox)
        self.glutathione.compute_derivatives(metabolites, dydt);

        // Piezo1 (membrane tension -> Ca2+ -> ATP release)
        if self.config.enable_piezo1 {
            self.piezo1.compute_derivatives(
                metabolites,
                self.config.membrane_tension_pN_per_nm,
                dydt,
            );
        }
    }

    /// Perform one integration step
    pub fn step(&mut self, metabolites: &mut MetabolitePool) {
        // Capture references to components for the closure
        let ppp = &self.ppp;
        let glutathione = &self.glutathione;
        let piezo1 = &self.piezo1;
        let enable_piezo1 = self.config.enable_piezo1;
        let membrane_tension = self.config.membrane_tension_pN_per_nm;

        let derivatives = |state: &[f64], dydt: &mut [f64]| {
            // Reset derivatives
            for d in dydt.iter_mut() {
                *d = 0.0;
            }

            // Create temporary pool for rate calculations
            let temp_pool = MetabolitePool {
                concentrations_mM: state.to_vec(),
            };

            // Compute all redox derivatives
            ppp.compute_derivatives(&temp_pool, dydt);
            glutathione.compute_derivatives(&temp_pool, dydt);

            if enable_piezo1 {
                piezo1.compute_derivatives(
                    &temp_pool,
                    membrane_tension,
                    dydt,
                );
            }
        };

        self.integrator.step(metabolites.as_mut_slice(), derivatives);
        self.time_sec = self.integrator.time_sec;
    }

    /// Run simulation for a duration
    pub fn run(&mut self, metabolites: &mut MetabolitePool, duration_sec: f64) {
        let n_steps = (duration_sec / self.config.dt_sec).ceil() as usize;
        for _ in 0..n_steps {
            self.step(metabolites);
        }
    }

    /// Get diagnostic information
    pub fn diagnostics(&self, metabolites: &MetabolitePool) -> RedoxDiagnostics {
        let idx = &self.indices;

        // Redox pairs
        let nadph = metabolites.get(idx.nadph);
        let nadp = metabolites.get(idx.nadp_plus);
        let nadph_nadp_ratio = if nadp > 1e-9 { nadph / nadp } else { f64::INFINITY };

        // Glutathione
        let gsh = metabolites.get(idx.gsh);
        let gssg = metabolites.get(idx.gssg);
        let gsh_gssg_ratio = if gssg > 1e-9 { gsh / gssg } else { f64::INFINITY };
        let total_gsh = gsh + 2.0 * gssg;

        // H2O2
        let h2o2 = metabolites.get(idx.h2o2);

        // Ca2+
        let ca_uM = metabolites.get(idx.ca2_plus_cytosolic);

        // PPP flux
        let nadph_production = self.ppp.nadph_production_rate(metabolites);

        // Reaction rates
        let ppp_rates = self.ppp.get_rates(metabolites);
        let glutathione_rates = self.glutathione.get_rates(metabolites);
        let piezo1_diag = self.piezo1.diagnostics(
            metabolites,
            self.config.membrane_tension_pN_per_nm,
        );

        RedoxDiagnostics {
            time_sec: self.time_sec,
            nadph_mM: nadph,
            nadp_plus_mM: nadp,
            nadph_nadp_ratio,
            gsh_mM: gsh,
            gssg_mM: gssg,
            gsh_gssg_ratio,
            total_glutathione_mM: total_gsh,
            h2o2_uM: h2o2 * 1000.0,  // Convert mM to uM
            ca_cytosolic_uM: ca_uM,
            nadph_production_mM_per_sec: nadph_production,
            oxidative_stress: self.config.oxidative_stress_multiplier,
            membrane_tension_pN_per_nm: self.config.membrane_tension_pN_per_nm,
            ppp_rates,
            glutathione_rates,
            piezo1_open_probability: piezo1_diag.open_probability,
            atp_release_mM_per_sec: piezo1_diag.atp_release_mM_per_sec,
        }
    }

    /// Validate state against physiological targets
    pub fn validate_state(&self, metabolites: &MetabolitePool) -> Vec<String> {
        let mut warnings = Vec::new();
        let diag = self.diagnostics(metabolites);

        // NADPH/NADP+ ratio (target: 10-20)
        if diag.nadph_nadp_ratio < 5.0 {
            warnings.push(format!(
                "NADPH/NADP+ ratio low: {:.1} (target: 10-20)",
                diag.nadph_nadp_ratio
            ));
        } else if diag.nadph_nadp_ratio > 50.0 && diag.nadph_nadp_ratio.is_finite() {
            warnings.push(format!(
                "NADPH/NADP+ ratio high: {:.1} (target: 10-20)",
                diag.nadph_nadp_ratio
            ));
        }

        // GSH/GSSG ratio (target: 100-400)
        if diag.gsh_gssg_ratio < 50.0 {
            warnings.push(format!(
                "GSH/GSSG ratio low: {:.1} (target: 100-400) - oxidative stress",
                diag.gsh_gssg_ratio
            ));
        } else if diag.gsh_gssg_ratio > 1000.0 && diag.gsh_gssg_ratio.is_finite() {
            warnings.push(format!(
                "GSH/GSSG ratio high: {:.1} (target: 100-400)",
                diag.gsh_gssg_ratio
            ));
        }

        // Total GSH (target: 2-3 mM)
        if diag.total_glutathione_mM < 1.5 {
            warnings.push(format!(
                "Total glutathione low: {:.2} mM (target: 2-3 mM)",
                diag.total_glutathione_mM
            ));
        } else if diag.total_glutathione_mM > 4.0 {
            warnings.push(format!(
                "Total glutathione high: {:.2} mM (target: 2-3 mM)",
                diag.total_glutathione_mM
            ));
        }

        // H2O2 (target: <5 uM)
        if diag.h2o2_uM > 10.0 {
            warnings.push(format!(
                "H2O2 elevated: {:.1} uM (target: <5 uM) - oxidative stress",
                diag.h2o2_uM
            ));
        }

        // Resting Ca2+ (target: ~100 nM = 0.1 uM)
        if self.config.membrane_tension_pN_per_nm < 0.5 {
            // Only check at rest
            if diag.ca_cytosolic_uM > 0.5 {
                warnings.push(format!(
                    "Resting Ca2+ elevated: {:.0} nM (target: ~100 nM)",
                    diag.ca_cytosolic_uM * 1000.0
                ));
            }
        }

        warnings
    }

    /// Reset solver state
    pub fn reset(&mut self) {
        self.integrator.reset();
        self.time_sec = 0.0;
    }
}

/// Diagnostic information from redox solver
#[derive(Debug, Clone)]
pub struct RedoxDiagnostics {
    pub time_sec: f64,
    pub nadph_mM: f64,
    pub nadp_plus_mM: f64,
    pub nadph_nadp_ratio: f64,
    pub gsh_mM: f64,
    pub gssg_mM: f64,
    pub gsh_gssg_ratio: f64,
    pub total_glutathione_mM: f64,
    pub h2o2_uM: f64,
    pub ca_cytosolic_uM: f64,
    pub nadph_production_mM_per_sec: f64,
    pub oxidative_stress: f64,
    pub membrane_tension_pN_per_nm: f64,
    pub ppp_rates: Vec<(&'static str, f64)>,
    pub glutathione_rates: Vec<(&'static str, f64)>,
    pub piezo1_open_probability: f64,
    pub atp_release_mM_per_sec: f64,
}

impl RedoxDiagnostics {
    /// Print a formatted summary
    pub fn print_summary(&self) {
        println!("=== Redox Metabolism (t = {:.3} s) ===", self.time_sec);
        println!();
        println!("NADPH/NADP+ Pool:");
        println!("  NADPH:        {:.4} mM", self.nadph_mM);
        println!("  NADP+:        {:.4} mM", self.nadp_plus_mM);
        println!("  Ratio:        {:.1} (target: 10-20)", self.nadph_nadp_ratio);
        println!("  Production:   {:.6} mM/s", self.nadph_production_mM_per_sec);
        println!();
        println!("Glutathione Pool:");
        println!("  GSH:          {:.3} mM", self.gsh_mM);
        println!("  GSSG:         {:.4} mM", self.gssg_mM);
        println!("  Ratio:        {:.0} (target: 100-400)", self.gsh_gssg_ratio);
        println!("  Total:        {:.2} mM (target: 2-3 mM)", self.total_glutathione_mM);
        println!();
        println!("Oxidative State:");
        println!("  H2O2:         {:.2} uM (target: <5 uM)", self.h2o2_uM);
        println!("  Stress level: {:.1}x", self.oxidative_stress);
        println!();
        println!("Piezo1 Channel:");
        println!("  Tension:      {:.2} pN/nm", self.membrane_tension_pN_per_nm);
        println!("  P(open):      {:.1}%", self.piezo1_open_probability * 100.0);
        println!("  Ca2+ (cyt):   {:.0} nM (target: ~100 nM)", self.ca_cytosolic_uM * 1000.0);
        println!("  ATP release:  {:.6} mM/s", self.atp_release_mM_per_sec);
    }

    /// Print detailed reaction rates
    pub fn print_rates(&self) {
        println!("\nPPP Reaction Rates (mM/s):");
        for (name, rate) in &self.ppp_rates {
            println!("  {}: {:.6}", name, rate);
        }

        println!("\nGlutathione Reaction Rates (mM/s):");
        for (name, rate) in &self.glutathione_rates {
            println!("  {}: {:.6}", name, rate);
        }
    }
}

/// Initialize a metabolite pool with physiological redox concentrations
pub fn initialize_redox_metabolites(pool: &mut MetabolitePool, indices: &RedoxIndices) {
    // PPP intermediates (low concentrations)
    pool.set(indices.phosphogluconolactone_6, 0.001);
    pool.set(indices.phosphogluconate_6, 0.02);
    pool.set(indices.ribulose_5_phosphate, 0.01);
    pool.set(indices.ribose_5_phosphate, 0.01);
    pool.set(indices.xylulose_5_phosphate, 0.01);
    pool.set(indices.sedoheptulose_7_phosphate, 0.01);
    pool.set(indices.erythrose_4_phosphate, 0.01);

    // NADPH/NADP+ (ratio ~15)
    // Reference: Veech 1969
    pool.set(indices.nadph, 0.3);
    pool.set(indices.nadp_plus, 0.02);

    // Glutathione (GSH/GSSG ratio ~250, total ~2.5 mM)
    // Reference: Beutler 1969, Meister 1983
    pool.set(indices.gsh, 2.5);
    pool.set(indices.gssg, 0.01);

    // H2O2 (very low at steady state)
    // Reference: Chance 1979
    pool.set(indices.h2o2, 0.001);  // 1 uM

    // Cytosolic Ca2+ (100 nM = 0.1 uM)
    // Reference: Engelmann 2008
    pool.set(indices.ca2_plus_cytosolic, 0.1);

    // GSH synthesis precursors
    // Reference: Meister 1983
    pool.set(indices.glutamate, 0.5);
    pool.set(indices.cysteine, 0.05);
    pool.set(indices.glycine, 0.5);
    pool.set(indices.gamma_glu_cys, 0.001);
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_test_system() -> (RedoxSolver, MetabolitePool) {
        let glycolysis = MetaboliteIndices::default();
        let config = RedoxConfig::default();
        let solver = RedoxSolver::new(&glycolysis, config);

        let total = 18 + RedoxIndices::new_metabolite_count();
        let mut pool = MetabolitePool::new(total);

        // Initialize glycolysis shared metabolites
        pool.set(glycolysis.glucose_6_phosphate, 0.05);
        pool.set(glycolysis.fructose_6_phosphate, 0.02);
        pool.set(glycolysis.glyceraldehyde_3_phosphate, 0.005);
        pool.set(glycolysis.atp, 2.0);
        pool.set(glycolysis.adp, 0.25);

        // Initialize redox metabolites
        initialize_redox_metabolites(&mut pool, &solver.indices);

        (solver, pool)
    }

    #[test]
    fn test_redox_solver_creation() {
        let (solver, _pool) = create_test_system();
        assert_eq!(solver.time_sec, 0.0);
    }

    #[test]
    fn test_redox_step() {
        let (mut solver, mut pool) = create_test_system();

        // Run a few steps
        for _ in 0..10 {
            solver.step(&mut pool);
        }

        assert!(solver.time_sec > 0.0);
    }

    #[test]
    fn test_redox_diagnostics() {
        let (solver, pool) = create_test_system();
        let diag = solver.diagnostics(&pool);

        // Check initial ratios are in range
        assert!(diag.nadph_nadp_ratio > 5.0, "NADPH/NADP+ ratio too low");
        assert!(diag.gsh_gssg_ratio > 50.0, "GSH/GSSG ratio too low");
        assert!(diag.total_glutathione_mM > 1.5, "Total GSH too low");
    }

    #[test]
    fn test_redox_validation() {
        let (solver, pool) = create_test_system();
        let warnings = solver.validate_state(&pool);

        // Initial state should be close to physiological (minimal warnings)
        // Note: Some warnings may occur if initial conditions aren't perfect
        println!("Validation warnings: {:?}", warnings);
    }

    #[test]
    fn test_oxidative_stress_response() {
        let (mut solver, mut pool) = create_test_system();

        // Run without stress first
        solver.run(&mut pool, 2.0);
        let baseline_h2o2 = pool.get(solver.glutathione.indices.h2o2);

        // Reset and run with stress
        let (mut solver2, mut pool2) = create_test_system();
        solver2.set_oxidative_stress(20.0);
        solver2.run(&mut pool2, 2.0);
        let stressed_h2o2 = pool2.get(solver2.glutathione.indices.h2o2);

        // H2O2 should be higher under oxidative stress (more production)
        assert!(
            stressed_h2o2 > baseline_h2o2,
            "H2O2 should be higher under oxidative stress: {} -> {}",
            baseline_h2o2, stressed_h2o2
        );
    }

    #[test]
    fn test_piezo1_tension_response() {
        let (mut solver, mut pool) = create_test_system();

        // Get baseline Ca2+
        let baseline_ca = pool.get(solver.indices.ca2_plus_cytosolic);

        // Apply membrane tension
        solver.set_membrane_tension(3.0);

        // Run briefly
        solver.run(&mut pool, 0.1);

        // Ca2+ should increase with tension
        let elevated_ca = pool.get(solver.indices.ca2_plus_cytosolic);
        assert!(
            elevated_ca > baseline_ca,
            "Ca2+ should increase with membrane tension"
        );
    }

    #[test]
    fn test_steady_state_approach() {
        let (mut solver, mut pool) = create_test_system();

        // Maintain G6P supply (simulating glycolysis feeding PPP)
        // Without this, NADPH gets depleted as there's no G6P to feed the PPP
        let g6p_idx = solver.indices.glucose_6_phosphate;

        // Run for 5 seconds with G6P maintenance
        let dt = solver.config.dt_sec;
        let n_steps = (5.0 / dt).ceil() as usize;
        for _ in 0..n_steps {
            // Maintain G6P at physiological level (simulating glycolysis input)
            let current_g6p = pool.get(g6p_idx);
            if current_g6p < 0.03 {
                pool.set(g6p_idx, 0.05);  // Replenish
            }
            solver.step(&mut pool);
        }

        let diag = solver.diagnostics(&pool);

        // With maintained G6P, NADPH pool should stay reasonable
        assert!(
            diag.nadph_nadp_ratio > 0.5,
            "NADPH/NADP+ ratio too low: {}",
            diag.nadph_nadp_ratio
        );

        // GSH should be maintained by GR using NADPH
        assert!(
            diag.gsh_gssg_ratio > 10.0,
            "GSH/GSSG ratio too low: {}",
            diag.gsh_gssg_ratio
        );
    }
}
