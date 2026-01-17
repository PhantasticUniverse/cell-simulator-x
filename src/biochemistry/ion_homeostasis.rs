//! Ion Homeostasis for RBC simulation.
//!
//! This module implements ion transport mechanisms that maintain physiological
//! ion gradients in red blood cells:
//! - Na+/K+-ATPase: 3 Na+ out, 2 K+ in, 1 ATP consumed
//! - Passive ion leaks (balanced with pump at steady state)
//!
//! Physiological targets:
//! - Cytosolic Na+: 5-15 mM (Bernstein 1954)
//! - Cytosolic K+: 140-150 mM (Bernstein 1954)
//! - Na/K-ATPase ATP consumption: 0.01-0.05 mM/s (Garrahan & Glynn 1967)
//!
//! References:
//! - Bernstein RE. J Clin Invest. 1954;33:1-6 (RBC ion concentrations)
//! - Garrahan PJ, Glynn IM. J Physiol. 1967;192:159-174 (Na/K pump kinetics)
//! - Hoffman JF. Am J Med. 1966;41:666-680 (RBC cation transport)
//! - Tosteson DC, Hoffman JF. J Gen Physiol. 1960;44:169-194 (Na/K pump stoichiometry)

use super::MetabolitePool;
use super::glycolysis::MetaboliteIndices;
use super::pentose_phosphate::RedoxIndices;

/// Indices for ion metabolites in the concentration vector
#[derive(Debug, Clone, Copy)]
pub struct IonIndices {
    // Ion concentrations (new metabolites starting at index 35)
    /// Cytosolic Na+ (mM)
    pub na_plus_cytosolic: usize,
    /// Cytosolic K+ (mM)
    pub k_plus_cytosolic: usize,
    /// Cytosolic Cl- (mM)
    pub cl_minus_cytosolic: usize,

    // Shared references from other modules
    /// Cytosolic Ca2+ (from RedoxIndices)
    pub ca2_plus_cytosolic: usize,
    /// ATP (from glycolysis)
    pub atp: usize,
    /// ADP (from glycolysis)
    pub adp: usize,
}

impl IonIndices {
    /// Create ion indices starting at the given index
    ///
    /// # Arguments
    /// * `glycolysis` - Glycolysis indices for ATP/ADP references
    /// * `redox` - Redox indices for Ca2+ reference
    /// * `start_idx` - Starting index for new ion metabolites
    pub fn new(glycolysis: &MetaboliteIndices, redox: &RedoxIndices, start_idx: usize) -> Self {
        Self {
            // New ion metabolites
            na_plus_cytosolic: start_idx,
            k_plus_cytosolic: start_idx + 1,
            cl_minus_cytosolic: start_idx + 2,

            // Shared references
            ca2_plus_cytosolic: redox.ca2_plus_cytosolic,
            atp: glycolysis.atp,
            adp: glycolysis.adp,
        }
    }

    /// Number of new metabolites added by ion module
    pub fn new_metabolite_count() -> usize {
        3  // Na+, K+, Cl-
    }
}

/// Configuration for ion homeostasis system
#[derive(Debug, Clone)]
pub struct IonHomeostasisConfig {
    // External ion concentrations (mM)
    /// External Na+ concentration (plasma ~140 mM)
    pub na_external_mM: f64,
    /// External K+ concentration (plasma ~5 mM)
    pub k_external_mM: f64,
    /// External Cl- concentration (plasma ~100 mM)
    pub cl_external_mM: f64,

    // Initial cytosolic concentrations (mM)
    /// Initial cytosolic Na+ (~10 mM)
    pub na_initial_mM: f64,
    /// Initial cytosolic K+ (~140 mM)
    pub k_initial_mM: f64,
    /// Initial cytosolic Cl- (~80 mM)
    pub cl_initial_mM: f64,

    // Passive leak conductances (mM/s per mM gradient)
    /// Na+ leak conductance
    pub g_na_per_sec: f64,
    /// K+ leak conductance
    pub g_k_per_sec: f64,
}

impl Default for IonHomeostasisConfig {
    fn default() -> Self {
        Self {
            // External concentrations (plasma values)
            // Reference: Bernstein 1954
            na_external_mM: 140.0,
            k_external_mM: 5.0,
            cl_external_mM: 100.0,

            // Initial cytosolic concentrations
            // Reference: Bernstein 1954
            na_initial_mM: 10.0,
            k_initial_mM: 140.0,
            cl_initial_mM: 80.0,

            // Passive leak conductances
            // Tuned to balance pump at steady state with pump rate ~0.01 mM/s
            // At Na+=10mM: pump efflux = 3*0.01 = 0.03 mM/s
            //              leak influx = g_na * (140-10) = g_na * 130
            //              For balance: g_na = 0.03/130 ≈ 0.00023
            // At K+=140mM: pump influx = 2*0.01 = 0.02 mM/s
            //              leak efflux = g_k * (140-5) = g_k * 135
            //              For balance: g_k = 0.02/135 ≈ 0.00015
            // Reference: Hoffman 1966
            g_na_per_sec: 0.00024,  // Low Na+ permeability (balanced with pump)
            g_k_per_sec: 0.00015,   // K+ permeability (balanced with pump)
        }
    }
}

/// Na+/K+-ATPase pump
///
/// Transports 3 Na+ out and 2 K+ in per ATP hydrolyzed.
/// This is the primary mechanism maintaining ion gradients in RBCs.
///
/// Rate equation uses Hill kinetics for substrate cooperativity:
/// v = Vmax * (Na_i^3 / (Km_Na^3 + Na_i^3)) * (K_e^2 / (Km_K^2 + K_e^2)) * (ATP / (Km_ATP + ATP))
///
/// Reference: Garrahan PJ, Glynn IM. J Physiol. 1967;192:159-174
#[derive(Debug, Clone)]
pub struct NaKATPase {
    /// Maximum pump rate (mM ATP/s)
    /// Reference: Garrahan & Glynn 1967, ~0.03 mM/s in RBCs
    pub vmax_mM_per_sec: f64,

    /// Km for intracellular Na+ (mM)
    /// Reference: Garrahan & Glynn 1967
    pub km_na_mM: f64,

    /// Km for extracellular K+ (mM)
    /// Reference: Garrahan & Glynn 1967
    pub km_k_mM: f64,

    /// Km for ATP (mM)
    /// Reference: Garrahan & Glynn 1967
    pub km_atp_mM: f64,

    /// Hill coefficient for Na+ (3 Na+ per cycle)
    pub hill_na: f64,

    /// Hill coefficient for K+ (2 K+ per cycle)
    pub hill_k: f64,

    /// Enable pump (set false to simulate ouabain inhibition)
    pub enabled: bool,
}

impl Default for NaKATPase {
    fn default() -> Self {
        Self {
            // Reference: Garrahan & Glynn 1967, Hoffman 1966
            // Vmax tuned to achieve pump rate of 0.01-0.05 mM/s at physiological conditions
            vmax_mM_per_sec: 0.055,  // Tunable to achieve target pump rate
            km_na_mM: 15.0,          // Km for cytosolic Na+
            km_k_mM: 1.5,            // Km for extracellular K+
            km_atp_mM: 0.2,          // Km for ATP
            hill_na: 3.0,            // Stoichiometry: 3 Na+
            hill_k: 2.0,             // Stoichiometry: 2 K+
            enabled: true,
        }
    }
}

impl NaKATPase {
    /// Calculate pump rate (mM ATP consumed per second)
    ///
    /// # Arguments
    /// * `na_cytosolic_mM` - Intracellular Na+ concentration
    /// * `k_external_mM` - Extracellular K+ concentration
    /// * `atp_mM` - ATP concentration
    ///
    /// # Returns
    /// Pump rate in mM/s (ATP consumption rate)
    pub fn pump_rate_mM_per_sec(
        &self,
        na_cytosolic_mM: f64,
        k_external_mM: f64,
        atp_mM: f64,
    ) -> f64 {
        if !self.enabled {
            return 0.0;
        }

        if na_cytosolic_mM <= 0.0 || k_external_mM <= 0.0 || atp_mM <= 0.0 {
            return 0.0;
        }

        // Na+ term with Hill kinetics
        let na_n = na_cytosolic_mM.powf(self.hill_na);
        let km_na_n = self.km_na_mM.powf(self.hill_na);
        let na_term = na_n / (km_na_n + na_n);

        // K+ term with Hill kinetics
        let k_n = k_external_mM.powf(self.hill_k);
        let km_k_n = self.km_k_mM.powf(self.hill_k);
        let k_term = k_n / (km_k_n + k_n);

        // ATP term (Michaelis-Menten)
        let atp_term = atp_mM / (self.km_atp_mM + atp_mM);

        self.vmax_mM_per_sec * na_term * k_term * atp_term
    }

    /// Compute derivatives for all affected metabolites
    ///
    /// Per pump cycle:
    /// - 3 Na+ exported (cytosolic Na+ decreases)
    /// - 2 K+ imported (cytosolic K+ increases)
    /// - 1 ATP consumed
    /// - 1 ADP produced
    ///
    /// # Arguments
    /// * `metabolites` - Current metabolite concentrations
    /// * `indices` - Ion indices
    /// * `k_external_mM` - Extracellular K+ concentration
    /// * `dydt` - Derivatives array to update
    pub fn compute_derivatives(
        &self,
        metabolites: &MetabolitePool,
        indices: &IonIndices,
        k_external_mM: f64,
        dydt: &mut [f64],
    ) {
        let na = metabolites.get(indices.na_plus_cytosolic);
        let atp = metabolites.get(indices.atp);

        let pump_rate = self.pump_rate_mM_per_sec(na, k_external_mM, atp);

        // Stoichiometry: 3 Na+ out, 2 K+ in, 1 ATP consumed
        dydt[indices.na_plus_cytosolic] -= 3.0 * pump_rate;
        dydt[indices.k_plus_cytosolic] += 2.0 * pump_rate;
        dydt[indices.atp] -= pump_rate;
        dydt[indices.adp] += pump_rate;
    }
}

/// Ion Homeostasis System
///
/// Combines Na+/K+-ATPase with passive leak channels to maintain
/// physiological ion gradients.
pub struct IonHomeostasisSystem {
    /// Configuration
    pub config: IonHomeostasisConfig,
    /// Na+/K+-ATPase pump
    pub na_k_pump: NaKATPase,
    /// Ion indices
    pub indices: IonIndices,
}

impl IonHomeostasisSystem {
    /// Create a new ion homeostasis system
    pub fn new(glycolysis: &MetaboliteIndices, redox: &RedoxIndices, start_idx: usize) -> Self {
        Self {
            config: IonHomeostasisConfig::default(),
            na_k_pump: NaKATPase::default(),
            indices: IonIndices::new(glycolysis, redox, start_idx),
        }
    }

    /// Create with custom configuration
    pub fn with_config(
        glycolysis: &MetaboliteIndices,
        redox: &RedoxIndices,
        start_idx: usize,
        config: IonHomeostasisConfig,
    ) -> Self {
        Self {
            config,
            na_k_pump: NaKATPase::default(),
            indices: IonIndices::new(glycolysis, redox, start_idx),
        }
    }

    /// Enable or disable the Na+/K+-ATPase (simulate ouabain treatment)
    pub fn set_pump_enabled(&mut self, enabled: bool) {
        self.na_k_pump.enabled = enabled;
    }

    /// Compute all derivatives for ion homeostasis
    ///
    /// Includes:
    /// - Na+/K+-ATPase pump (active transport)
    /// - Passive Na+ leak (inward)
    /// - Passive K+ leak (outward)
    pub fn compute_derivatives(&self, metabolites: &MetabolitePool, dydt: &mut [f64]) {
        let na_cyt = metabolites.get(self.indices.na_plus_cytosolic);
        let k_cyt = metabolites.get(self.indices.k_plus_cytosolic);

        // Na+/K+-ATPase pump
        self.na_k_pump.compute_derivatives(
            metabolites,
            &self.indices,
            self.config.k_external_mM,
            dydt,
        );

        // Passive Na+ leak (inward, down concentration gradient)
        // Na+ flows in because [Na+]ext > [Na+]cyt
        let na_gradient = self.config.na_external_mM - na_cyt;
        let na_leak = self.config.g_na_per_sec * na_gradient;
        dydt[self.indices.na_plus_cytosolic] += na_leak;

        // Passive K+ leak (outward, down concentration gradient)
        // K+ flows out because [K+]cyt > [K+]ext
        let k_gradient = k_cyt - self.config.k_external_mM;
        let k_leak = self.config.g_k_per_sec * k_gradient;
        dydt[self.indices.k_plus_cytosolic] -= k_leak;
    }

    /// Get diagnostic information
    pub fn diagnostics(&self, metabolites: &MetabolitePool) -> IonDiagnostics {
        let na_cyt = metabolites.get(self.indices.na_plus_cytosolic);
        let k_cyt = metabolites.get(self.indices.k_plus_cytosolic);
        let cl_cyt = metabolites.get(self.indices.cl_minus_cytosolic);
        let atp = metabolites.get(self.indices.atp);

        let pump_rate = self.na_k_pump.pump_rate_mM_per_sec(
            na_cyt,
            self.config.k_external_mM,
            atp,
        );

        // Calculate leak rates
        let na_gradient = self.config.na_external_mM - na_cyt;
        let na_leak = self.config.g_na_per_sec * na_gradient;

        let k_gradient = k_cyt - self.config.k_external_mM;
        let k_leak = self.config.g_k_per_sec * k_gradient;

        IonDiagnostics {
            na_cytosolic_mM: na_cyt,
            k_cytosolic_mM: k_cyt,
            cl_cytosolic_mM: cl_cyt,
            pump_rate_mM_per_sec: pump_rate,
            pump_atp_consumption_mM_per_sec: pump_rate,
            na_leak_mM_per_sec: na_leak,
            k_leak_mM_per_sec: k_leak,
            pump_enabled: self.na_k_pump.enabled,
        }
    }

    /// Validate ion concentrations against physiological targets
    pub fn validate_state(&self, metabolites: &MetabolitePool) -> Vec<String> {
        let mut warnings = Vec::new();
        let diag = self.diagnostics(metabolites);

        // Na+ (target: 5-15 mM)
        if diag.na_cytosolic_mM < 5.0 {
            warnings.push(format!(
                "Na+ too low: {:.1} mM (target: 5-15 mM)",
                diag.na_cytosolic_mM
            ));
        } else if diag.na_cytosolic_mM > 15.0 {
            warnings.push(format!(
                "Na+ too high: {:.1} mM (target: 5-15 mM)",
                diag.na_cytosolic_mM
            ));
        }

        // K+ (target: 140-150 mM)
        if diag.k_cytosolic_mM < 140.0 {
            warnings.push(format!(
                "K+ too low: {:.1} mM (target: 140-150 mM)",
                diag.k_cytosolic_mM
            ));
        } else if diag.k_cytosolic_mM > 150.0 {
            warnings.push(format!(
                "K+ too high: {:.1} mM (target: 140-150 mM)",
                diag.k_cytosolic_mM
            ));
        }

        // Na/K pump ATP consumption (target: 0.01-0.05 mM/s)
        if self.na_k_pump.enabled {
            if diag.pump_atp_consumption_mM_per_sec < 0.01 {
                warnings.push(format!(
                    "Na/K pump rate low: {:.4} mM/s (target: 0.01-0.05 mM/s)",
                    diag.pump_atp_consumption_mM_per_sec
                ));
            } else if diag.pump_atp_consumption_mM_per_sec > 0.05 {
                warnings.push(format!(
                    "Na/K pump rate high: {:.4} mM/s (target: 0.01-0.05 mM/s)",
                    diag.pump_atp_consumption_mM_per_sec
                ));
            }
        }

        warnings
    }
}

/// Diagnostic information for ion homeostasis
#[derive(Debug, Clone)]
pub struct IonDiagnostics {
    /// Cytosolic Na+ concentration (mM)
    pub na_cytosolic_mM: f64,
    /// Cytosolic K+ concentration (mM)
    pub k_cytosolic_mM: f64,
    /// Cytosolic Cl- concentration (mM)
    pub cl_cytosolic_mM: f64,
    /// Na+/K+-ATPase pump rate (mM/s)
    pub pump_rate_mM_per_sec: f64,
    /// ATP consumption by pump (mM/s)
    pub pump_atp_consumption_mM_per_sec: f64,
    /// Na+ leak rate (mM/s, positive = inward)
    pub na_leak_mM_per_sec: f64,
    /// K+ leak rate (mM/s, positive = outward)
    pub k_leak_mM_per_sec: f64,
    /// Whether pump is enabled
    pub pump_enabled: bool,
}

impl IonDiagnostics {
    /// Print a formatted summary
    pub fn print_summary(&self) {
        println!("=== Ion Homeostasis ===");
        println!("Cytosolic concentrations:");
        println!("  Na+:  {:.1} mM (target: 5-15 mM)", self.na_cytosolic_mM);
        println!("  K+:   {:.1} mM (target: 140-150 mM)", self.k_cytosolic_mM);
        println!("  Cl-:  {:.1} mM", self.cl_cytosolic_mM);
        println!();
        println!("Na+/K+-ATPase:");
        println!("  Enabled:          {}", if self.pump_enabled { "Yes" } else { "No (ouabain)" });
        println!("  Pump rate:        {:.4} mM/s", self.pump_rate_mM_per_sec);
        println!("  ATP consumption:  {:.4} mM/s (target: 0.01-0.05 mM/s)", self.pump_atp_consumption_mM_per_sec);
        println!();
        println!("Passive leaks:");
        println!("  Na+ leak (in):    {:.4} mM/s", self.na_leak_mM_per_sec);
        println!("  K+ leak (out):    {:.4} mM/s", self.k_leak_mM_per_sec);
    }
}

/// Initialize ion metabolites in the pool
pub fn initialize_ion_metabolites(metabolites: &mut MetabolitePool, indices: &IonIndices, config: &IonHomeostasisConfig) {
    metabolites.set(indices.na_plus_cytosolic, config.na_initial_mM);
    metabolites.set(indices.k_plus_cytosolic, config.k_initial_mM);
    metabolites.set(indices.cl_minus_cytosolic, config.cl_initial_mM);
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_test_system() -> (IonHomeostasisSystem, MetabolitePool) {
        let glycolysis = MetaboliteIndices::default();
        let redox = RedoxIndices::new(&glycolysis, 18);
        let system = IonHomeostasisSystem::new(&glycolysis, &redox, 35);

        // Create pool with all metabolites (38 total)
        let mut pool = MetabolitePool::new(38);

        // Set ATP/ADP
        pool.set(glycolysis.atp, 2.0);
        pool.set(glycolysis.adp, 0.25);

        // Initialize ion concentrations
        initialize_ion_metabolites(&mut pool, &system.indices, &system.config);

        (system, pool)
    }

    #[test]
    fn test_ion_indices() {
        let glycolysis = MetaboliteIndices::default();
        let redox = RedoxIndices::new(&glycolysis, 18);
        let indices = IonIndices::new(&glycolysis, &redox, 35);

        assert_eq!(indices.na_plus_cytosolic, 35);
        assert_eq!(indices.k_plus_cytosolic, 36);
        assert_eq!(indices.cl_minus_cytosolic, 37);
        assert_eq!(IonIndices::new_metabolite_count(), 3);
    }

    #[test]
    fn test_na_k_pump_basic() {
        let pump = NaKATPase::default();

        // Test with physiological concentrations
        let rate = pump.pump_rate_mM_per_sec(10.0, 5.0, 2.0);
        assert!(rate > 0.0, "Pump should be active");
        assert!(rate < pump.vmax_mM_per_sec, "Rate should be below Vmax");
    }

    #[test]
    fn test_na_k_pump_na_dependence() {
        let pump = NaKATPase::default();

        let rate_low_na = pump.pump_rate_mM_per_sec(5.0, 5.0, 2.0);
        let rate_high_na = pump.pump_rate_mM_per_sec(20.0, 5.0, 2.0);

        assert!(rate_high_na > rate_low_na, "Pump rate should increase with Na+");
    }

    #[test]
    fn test_na_k_pump_disabled() {
        let mut pump = NaKATPase::default();
        pump.enabled = false;

        let rate = pump.pump_rate_mM_per_sec(10.0, 5.0, 2.0);
        assert_eq!(rate, 0.0, "Disabled pump should have zero rate");
    }

    #[test]
    fn test_system_derivatives() {
        let (system, pool) = create_test_system();
        let mut dydt = vec![0.0; 38];

        system.compute_derivatives(&pool, &mut dydt);

        // Na+ should decrease (pump > leak at 10 mM)
        // K+ should change (pump adds, leak removes)
        // ATP should decrease
        assert!(dydt[system.indices.atp] < 0.0, "ATP should be consumed");
        assert!(dydt[system.indices.adp] > 0.0, "ADP should be produced");
    }

    #[test]
    fn test_diagnostics() {
        let (system, pool) = create_test_system();
        let diag = system.diagnostics(&pool);

        assert!((diag.na_cytosolic_mM - 10.0).abs() < 0.1);
        assert!((diag.k_cytosolic_mM - 140.0).abs() < 0.1);
        assert!(diag.pump_rate_mM_per_sec > 0.0);
    }

    #[test]
    fn test_ouabain_simulation() {
        let (mut system, pool) = create_test_system();

        // Normal pump rate
        let diag_normal = system.diagnostics(&pool);
        assert!(diag_normal.pump_rate_mM_per_sec > 0.0);

        // Disable pump (simulate ouabain)
        system.set_pump_enabled(false);
        let diag_ouabain = system.diagnostics(&pool);
        assert_eq!(diag_ouabain.pump_rate_mM_per_sec, 0.0);
    }

    #[test]
    fn test_validation() {
        let (system, pool) = create_test_system();
        let warnings = system.validate_state(&pool);

        // Initial state should be within targets
        assert!(
            warnings.is_empty(),
            "Initial state should be valid: {:?}",
            warnings
        );
    }
}
