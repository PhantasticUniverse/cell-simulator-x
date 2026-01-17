//! Sickle Cell Disease Model.
//!
//! Simulates the effects of hemoglobin S (HbS) on RBC oxygen binding and mechanics:
//! - Modified Adair constants (shifted P50)
//! - Polymerization kinetics at low O2 saturation
//! - Stiffness increase during sickling (future: mechanics coupling)
//! - Chronic oxidative stress from repeated sickling
//!
//! Genotype effects:
//! - HbAA (normal): P50 = 26.8 mmHg, no polymerization
//! - HbAS (trait): P50 ~28 mmHg, minimal polymerization
//! - HbSS (disease): P50 = 31 mmHg, significant polymerization at low pO2
//!
//! Polymerization occurs when:
//! - O2 saturation drops below critical threshold (~35%)
//! - Delay time inversely related to HbS concentration (15th power dependence)
//!
//! Validation targets (Eaton 1987):
//! | Genotype | P50 (mmHg) | Polymerization at pO2=40 |
//! |----------|------------|--------------------------|
//! | HbAA     | 26.8       | None                     |
//! | HbAS     | 28.0       | Minimal                  |
//! | HbSS     | 31.0       | Significant              |
//!
//! References:
//! - Eaton WA, Hofrichter J. Blood. 1987;70:1245-1266
//! - Bunn HF. N Engl J Med. 1997;337:762-769
//! - Ferrone FA. J Mol Biol. 2004;341:843-855

use super::{DiseaseModel, DiseaseDiagnostics};
use crate::biochemistry::{
    FullyIntegratedConfig, FullyIntegratedSolver, FullyIntegratedIndices,
    MetabolitePool, AdairConstants,
};

/// Configuration for sickle cell model.
#[derive(Debug, Clone)]
pub struct SickleCellConfig {
    /// P50 for HbS at standard conditions (mmHg).
    /// Reference: Eaton 1987 - HbS has lower O2 affinity
    pub hbs_p50_mmHg: f64,

    /// Hill coefficient for HbS.
    /// Reference: Slightly lower cooperativity than HbA
    pub hbs_hill_coefficient: f64,

    /// P50 for HbA (normal) for comparison
    pub hba_p50_mmHg: f64,

    /// Critical O2 saturation below which polymerization occurs.
    /// Reference: Eaton 1987 - sickling begins below ~35% saturation
    pub polymerization_critical_saturation: f64,

    /// Polymerization rate constant (per second).
    /// Rate of polymer formation once nucleation occurs
    pub polymerization_rate_per_sec: f64,

    /// Depolymerization rate constant (per second).
    /// Rate of polymer dissolution upon reoxygenation
    pub depolymerization_rate_per_sec: f64,

    /// Delay time exponent (cooperative nucleation).
    /// Reference: Ferrone 2004 - 15th power dependence on HbS concentration
    pub delay_time_exponent: f64,

    /// Stiffness multiplier at full polymerization.
    /// For future mechanics coupling (sickled cells are ~10x stiffer)
    pub stiffness_multiplier_full_polymer: f64,

    /// Chronic oxidative stress from sickling cycles.
    /// Reference: Bunn 1997 - repeated sickling damages membrane
    pub chronic_oxidative_stress: f64,

    /// Additional NADPH consumption from membrane repair
    pub nadph_consumption_sickling: f64,
}

impl Default for SickleCellConfig {
    fn default() -> Self {
        Self {
            hbs_p50_mmHg: 31.0,     // Eaton 1987
            hbs_hill_coefficient: 2.6,  // Slightly reduced cooperativity
            hba_p50_mmHg: 26.8,    // Normal HbA
            polymerization_critical_saturation: 0.35,  // Eaton 1987
            polymerization_rate_per_sec: 0.1,
            depolymerization_rate_per_sec: 0.5,  // Faster recovery
            delay_time_exponent: 15.0,  // Ferrone 2004
            stiffness_multiplier_full_polymer: 10.0,
            chronic_oxidative_stress: 1.8,  // Bunn 1997
            nadph_consumption_sickling: 0.0002,
        }
    }
}

/// Sickle Cell Disease Model.
///
/// Models hemoglobin S effects on oxygen binding and cell mechanics.
pub struct SickleCellModel {
    /// HbS fraction (0 = HbAA, 0.5 = HbAS, 1.0 = HbSS)
    pub hbs_fraction: f64,
    /// Current polymer fraction (0-1)
    pub polymer_fraction: f64,
    /// Configuration parameters
    pub config: SickleCellConfig,
    /// Whether currently in sickling episode
    is_sickling: bool,
    /// Cumulative sickling events (for chronic damage tracking)
    sickling_events: u32,
}

impl SickleCellModel {
    /// Create a new sickle cell model with given HbS fraction.
    ///
    /// # Arguments
    /// * `hbs_fraction` - Fraction of hemoglobin that is HbS
    ///   - 0.0 = HbAA (normal)
    ///   - 0.3-0.5 = HbAS (sickle trait)
    ///   - 0.8-1.0 = HbSS (sickle cell disease)
    pub fn new(hbs_fraction: f64) -> Self {
        Self {
            hbs_fraction: hbs_fraction.clamp(0.0, 1.0),
            polymer_fraction: 0.0,
            config: SickleCellConfig::default(),
            is_sickling: false,
            sickling_events: 0,
        }
    }

    /// Create with custom configuration
    pub fn with_config(hbs_fraction: f64, config: SickleCellConfig) -> Self {
        Self {
            hbs_fraction: hbs_fraction.clamp(0.0, 1.0),
            polymer_fraction: 0.0,
            config,
            is_sickling: false,
            sickling_events: 0,
        }
    }

    /// Calculate effective P50 based on HbS fraction.
    ///
    /// Linear interpolation between HbA and HbS P50 values.
    pub fn effective_p50_mmHg(&self) -> f64 {
        let p50_a = self.config.hba_p50_mmHg;
        let p50_s = self.config.hbs_p50_mmHg;
        p50_a + (p50_s - p50_a) * self.hbs_fraction
    }

    /// Calculate effective Hill coefficient.
    pub fn effective_hill_coefficient(&self) -> f64 {
        // Linear interpolation: 2.7 for HbA, 2.6 for HbS
        2.7 + (self.config.hbs_hill_coefficient - 2.7) * self.hbs_fraction
    }

    /// Get Adair constants for the mixed Hb population.
    pub fn effective_adair_constants(&self) -> AdairConstants {
        AdairConstants::from_p50_and_hill(
            self.effective_p50_mmHg(),
            self.effective_hill_coefficient(),
        )
    }

    /// Calculate disease severity based on genotype.
    ///
    /// 0 = normal, 0.3 = trait carrier, 1.0 = homozygous disease
    pub fn severity(&self) -> f64 {
        // Nonlinear: trait (0.5 HbS) is relatively mild
        if self.hbs_fraction < 0.5 {
            self.hbs_fraction * 0.6  // Trait carriers have milder disease
        } else {
            0.3 + (self.hbs_fraction - 0.5) * 1.4  // Disease severity increases
        }
    }

    /// Get genotype classification
    pub fn genotype(&self) -> &'static str {
        if self.hbs_fraction < 0.1 {
            "HbAA (Normal)"
        } else if self.hbs_fraction < 0.6 {
            "HbAS (Sickle Trait)"
        } else if self.hbs_fraction < 0.95 {
            "HbSC or HbS/beta-thal"
        } else {
            "HbSS (Sickle Cell Disease)"
        }
    }

    /// Calculate polymerization propensity based on saturation.
    ///
    /// Returns a value from 0 (no polymerization) to 1 (full polymerization).
    pub fn polymerization_propensity(&self, saturation: f64) -> f64 {
        if self.hbs_fraction < 0.3 {
            // Trait carriers don't sickle under normal conditions
            return 0.0;
        }

        if saturation > self.config.polymerization_critical_saturation {
            return 0.0;
        }

        // Below critical saturation, polymerization depends on:
        // 1. How far below threshold (more deoxygenation = more polymer)
        // 2. HbS concentration (15th power dependence on delay time)
        let desat_factor = (self.config.polymerization_critical_saturation - saturation)
            / self.config.polymerization_critical_saturation;

        // Cooperative nucleation - higher HbS = faster polymerization
        let hbs_factor = self.hbs_fraction.powf(2.0);  // Simplified from 15th power for stability

        (desat_factor * hbs_factor).min(1.0)
    }

    /// Calculate current stiffness multiplier based on polymer fraction.
    pub fn stiffness_multiplier(&self) -> f64 {
        1.0 + (self.config.stiffness_multiplier_full_polymer - 1.0) * self.polymer_fraction
    }

    /// Update polymerization state based on current oxygen saturation.
    pub fn update_polymerization(&mut self, saturation: f64, dt_sec: f64) {
        let propensity = self.polymerization_propensity(saturation);

        if propensity > 0.0 {
            // Polymerize
            let rate = self.config.polymerization_rate_per_sec * propensity;
            let new_polymer = (self.polymer_fraction + rate * dt_sec).min(1.0);

            // Track sickling event
            if !self.is_sickling && new_polymer > 0.1 {
                self.is_sickling = true;
                self.sickling_events += 1;
            }

            self.polymer_fraction = new_polymer;
        } else {
            // Depolymerize (reoxygenation)
            let rate = self.config.depolymerization_rate_per_sec;
            self.polymer_fraction = (self.polymer_fraction - rate * dt_sec).max(0.0);

            if self.polymer_fraction < 0.05 {
                self.is_sickling = false;
            }
        }
    }

    /// Calculate effective oxidative stress.
    ///
    /// Chronic stress from repeated sickling cycles.
    pub fn effective_oxidative_stress(&self) -> f64 {
        // Base chronic stress from HbS
        let chronic = 1.0 + (self.config.chronic_oxidative_stress - 1.0) * self.hbs_fraction;

        // Additional acute stress during sickling
        let acute = if self.is_sickling { 0.5 } else { 0.0 };

        chronic + acute * self.polymer_fraction
    }
}

impl DiseaseModel for SickleCellModel {
    fn name(&self) -> &'static str {
        "Sickle Cell"
    }

    fn description(&self) -> String {
        format!(
            "{}: {:.0}% HbS, P50={:.1}mmHg",
            self.genotype(),
            self.hbs_fraction * 100.0,
            self.effective_p50_mmHg()
        )
    }

    fn modify_config(&self, config: &mut FullyIntegratedConfig) {
        // Set oxidative stress
        config.oxidative_stress_multiplier *= self.effective_oxidative_stress();
    }

    fn apply_time_effects(
        &mut self,
        solver: &mut FullyIntegratedSolver,
        _metabolites: &mut MetabolitePool,
        _elapsed_sec: f64,
    ) {
        // Get current saturation from hemoglobin state
        let saturation = solver.hb_state.saturation;

        // Update polymerization state
        self.update_polymerization(saturation, solver.config.dt_sec);

        // Update hemoglobin solver with modified constants
        // This shifts the oxygen equilibrium curve
        let effective_p50 = self.effective_p50_mmHg();
        let effective_hill = self.effective_hill_coefficient();

        // Modify the hemoglobin solver's base P50
        solver.hemoglobin.base_p50_mmHg = effective_p50;
        solver.hemoglobin.base_hill_n = effective_hill;
        solver.hemoglobin.base_constants = self.effective_adair_constants();

        // Update oxidative stress
        solver.glutathione.set_oxidative_stress(self.effective_oxidative_stress());
    }

    fn modify_derivatives(
        &self,
        state: &[f64],
        dydt: &mut [f64],
        indices: &FullyIntegratedIndices,
    ) {
        if self.hbs_fraction < 0.1 {
            return;  // No effect for normal genotype
        }

        // Additional NADPH consumption from membrane repair during sickling
        if self.polymer_fraction > 0.0 {
            let nadph = state[indices.redox.nadph];
            let consumption = self.config.nadph_consumption_sickling
                * self.polymer_fraction
                * nadph / (0.1 + nadph);
            dydt[indices.redox.nadph] -= consumption;
            dydt[indices.redox.nadp_plus] += consumption;
        }

        // Enhanced GSH oxidation during sickling (membrane damage)
        if self.is_sickling {
            let gsh = state[indices.redox.gsh];
            let oxidation = 0.0005 * self.polymer_fraction * gsh / (1.0 + gsh);
            dydt[indices.redox.gsh] -= 2.0 * oxidation;
            dydt[indices.redox.gssg] += oxidation;
        }

        // Chronic ATP consumption from membrane repair
        let atp = state[indices.glycolysis.atp];
        let chronic_consumption = 0.00005 * self.hbs_fraction * atp / (0.5 + atp);
        dydt[indices.glycolysis.atp] -= chronic_consumption;
        dydt[indices.glycolysis.adp] += chronic_consumption;
    }

    fn diagnostics(&self, metabolites: &MetabolitePool) -> DiseaseDiagnostics {
        use crate::biochemistry::FullyIntegratedIndices;
        let indices = FullyIntegratedIndices::new();

        let mut diag = DiseaseDiagnostics::new(self.name());
        diag.severity = self.severity();

        // Current values
        let atp = metabolites.get(indices.glycolysis.atp);
        let nadph = metabolites.get(indices.redox.nadph);
        let nadp = metabolites.get(indices.redox.nadp_plus);
        let gsh = metabolites.get(indices.redox.gsh);
        let gssg = metabolites.get(indices.redox.gssg);

        let nadph_nadp_ratio = if nadp > 1e-9 { nadph / nadp } else { f64::INFINITY };
        let gsh_gssg_ratio = if gssg > 1e-9 { gsh / gssg } else { f64::INFINITY };

        // Add metrics
        diag.add_metric("hbs_fraction", self.hbs_fraction);
        diag.add_metric("effective_p50_mmHg", self.effective_p50_mmHg());
        diag.add_metric("effective_hill_n", self.effective_hill_coefficient());
        diag.add_metric("polymer_fraction", self.polymer_fraction);
        diag.add_metric("stiffness_multiplier", self.stiffness_multiplier());
        diag.add_metric("sickling_events", self.sickling_events as f64);
        diag.add_metric("atp_mM", atp);
        diag.add_metric("nadph_nadp_ratio", nadph_nadp_ratio);
        diag.add_metric("gsh_gssg_ratio", gsh_gssg_ratio);
        diag.add_metric("oxidative_stress", self.effective_oxidative_stress());

        // Status
        diag.add_status(self.genotype());
        diag.add_status(&format!("P50: {:.1} mmHg (normal: 26.8)", self.effective_p50_mmHg()));

        if self.is_sickling {
            diag.add_status(&format!(
                "SICKLING IN PROGRESS - {:.0}% polymerized",
                self.polymer_fraction * 100.0
            ));
        } else if self.polymer_fraction > 0.0 {
            diag.add_status(&format!(
                "Recovering - {:.0}% residual polymer",
                self.polymer_fraction * 100.0
            ));
        } else {
            diag.add_status("No active sickling");
        }

        // Warnings
        if self.hbs_fraction >= 0.5 {
            if self.polymer_fraction > 0.5 {
                diag.add_warning("Severe sickling - vaso-occlusive crisis risk");
            }
            if self.sickling_events > 10 {
                diag.add_warning("Multiple sickling cycles - cumulative membrane damage");
            }
        }

        if nadph_nadp_ratio < 5.0 {
            diag.add_warning("NADPH depleted - oxidative damage increasing");
        }
        if gsh_gssg_ratio < 50.0 {
            diag.add_warning("GSH depleted - antioxidant defense compromised");
        }

        diag
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sickle_cell_creation() {
        let model = SickleCellModel::new(1.0);
        assert!((model.hbs_fraction - 1.0).abs() < 0.01);
        assert_eq!(model.name(), "Sickle Cell");
    }

    #[test]
    fn test_hbs_fraction_clamping() {
        let low = SickleCellModel::new(-0.5);
        let high = SickleCellModel::new(1.5);

        assert_eq!(low.hbs_fraction, 0.0);
        assert_eq!(high.hbs_fraction, 1.0);
    }

    #[test]
    fn test_effective_p50() {
        let hbaa = SickleCellModel::new(0.0);
        let hbas = SickleCellModel::new(0.5);
        let hbss = SickleCellModel::new(1.0);

        assert!((hbaa.effective_p50_mmHg() - 26.8).abs() < 0.1);
        assert!((hbas.effective_p50_mmHg() - 28.9).abs() < 0.5);  // ~halfway
        assert!((hbss.effective_p50_mmHg() - 31.0).abs() < 0.1);
    }

    #[test]
    fn test_genotype_classification() {
        let hbaa = SickleCellModel::new(0.0);
        let hbas = SickleCellModel::new(0.4);
        let hbss = SickleCellModel::new(1.0);

        assert!(hbaa.genotype().contains("Normal"));
        assert!(hbas.genotype().contains("Trait"));
        assert!(hbss.genotype().contains("Disease"));
    }

    #[test]
    fn test_polymerization_propensity_normal() {
        let hbaa = SickleCellModel::new(0.0);

        // Normal Hb should never polymerize
        assert_eq!(hbaa.polymerization_propensity(0.95), 0.0);
        assert_eq!(hbaa.polymerization_propensity(0.30), 0.0);
        assert_eq!(hbaa.polymerization_propensity(0.10), 0.0);
    }

    #[test]
    fn test_polymerization_propensity_trait() {
        let hbas = SickleCellModel::new(0.4);

        // Trait carriers rarely sickle - at normal oxygenation, zero
        assert_eq!(hbas.polymerization_propensity(0.95), 0.0);
        // At very low saturation, may have minimal propensity (< 5%)
        assert!(hbas.polymerization_propensity(0.30) < 0.05);
    }

    #[test]
    fn test_polymerization_propensity_disease() {
        let hbss = SickleCellModel::new(1.0);

        // Above threshold - no polymerization
        assert_eq!(hbss.polymerization_propensity(0.50), 0.0);

        // Below threshold - polymerization occurs
        let prop_low = hbss.polymerization_propensity(0.20);
        assert!(prop_low > 0.0);
    }

    #[test]
    fn test_polymerization_update() {
        let mut model = SickleCellModel::new(1.0);

        // Start oxygenated
        assert_eq!(model.polymer_fraction, 0.0);

        // Deoxygenate
        for _ in 0..100 {
            model.update_polymerization(0.20, 0.1);
        }
        assert!(model.polymer_fraction > 0.0);

        // Reoxygenate
        for _ in 0..100 {
            model.update_polymerization(0.95, 0.1);
        }
        assert!(model.polymer_fraction < 0.1);
    }

    #[test]
    fn test_stiffness_multiplier() {
        let mut model = SickleCellModel::new(1.0);

        // No polymer - normal stiffness
        assert!((model.stiffness_multiplier() - 1.0).abs() < 0.01);

        // Full polymer - max stiffness
        model.polymer_fraction = 1.0;
        assert!((model.stiffness_multiplier() - 10.0).abs() < 0.1);
    }

    #[test]
    fn test_severity() {
        let hbaa = SickleCellModel::new(0.0);
        let hbas = SickleCellModel::new(0.4);
        let hbss = SickleCellModel::new(1.0);

        assert!((hbaa.severity() - 0.0).abs() < 0.01);
        assert!(hbas.severity() > 0.0 && hbas.severity() < 0.5);
        assert!(hbss.severity() > 0.9);
    }

    #[test]
    fn test_adair_constants() {
        let hbss = SickleCellModel::new(1.0);
        let constants = hbss.effective_adair_constants();

        // Should have valid constants
        assert!(constants.k1_per_mmHg > 0.0);
        assert!(constants.k4_per_mmHg > constants.k1_per_mmHg);
    }
}
