//! Malaria (Plasmodium falciparum) Infection Model.
//!
//! Simulates the metabolic effects of P. falciparum infection in RBCs:
//! - Parasite glucose competition (high glycolytic rate)
//! - Direct lactate production (bypassing host glycolysis ATP)
//! - Severe oxidative stress from hemoglobin digestion
//! - pH drop from lactate accumulation
//!
//! Parasite stages have different metabolic rates:
//! - Ring: 0.2x baseline (early, low metabolic activity)
//! - Trophozoite: 1.0x baseline (peak growth, highest metabolism)
//! - Schizont: 0.7x baseline (pre-rupture, decreased metabolism)
//!
//! Validation targets (Roth 1990, Sherman 1979):
//! | Parasitemia | ATP (mM) | Lactate (mM) | pH   |
//! |-------------|----------|--------------|------|
//! | 0%          | 2.0      | 1.5          | 7.4  |
//! | 1%          | 1.8      | 3.0          | 7.3  |
//! | 5%          | 1.5      | 8.0          | 7.1  |
//! | 10%         | 1.0      | 15.0         | 6.9  |
//!
//! References:
//! - Roth EF Jr et al. Blood. 1990;76:1151-1158
//! - Sherman IW. Microbiol Rev. 1979;43:453-495
//! - Zolg JW et al. J Biol Chem. 1984;259:8545-8551

use super::{DiseaseModel, DiseaseDiagnostics};
use crate::biochemistry::{
    FullyIntegratedConfig, FullyIntegratedSolver, FullyIntegratedIndices,
    MetabolitePool,
};

/// Parasite developmental stage.
///
/// Different stages have different metabolic rates relative to trophozoite peak.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ParasiteStage {
    /// Ring stage (0-20 hrs post-invasion): Low metabolism
    Ring,
    /// Trophozoite stage (20-36 hrs): Peak metabolism
    Trophozoite,
    /// Schizont stage (36-48 hrs): Preparing for rupture
    Schizont,
}

impl ParasiteStage {
    /// Metabolic rate multiplier relative to trophozoite (1.0)
    pub fn metabolic_multiplier(&self) -> f64 {
        match self {
            ParasiteStage::Ring => 0.2,
            ParasiteStage::Trophozoite => 1.0,
            ParasiteStage::Schizont => 0.7,
        }
    }

    /// Name for display
    pub fn name(&self) -> &'static str {
        match self {
            ParasiteStage::Ring => "Ring",
            ParasiteStage::Trophozoite => "Trophozoite",
            ParasiteStage::Schizont => "Schizont",
        }
    }
}

impl Default for ParasiteStage {
    fn default() -> Self {
        ParasiteStage::Trophozoite
    }
}

/// Configuration for malaria model.
#[derive(Debug, Clone)]
pub struct MalariaConfig {
    /// Parasite glucose consumption rate at 100% parasitemia, trophozoite stage (mM/s).
    /// Reference: Roth 1990 - parasites consume glucose ~100x faster than normal RBC
    pub glucose_consumption_mM_per_sec: f64,

    /// Parasite lactate production rate at 100% parasitemia (mM/s).
    /// Reference: Sherman 1979 - produces ~2 lactate per glucose
    pub lactate_production_mM_per_sec: f64,

    /// Oxidative stress multiplier for infected cells.
    /// Reference: Hemoglobin digestion releases free heme, generating ROS
    pub oxidative_stress_multiplier: f64,

    /// ATP steal factor (fraction of ATP consumed by parasite per parasitemia).
    /// Parasite uses host ATP for invasion and replication
    pub atp_consumption_factor: f64,

    /// Hemoglobin destruction rate per parasitemia (fraction/s).
    /// Parasites digest Hb for amino acids
    pub hb_destruction_rate_per_sec: f64,
}

impl Default for MalariaConfig {
    fn default() -> Self {
        Self {
            // Parasite glycolysis is ~100x faster than normal RBC
            // At 100% parasitemia, glucose consumption ~0.5 mM/s
            glucose_consumption_mM_per_sec: 0.5,  // Roth 1990
            // Lactate production ~2x glucose consumption (no mitochondria)
            lactate_production_mM_per_sec: 1.0,   // Sherman 1979
            // Significant oxidative stress from heme release
            oxidative_stress_multiplier: 2.0,
            // Parasite ATP consumption
            atp_consumption_factor: 0.1,  // 10% of ATP per 100% parasitemia
            // Hb destruction (eventual cell lysis)
            hb_destruction_rate_per_sec: 0.001,
        }
    }
}

/// Malaria Disease Model.
///
/// Models P. falciparum infection effects on RBC metabolism.
pub struct MalariaModel {
    /// Parasitemia fraction (0.0-1.0, typically 0.01-0.10)
    pub parasitemia_fraction: f64,
    /// Current parasite stage
    pub stage: ParasiteStage,
    /// Configuration parameters
    pub config: MalariaConfig,
    /// Cumulative hemoglobin destroyed (fraction)
    hb_destroyed_fraction: f64,
}

impl MalariaModel {
    /// Create a new malaria model with given parasitemia.
    ///
    /// # Arguments
    /// * `parasitemia_fraction` - Fraction of infected cells (0.01-0.10 typical)
    pub fn new(parasitemia_fraction: f64) -> Self {
        Self {
            parasitemia_fraction: parasitemia_fraction.clamp(0.0, 1.0),
            stage: ParasiteStage::Trophozoite,  // Default to peak metabolism
            config: MalariaConfig::default(),
            hb_destroyed_fraction: 0.0,
        }
    }

    /// Create with specific parasite stage
    pub fn with_stage(parasitemia_fraction: f64, stage: ParasiteStage) -> Self {
        Self {
            parasitemia_fraction: parasitemia_fraction.clamp(0.0, 1.0),
            stage,
            config: MalariaConfig::default(),
            hb_destroyed_fraction: 0.0,
        }
    }

    /// Create with custom configuration
    pub fn with_config(parasitemia_fraction: f64, config: MalariaConfig) -> Self {
        Self {
            parasitemia_fraction: parasitemia_fraction.clamp(0.0, 1.0),
            stage: ParasiteStage::Trophozoite,
            config,
            hb_destroyed_fraction: 0.0,
        }
    }

    /// Calculate effective parasite metabolic rate.
    ///
    /// Scales with parasitemia and developmental stage.
    pub fn effective_metabolic_rate(&self) -> f64 {
        self.parasitemia_fraction * self.stage.metabolic_multiplier()
    }

    /// Calculate disease severity.
    ///
    /// Based on parasitemia (hyperparasitemia > 5% is severe)
    pub fn severity(&self) -> f64 {
        // 0% = 0.0, 5% = 0.5, 10%+ = 1.0
        (self.parasitemia_fraction / 0.10).min(1.0)
    }

    /// Calculate effective oxidative stress.
    pub fn effective_oxidative_stress(&self) -> f64 {
        1.0 + (self.config.oxidative_stress_multiplier - 1.0) * self.parasitemia_fraction * 10.0
    }

    /// Get clinical classification based on parasitemia
    pub fn clinical_classification(&self) -> &'static str {
        if self.parasitemia_fraction < 0.01 {
            "Mild (parasitemia <1%)"
        } else if self.parasitemia_fraction < 0.05 {
            "Moderate (parasitemia 1-5%)"
        } else if self.parasitemia_fraction < 0.10 {
            "Severe (parasitemia 5-10%)"
        } else {
            "Critical (hyperparasitemia >10%)"
        }
    }
}

impl DiseaseModel for MalariaModel {
    fn name(&self) -> &'static str {
        "Malaria"
    }

    fn description(&self) -> String {
        format!(
            "P. falciparum: {:.1}% parasitemia, {} stage",
            self.parasitemia_fraction * 100.0,
            self.stage.name()
        )
    }

    fn modify_config(&self, config: &mut FullyIntegratedConfig) {
        // Increase oxidative stress
        config.oxidative_stress_multiplier *= self.effective_oxidative_stress();
    }

    fn apply_time_effects(
        &mut self,
        solver: &mut FullyIntegratedSolver,
        metabolites: &mut MetabolitePool,
        _elapsed_sec: f64,
    ) {
        // Track hemoglobin destruction
        let hb_destruction = self.config.hb_destruction_rate_per_sec
            * self.parasitemia_fraction
            * solver.config.dt_sec;
        self.hb_destroyed_fraction = (self.hb_destroyed_fraction + hb_destruction).min(1.0);

        // Update hemoglobin total if it drops significantly
        // (affects oxygen carrying capacity)
        if self.hb_destroyed_fraction > 0.1 {
            let remaining_hb = 5.0 * (1.0 - self.hb_destroyed_fraction);
            solver.hb_state.total_hb_mM = remaining_hb.max(1.0);
        }

        // Update oxidative stress based on current state
        solver.glutathione.set_oxidative_stress(self.effective_oxidative_stress());

        // Note: Direct metabolite modifications happen in modify_derivatives
        let _ = metabolites;  // Silence unused warning
    }

    fn modify_derivatives(
        &self,
        state: &[f64],
        dydt: &mut [f64],
        indices: &FullyIntegratedIndices,
    ) {
        let metabolic_rate = self.effective_metabolic_rate();
        if metabolic_rate <= 0.0 {
            return;
        }

        // Parasite glucose consumption
        // Parasites consume glucose at very high rates
        let glucose = state[indices.glycolysis.glucose];
        let glucose_consumption = self.config.glucose_consumption_mM_per_sec
            * metabolic_rate
            * glucose / (0.5 + glucose);  // Saturation kinetics
        dydt[indices.glycolysis.glucose] -= glucose_consumption;

        // Direct lactate production by parasite
        // Parasite glycolysis produces lactate but NOT through host RBC enzymes
        // (so no ATP benefit to host)
        let lactate_production = self.config.lactate_production_mM_per_sec * metabolic_rate;
        dydt[indices.glycolysis.lactate] += lactate_production;

        // Parasite ATP consumption (uses host ATP)
        let atp = state[indices.glycolysis.atp];
        let atp_consumption = self.config.atp_consumption_factor
            * metabolic_rate
            * atp / (0.5 + atp);
        dydt[indices.glycolysis.atp] -= atp_consumption;
        dydt[indices.glycolysis.adp] += atp_consumption;

        // Enhanced H2O2 production from heme release
        // Free heme from Hb digestion catalyzes ROS formation
        let h2o2_production = 0.0001 * metabolic_rate;
        dydt[indices.redox.h2o2] += h2o2_production;

        // GSH consumption from oxidative stress
        let gsh = state[indices.redox.gsh];
        let gsh_oxidation = 0.001 * metabolic_rate * gsh / (1.0 + gsh);
        dydt[indices.redox.gsh] -= 2.0 * gsh_oxidation;
        dydt[indices.redox.gssg] += gsh_oxidation;
    }

    fn diagnostics(&self, metabolites: &MetabolitePool) -> DiseaseDiagnostics {
        use crate::biochemistry::FullyIntegratedIndices;
        let indices = FullyIntegratedIndices::new();

        let mut diag = DiseaseDiagnostics::new(self.name());
        diag.severity = self.severity();

        // Current values
        let atp = metabolites.get(indices.glycolysis.atp);
        let glucose = metabolites.get(indices.glycolysis.glucose);
        let lactate = metabolites.get(indices.glycolysis.lactate);
        let gsh = metabolites.get(indices.redox.gsh);
        let gssg = metabolites.get(indices.redox.gssg);
        let h2o2 = metabolites.get(indices.redox.h2o2) * 1000.0;  // uM

        let gsh_gssg_ratio = if gssg > 1e-9 { gsh / gssg } else { f64::INFINITY };

        // Add metrics
        diag.add_metric("parasitemia_percent", self.parasitemia_fraction * 100.0);
        diag.add_metric("metabolic_rate", self.effective_metabolic_rate());
        diag.add_metric("hb_destroyed_fraction", self.hb_destroyed_fraction);
        diag.add_metric("atp_mM", atp);
        diag.add_metric("glucose_mM", glucose);
        diag.add_metric("lactate_mM", lactate);
        diag.add_metric("gsh_gssg_ratio", gsh_gssg_ratio);
        diag.add_metric("h2o2_uM", h2o2);
        diag.add_metric("oxidative_stress", self.effective_oxidative_stress());

        // Status
        diag.add_status(self.clinical_classification());
        diag.add_status(&format!("Parasite stage: {}", self.stage.name()));
        diag.add_status(&format!(
            "Metabolic activity: {:.0}% of peak",
            self.stage.metabolic_multiplier() * 100.0
        ));

        // Warnings
        if self.parasitemia_fraction > 0.05 {
            diag.add_warning("Hyperparasitemia - severe malaria risk");
        }
        if atp < 1.0 {
            diag.add_warning("ATP critically low - metabolic crisis");
        }
        if lactate > 10.0 {
            diag.add_warning("Severe lactic acidosis");
        }
        if self.hb_destroyed_fraction > 0.3 {
            diag.add_warning("Significant hemoglobin loss - anemia");
        }
        if gsh_gssg_ratio < 50.0 {
            diag.add_warning("Severe oxidative stress");
        }

        diag
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_malaria_creation() {
        let model = MalariaModel::new(0.05);
        assert!((model.parasitemia_fraction - 0.05).abs() < 0.001);
        assert_eq!(model.name(), "Malaria");
    }

    #[test]
    fn test_parasitemia_clamping() {
        let low = MalariaModel::new(-0.1);
        let high = MalariaModel::new(1.5);

        assert_eq!(low.parasitemia_fraction, 0.0);
        assert_eq!(high.parasitemia_fraction, 1.0);
    }

    #[test]
    fn test_parasite_stages() {
        let ring = ParasiteStage::Ring;
        let troph = ParasiteStage::Trophozoite;
        let schiz = ParasiteStage::Schizont;

        assert_eq!(ring.metabolic_multiplier(), 0.2);
        assert_eq!(troph.metabolic_multiplier(), 1.0);
        assert_eq!(schiz.metabolic_multiplier(), 0.7);
    }

    #[test]
    fn test_effective_metabolic_rate() {
        let model_troph = MalariaModel::with_stage(0.1, ParasiteStage::Trophozoite);
        let model_ring = MalariaModel::with_stage(0.1, ParasiteStage::Ring);

        assert!((model_troph.effective_metabolic_rate() - 0.1).abs() < 0.001);
        assert!((model_ring.effective_metabolic_rate() - 0.02).abs() < 0.001);
    }

    #[test]
    fn test_severity() {
        let mild = MalariaModel::new(0.005);  // 0.5%
        let moderate = MalariaModel::new(0.03);  // 3%
        let severe = MalariaModel::new(0.08);  // 8%
        let critical = MalariaModel::new(0.15);  // 15%

        assert!(mild.severity() < 0.1);
        assert!(moderate.severity() > 0.2 && moderate.severity() < 0.5);
        assert!(severe.severity() > 0.7 && severe.severity() < 1.0);
        assert!((critical.severity() - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_clinical_classification() {
        let mild = MalariaModel::new(0.005);
        let severe = MalariaModel::new(0.08);

        assert!(mild.clinical_classification().contains("Mild"));
        assert!(severe.clinical_classification().contains("Severe"));
    }

    #[test]
    fn test_oxidative_stress_scaling() {
        let zero = MalariaModel::new(0.0);
        let infected = MalariaModel::new(0.05);

        assert!((zero.effective_oxidative_stress() - 1.0).abs() < 0.01);
        assert!(infected.effective_oxidative_stress() > 1.0);
    }
}
