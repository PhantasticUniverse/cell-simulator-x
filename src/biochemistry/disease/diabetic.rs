//! Diabetic RBC Model for hyperglycemia effects.
//!
//! Simulates the effects of elevated blood glucose on RBC metabolism:
//! - Elevated external glucose (10-15 mM vs normal 5 mM)
//! - Increased oxidative stress
//! - Enhanced NADPH consumption
//! - HbA1c formation (glycation tracking)
//!
//! Validation targets (Giugliano 1996):
//! | Glucose (mM) | GSH/GSSG | H2O2 (uM) | NADPH/NADP+ |
//! |--------------|----------|-----------|-------------|
//! | 5 (normal)   | >200     | <5        | 10-20       |
//! | 10           | 100-200  | 5-10      | 8-12        |
//! | 15           | <100     | >10       | 5-10        |
//!
//! References:
//! - Giugliano D et al. Diabetes Care. 1996;19:257-267
//! - Bunn HF et al. J Clin Invest. 1981;67:1361-1369
//! - Wolff SP. Br Med Bull. 1993;49:642-652

use super::{DiseaseModel, DiseaseDiagnostics};
use crate::biochemistry::{
    FullyIntegratedConfig, FullyIntegratedSolver, FullyIntegratedIndices,
    MetabolitePool,
};

/// Configuration for diabetic RBC model.
#[derive(Debug, Clone)]
pub struct DiabeticConfig {
    /// Oxidative stress multiplier for diabetic state.
    /// Reference: Giugliano 1996 - ~1.5-2x increase in oxidative markers
    pub oxidative_stress_multiplier: f64,

    /// Glycation rate per mM glucose per day (fraction of Hb glycated).
    /// Reference: Bunn 1981 - ~0.01% per day per mM excess glucose
    pub glycation_rate_per_mM_glucose_per_day: f64,

    /// Additional NADPH consumption rate (mM/s per mM excess glucose).
    /// Represents increased demand from polyol pathway and non-enzymatic glycation
    pub nadph_consumption_per_mM_excess_glucose: f64,

    /// GSH oxidation rate multiplier at elevated glucose.
    /// Reference: Wolff 1993 - glucose autoxidation generates ROS
    pub gsh_oxidation_multiplier: f64,

    /// Normal (euglycemic) glucose concentration (mM)
    pub normal_glucose_mM: f64,
}

impl Default for DiabeticConfig {
    fn default() -> Self {
        Self {
            oxidative_stress_multiplier: 1.5,  // Giugliano 1996
            glycation_rate_per_mM_glucose_per_day: 0.0001,  // Bunn 1981
            nadph_consumption_per_mM_excess_glucose: 0.0001,  // Estimated
            gsh_oxidation_multiplier: 1.3,  // Moderate increase
            normal_glucose_mM: 5.0,
        }
    }
}

/// Diabetic RBC Disease Model.
///
/// Models the effects of chronic hyperglycemia on RBC redox metabolism.
pub struct DiabeticModel {
    /// External glucose concentration (mM)
    pub external_glucose_mM: f64,
    /// Duration of diabetic exposure (days) - for HbA1c tracking
    pub duration_days: f64,
    /// Configuration parameters
    pub config: DiabeticConfig,
    /// Estimated HbA1c fraction (0-0.2 typical range)
    hba1c_fraction: f64,
}

impl DiabeticModel {
    /// Create a new diabetic model with given external glucose.
    ///
    /// # Arguments
    /// * `external_glucose_mM` - External glucose concentration (5-20 mM)
    pub fn new(external_glucose_mM: f64) -> Self {
        let config = DiabeticConfig::default();
        Self {
            external_glucose_mM: external_glucose_mM.max(config.normal_glucose_mM),
            duration_days: 0.0,
            config,
            hba1c_fraction: 0.05,  // Normal ~5%
        }
    }

    /// Create with explicit duration for HbA1c estimation
    pub fn with_duration(external_glucose_mM: f64, duration_days: f64) -> Self {
        let config = DiabeticConfig::default();
        let mut model = Self {
            external_glucose_mM: external_glucose_mM.max(config.normal_glucose_mM),
            duration_days: duration_days.max(0.0),
            config,
            hba1c_fraction: 0.05,
        };
        model.update_hba1c();
        model
    }

    /// Create with custom configuration
    pub fn with_config(external_glucose_mM: f64, config: DiabeticConfig) -> Self {
        Self {
            external_glucose_mM: external_glucose_mM.max(config.normal_glucose_mM),
            duration_days: 0.0,
            config,
            hba1c_fraction: 0.05,
        }
    }

    /// Update HbA1c estimate based on glucose and duration
    fn update_hba1c(&mut self) {
        let excess_glucose = (self.external_glucose_mM - self.config.normal_glucose_mM).max(0.0);
        // HbA1c accumulates based on excess glucose and time
        // Reference: Bunn 1981 - linear relationship with average glucose
        self.hba1c_fraction = 0.05 +
            self.config.glycation_rate_per_mM_glucose_per_day * excess_glucose * self.duration_days;
        self.hba1c_fraction = self.hba1c_fraction.min(0.20);  // Cap at 20% (extreme)
    }

    /// Calculate excess glucose above normal
    pub fn excess_glucose_mM(&self) -> f64 {
        (self.external_glucose_mM - self.config.normal_glucose_mM).max(0.0)
    }

    /// Calculate disease severity based on glucose level
    /// 0 = normal (5 mM), 1 = severe diabetic (20 mM)
    pub fn severity(&self) -> f64 {
        ((self.external_glucose_mM - self.config.normal_glucose_mM) / 15.0).clamp(0.0, 1.0)
    }

    /// Calculate effective oxidative stress multiplier.
    ///
    /// Scales with glucose level above normal.
    pub fn effective_oxidative_stress(&self) -> f64 {
        // Linear interpolation: 1.0 at 5mM, config.oxidative_stress_multiplier at 15mM
        let excess_factor = self.excess_glucose_mM() / 10.0;
        1.0 + (self.config.oxidative_stress_multiplier - 1.0) * excess_factor.min(1.0)
    }

    /// Estimate HbA1c percentage
    pub fn hba1c_percent(&self) -> f64 {
        self.hba1c_fraction * 100.0
    }
}

impl DiseaseModel for DiabeticModel {
    fn name(&self) -> &'static str {
        "Diabetic RBC"
    }

    fn description(&self) -> String {
        format!(
            "Diabetic: glucose={:.1}mM, HbA1c={:.1}%",
            self.external_glucose_mM,
            self.hba1c_percent()
        )
    }

    fn modify_config(&self, config: &mut FullyIntegratedConfig) {
        // Set elevated external glucose
        config.external_glucose_mM = self.external_glucose_mM;

        // Increase oxidative stress
        config.oxidative_stress_multiplier *= self.effective_oxidative_stress();
    }

    fn apply_time_effects(
        &mut self,
        solver: &mut FullyIntegratedSolver,
        _metabolites: &mut MetabolitePool,
        elapsed_sec: f64,
    ) {
        // Update duration (convert seconds to days for HbA1c tracking)
        // Note: In reality HbA1c forms over weeks, but we use accelerated kinetics
        self.duration_days += elapsed_sec / 86400.0;  // seconds to days
        self.update_hba1c();

        // Update oxidative stress in glutathione cycle
        solver.glutathione.set_oxidative_stress(self.effective_oxidative_stress());
    }

    fn modify_derivatives(
        &self,
        state: &[f64],
        dydt: &mut [f64],
        indices: &FullyIntegratedIndices,
    ) {
        let excess_glucose = self.excess_glucose_mM();
        if excess_glucose <= 0.0 {
            return;
        }

        // Additional NADPH consumption from polyol pathway and glycation stress
        // The polyol pathway (aldose reductase) consumes NADPH when glucose is high
        let nadph = state[indices.redox.nadph];
        let extra_nadph_consumption = self.config.nadph_consumption_per_mM_excess_glucose
            * excess_glucose
            * nadph / (0.1 + nadph);
        dydt[indices.redox.nadph] -= extra_nadph_consumption;
        dydt[indices.redox.nadp_plus] += extra_nadph_consumption;

        // Enhanced GSH oxidation from glucose autoxidation
        // Reference: Wolff 1993 - glucose generates ROS non-enzymatically
        let gsh = state[indices.redox.gsh];
        let base_oxidation = 0.0005 * gsh;  // Baseline rate
        let glucose_induced = base_oxidation * (self.config.gsh_oxidation_multiplier - 1.0)
            * (excess_glucose / 10.0);
        dydt[indices.redox.gsh] -= 2.0 * glucose_induced;
        dydt[indices.redox.gssg] += glucose_induced;

        // Slight increase in H2O2 production from glucose autoxidation
        let h2o2_production = 0.00001 * excess_glucose;
        dydt[indices.redox.h2o2] += h2o2_production;
    }

    fn diagnostics(&self, metabolites: &MetabolitePool) -> DiseaseDiagnostics {
        use crate::biochemistry::FullyIntegratedIndices;
        let indices = FullyIntegratedIndices::new();

        let mut diag = DiseaseDiagnostics::new(self.name());
        diag.severity = self.severity();

        // Current values
        let glucose = metabolites.get(indices.glycolysis.glucose);
        let nadph = metabolites.get(indices.redox.nadph);
        let nadp = metabolites.get(indices.redox.nadp_plus);
        let gsh = metabolites.get(indices.redox.gsh);
        let gssg = metabolites.get(indices.redox.gssg);
        let h2o2 = metabolites.get(indices.redox.h2o2) * 1000.0;  // Convert to uM

        let nadph_nadp_ratio = if nadp > 1e-9 { nadph / nadp } else { f64::INFINITY };
        let gsh_gssg_ratio = if gssg > 1e-9 { gsh / gssg } else { f64::INFINITY };

        // Add metrics
        diag.add_metric("external_glucose_mM", self.external_glucose_mM);
        diag.add_metric("internal_glucose_mM", glucose);
        diag.add_metric("excess_glucose_mM", self.excess_glucose_mM());
        diag.add_metric("hba1c_percent", self.hba1c_percent());
        diag.add_metric("oxidative_stress_multiplier", self.effective_oxidative_stress());
        diag.add_metric("nadph_nadp_ratio", nadph_nadp_ratio);
        diag.add_metric("gsh_gssg_ratio", gsh_gssg_ratio);
        diag.add_metric("h2o2_uM", h2o2);

        // Status
        let glycemic_status = if self.external_glucose_mM < 7.0 {
            "Euglycemic"
        } else if self.external_glucose_mM < 11.0 {
            "Mildly hyperglycemic"
        } else if self.external_glucose_mM < 15.0 {
            "Moderately hyperglycemic"
        } else {
            "Severely hyperglycemic"
        };
        diag.add_status(glycemic_status);
        diag.add_status(&format!("Estimated HbA1c: {:.1}%", self.hba1c_percent()));

        // Warnings based on diabetic targets
        if nadph_nadp_ratio < 5.0 {
            diag.add_warning("NADPH/NADP+ ratio low - antioxidant capacity compromised");
        }
        if gsh_gssg_ratio < 100.0 {
            diag.add_warning("GSH/GSSG ratio low - oxidative stress elevated");
        }
        if h2o2 > 10.0 {
            diag.add_warning("H2O2 elevated - oxidative damage risk");
        }
        if self.hba1c_percent() > 8.0 {
            diag.add_warning("HbA1c elevated - poor glycemic control");
        }

        diag
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_diabetic_creation() {
        let model = DiabeticModel::new(12.0);
        assert_eq!(model.external_glucose_mM, 12.0);
        assert_eq!(model.name(), "Diabetic RBC");
    }

    #[test]
    fn test_excess_glucose() {
        let normal = DiabeticModel::new(5.0);
        let diabetic = DiabeticModel::new(15.0);

        assert!((normal.excess_glucose_mM() - 0.0).abs() < 0.01);
        assert!((diabetic.excess_glucose_mM() - 10.0).abs() < 0.01);
    }

    #[test]
    fn test_severity() {
        let normal = DiabeticModel::new(5.0);
        let mild = DiabeticModel::new(10.0);
        let severe = DiabeticModel::new(20.0);

        assert!((normal.severity() - 0.0).abs() < 0.01);
        assert!(mild.severity() > 0.0 && mild.severity() < 1.0);
        assert!((severe.severity() - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_oxidative_stress_scaling() {
        let normal = DiabeticModel::new(5.0);
        let diabetic = DiabeticModel::new(15.0);

        assert!((normal.effective_oxidative_stress() - 1.0).abs() < 0.01);
        assert!(diabetic.effective_oxidative_stress() > 1.0);
        assert!(diabetic.effective_oxidative_stress() <= 1.5);
    }

    #[test]
    fn test_minimum_glucose() {
        // Should clamp to normal glucose level
        let model = DiabeticModel::new(2.0);
        assert_eq!(model.external_glucose_mM, 5.0);
    }

    #[test]
    fn test_hba1c_estimation() {
        let model = DiabeticModel::with_duration(15.0, 60.0);  // 15mM for 60 days
        assert!(model.hba1c_percent() > 5.0);
        assert!(model.hba1c_percent() < 20.0);
    }
}
