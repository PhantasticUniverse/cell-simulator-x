//! Disease models for RBC pathophysiology simulation.
//!
//! This module implements various disease states that modify RBC metabolism:
//! - Storage Lesion: Time-dependent decay during blood storage
//! - Diabetic RBC: Hyperglycemia and glycation effects
//! - Malaria: Plasmodium falciparum infection
//! - Sickle Cell: HbS polymerization and mechanics
//!
//! Each disease model implements the `DiseaseModel` trait, allowing it to:
//! - Modify solver configuration at initialization
//! - Apply time-dependent effects during simulation
//! - Modify ODE derivatives for metabolite dynamics
//! - Provide disease-specific diagnostics
//!
//! References:
//! - Hess JR. Transfusion. 2010;50:2200-2214 (Storage lesion)
//! - Giugliano D et al. Diabetes Care. 1996;19:257-267 (Diabetic oxidative stress)
//! - Roth EF Jr. Blood. 1990;76:1151-1158 (Malaria metabolism)
//! - Eaton WA, Hofrichter J. Blood. 1987;70:1245-1266 (Sickle cell polymerization)

pub mod storage_lesion;
pub mod diabetic;
pub mod malaria;
pub mod sickle_cell;

pub use storage_lesion::{StorageLesionModel, StorageLesionConfig};
pub use diabetic::{DiabeticModel, DiabeticConfig};
pub use malaria::{MalariaModel, MalariaConfig, ParasiteStage};
pub use sickle_cell::{SickleCellModel, SickleCellConfig};

use super::{
    FullyIntegratedConfig, FullyIntegratedSolver, FullyIntegratedIndices,
    MetabolitePool,
};

/// Trait for disease models that modify RBC metabolism.
///
/// Disease models can hook into the simulation at multiple points:
/// 1. Configuration modification (before simulation starts)
/// 2. Time-dependent effects (each timestep)
/// 3. Derivative modification (during ODE integration)
/// 4. Diagnostics (for reporting)
pub trait DiseaseModel: Send + Sync {
    /// Name of the disease model for display
    fn name(&self) -> &'static str;

    /// Short description of the disease state
    fn description(&self) -> String;

    /// Modify the solver configuration before simulation starts.
    ///
    /// Use this to set baseline parameters like oxidative stress multipliers,
    /// external glucose levels, etc.
    fn modify_config(&self, config: &mut FullyIntegratedConfig);

    /// Apply time-dependent effects that aren't captured in the ODE derivatives.
    ///
    /// Called at each timestep after integration. Use for:
    /// - Discrete events (e.g., cell lysis threshold)
    /// - Gradual parameter changes (e.g., storage-time decay)
    /// - Direct metabolite modifications
    ///
    /// # Arguments
    /// * `solver` - Mutable reference to modify solver state
    /// * `metabolites` - Mutable reference to metabolite pool
    /// * `elapsed_sec` - Total elapsed simulation time
    fn apply_time_effects(
        &mut self,
        solver: &mut FullyIntegratedSolver,
        metabolites: &mut MetabolitePool,
        elapsed_sec: f64,
    );

    /// Modify ODE derivatives during integration.
    ///
    /// Called within the derivative computation closure. Use for:
    /// - Additional reaction fluxes (e.g., parasite metabolism)
    /// - Modified enzyme kinetics
    /// - Additional sources/sinks
    ///
    /// # Arguments
    /// * `state` - Current metabolite concentrations
    /// * `dydt` - Derivatives to modify
    /// * `indices` - Metabolite indices for accessing specific metabolites
    fn modify_derivatives(
        &self,
        state: &[f64],
        dydt: &mut [f64],
        indices: &FullyIntegratedIndices,
    );

    /// Get disease-specific diagnostic information.
    ///
    /// Returns structured data about the disease state for reporting.
    fn diagnostics(&self, metabolites: &MetabolitePool) -> DiseaseDiagnostics;
}

/// Disease-specific diagnostic information.
#[derive(Debug, Clone)]
pub struct DiseaseDiagnostics {
    /// Disease model name
    pub disease_name: String,
    /// Current disease severity (0.0 = normal, 1.0 = severe)
    pub severity: f64,
    /// Disease-specific metrics as key-value pairs
    pub metrics: Vec<(String, f64)>,
    /// Status messages
    pub status: Vec<String>,
    /// Warning messages
    pub warnings: Vec<String>,
}

impl DiseaseDiagnostics {
    /// Create a new diagnostics instance
    pub fn new(disease_name: &str) -> Self {
        Self {
            disease_name: disease_name.to_string(),
            severity: 0.0,
            metrics: Vec::new(),
            status: Vec::new(),
            warnings: Vec::new(),
        }
    }

    /// Add a metric
    pub fn add_metric(&mut self, name: &str, value: f64) {
        self.metrics.push((name.to_string(), value));
    }

    /// Add a status message
    pub fn add_status(&mut self, msg: &str) {
        self.status.push(msg.to_string());
    }

    /// Add a warning
    pub fn add_warning(&mut self, msg: &str) {
        self.warnings.push(msg.to_string());
    }

    /// Print a formatted summary
    pub fn print_summary(&self) {
        println!("=== Disease Model: {} ===", self.disease_name);
        println!("Severity: {:.1}%", self.severity * 100.0);
        println!();

        if !self.metrics.is_empty() {
            println!("Metrics:");
            for (name, value) in &self.metrics {
                println!("  {}: {:.4}", name, value);
            }
            println!();
        }

        if !self.status.is_empty() {
            println!("Status:");
            for msg in &self.status {
                println!("  {}", msg);
            }
            println!();
        }

        if !self.warnings.is_empty() {
            println!("Warnings:");
            for msg in &self.warnings {
                println!("  {}", msg);
            }
        }
    }
}

/// Registry for available disease models.
///
/// Provides factory methods for creating disease models from CLI arguments.
pub struct DiseaseRegistry;

impl DiseaseRegistry {
    /// List available disease models
    pub fn list_models() -> Vec<&'static str> {
        vec!["storage", "diabetic", "malaria", "sickle"]
    }

    /// Create a disease model by name with a parameter value.
    ///
    /// # Arguments
    /// * `name` - Disease model name ("storage", "diabetic", "malaria", "sickle")
    /// * `param` - Disease-specific parameter:
    ///   - storage: storage days (0-42)
    ///   - diabetic: external glucose (mM, 5-20)
    ///   - malaria: parasitemia fraction (0.0-0.1)
    ///   - sickle: HbS fraction (0.0-1.0)
    pub fn create(name: &str, param: f64) -> Option<Box<dyn DiseaseModel>> {
        match name.to_lowercase().as_str() {
            "storage" => Some(Box::new(StorageLesionModel::new(param))),
            "diabetic" => Some(Box::new(DiabeticModel::new(param))),
            "malaria" => Some(Box::new(MalariaModel::new(param))),
            "sickle" => Some(Box::new(SickleCellModel::new(param))),
            _ => None,
        }
    }

    /// Get help text for a disease model
    pub fn help(name: &str) -> Option<&'static str> {
        match name.to_lowercase().as_str() {
            "storage" => Some(
                "Storage Lesion Model\n\
                 Parameter: storage days (0-42)\n\
                 Effects: ATP decay, 2,3-DPG depletion, ion gradient collapse\n\
                 Reference: Hess JR. Transfusion. 2010"
            ),
            "diabetic" => Some(
                "Diabetic RBC Model\n\
                 Parameter: external glucose (mM, normal=5, diabetic=10-15)\n\
                 Effects: Oxidative stress, GSH depletion, elevated H2O2\n\
                 Reference: Giugliano D et al. Diabetes Care. 1996"
            ),
            "malaria" => Some(
                "Malaria (P. falciparum) Model\n\
                 Parameter: parasitemia fraction (0.01-0.10 = 1-10%)\n\
                 Effects: Glucose competition, lactate elevation, pH drop\n\
                 Reference: Roth EF Jr. Blood. 1990"
            ),
            "sickle" => Some(
                "Sickle Cell Model\n\
                 Parameter: HbS fraction (0=HbAA, 0.5=HbAS, 1.0=HbSS)\n\
                 Effects: P50 shift, polymerization at low pO2\n\
                 Reference: Eaton WA, Hofrichter J. Blood. 1987"
            ),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_disease_registry_list() {
        let models = DiseaseRegistry::list_models();
        assert_eq!(models.len(), 4);
        assert!(models.contains(&"storage"));
        assert!(models.contains(&"diabetic"));
        assert!(models.contains(&"malaria"));
        assert!(models.contains(&"sickle"));
    }

    #[test]
    fn test_disease_registry_create() {
        let storage = DiseaseRegistry::create("storage", 21.0);
        assert!(storage.is_some());
        assert_eq!(storage.unwrap().name(), "Storage Lesion");

        let diabetic = DiseaseRegistry::create("diabetic", 12.0);
        assert!(diabetic.is_some());
        assert_eq!(diabetic.unwrap().name(), "Diabetic RBC");

        let malaria = DiseaseRegistry::create("malaria", 0.05);
        assert!(malaria.is_some());
        assert_eq!(malaria.unwrap().name(), "Malaria");

        let sickle = DiseaseRegistry::create("sickle", 1.0);
        assert!(sickle.is_some());
        assert_eq!(sickle.unwrap().name(), "Sickle Cell");

        let unknown = DiseaseRegistry::create("unknown", 0.0);
        assert!(unknown.is_none());
    }

    #[test]
    fn test_diagnostics_builder() {
        let mut diag = DiseaseDiagnostics::new("Test Disease");
        diag.severity = 0.5;
        diag.add_metric("atp_fraction", 0.75);
        diag.add_status("Day 21 of storage");
        diag.add_warning("ATP critically low");

        assert_eq!(diag.disease_name, "Test Disease");
        assert_eq!(diag.severity, 0.5);
        assert_eq!(diag.metrics.len(), 1);
        assert_eq!(diag.status.len(), 1);
        assert_eq!(diag.warnings.len(), 1);
    }
}
