//! # Cell Simulator X
//!
//! GPU-accelerated human red blood cell simulation engine integrating
//! mechanics, metabolism, and oxygen transport.
//!
//! ## Overview
//!
//! Cell Simulator X provides a comprehensive computational model of the human
//! red blood cell (RBC), simulating the complex interplay between:
//!
//! - **Membrane Mechanics**: Spectrin network elasticity, Skalak membrane model
//! - **Core Metabolism**: Glycolysis, Rapoport-Luebering shunt (2,3-DPG)
//! - **Oxygen Transport**: Adair 4-site hemoglobin model with allosteric effects
//! - **Redox Homeostasis**: Pentose phosphate pathway, glutathione cycle
//! - **Ion Homeostasis**: Na+/K+-ATPase, Ca2+ regulation via Piezo1
//! - **Disease Models**: Storage lesion, diabetic, malaria, sickle cell
//!
//! ## Quick Start
//!
//! ```rust,ignore
//! use cell_simulator_x::{
//!     FullyIntegratedSolver, FullyIntegratedConfig, MetabolitePool
//! };
//!
//! // Create solver with default configuration
//! let config = FullyIntegratedConfig::default();
//! let mut solver = FullyIntegratedSolver::new(config);
//!
//! // Initialize metabolite pool with physiological concentrations
//! let mut metabolites = MetabolitePool::default_fully_integrated();
//!
//! // Run 60 seconds of simulation
//! let atp_consumption = 0.001; // mM/s external ATP demand
//! solver.run(&mut metabolites, 60.0, atp_consumption);
//!
//! // Get diagnostics
//! let diag = solver.diagnostics(&metabolites);
//! println!("ATP: {:.2} mM", diag.atp_mM);
//! println!("GSH/GSSG: {:.0}", diag.gsh_gssg_ratio);
//! ```
//!
//! ## Module Structure
//!
//! - [`biochemistry`] - Metabolic pathways: glycolysis, PPP, glutathione, hemoglobin
//! - [`coupling`] - Mechano-metabolic coupling: tension ↔ metabolism feedback
//! - [`geometry`] - RBC morphology: Fung-Tong shape, spectrin network
//! - [`physics`] - Membrane mechanics: DPD, WLC, Skalak model
//! - [`render`] - WebGPU/Metal visualization
//! - [`config`] - Parameter loading from JSON
//! - [`state`] - Cell state management
//!
//! ## Key Types
//!
//! | Type | Purpose |
//! |------|---------|
//! | [`FullyIntegratedSolver`] | Main solver integrating all 38 metabolites |
//! | [`CoupledSolver`] | Mechano-metabolic coupling (physics + biochemistry) |
//! | [`MetabolitePool`] | Metabolite concentration storage |
//! | [`HemoglobinSolver`] | Oxygen binding with allosteric effects |
//! | [`PhysicsSolver`] | Membrane mechanics (DPD, WLC, Skalak) |
//!
//! ## Disease Models
//!
//! Pathophysiological models accessible via the [`DiseaseRegistry`]:
//!
//! - [`StorageLesionModel`] - Blood storage aging (ATP decay, 2,3-DPG depletion)
//! - [`DiabeticModel`] - Hyperglycemia effects (elevated glucose, oxidative stress)
//! - [`MalariaModel`] - P. falciparum infection (glucose competition, oxidative stress)
//! - [`SickleCellModel`] - HbS polymerization (P50 shift, chronic oxidative stress)
//!
//! ## References
//!
//! - Beutler E (1984) *Red Cell Metabolism: A Manual of Biochemical Methods*
//! - Imai K (1982) *Allosteric Effects in Haemoglobin*
//! - Evans EA, Hochmuth RM (1977) Membrane viscoelasticity. *Biophys J* 16:1

// Allow non-snake-case for unit suffixes in field names (mM, mmHg, μN, etc.)
// This follows the project convention of including units in names.
#![allow(non_snake_case)]

pub mod biochemistry;
pub mod config;
pub mod coupling;
pub mod export;
pub mod geometry;
pub mod physics;
pub mod render;
pub mod state;

pub use biochemistry::{
    MetabolismSolver, MetabolismConfig, MetabolitePool, MetabolismDiagnostics,
    ExtendedMetaboliteIndices, MetaboliteIndices,
    HemoglobinSolver, HemoglobinState, OxygenDiagnostics, AdairConstants,
    STANDARD_PH, STANDARD_TEMPERATURE_K, STANDARD_DPG_MM, STANDARD_PCO2_MMHG,
    // Phase 5: Metabolism-Oxygen Integration
    IntegratedSolver, IntegratedDiagnostics, IntegratedEnvironment, PhBufferModel,
    run_integrated_diagnostics,
    // Phase 6: Redox Metabolism (PPP, Glutathione, Piezo1)
    RedoxSolver, RedoxConfig, RedoxDiagnostics, RedoxIndices,
    PentosePhosphatePathway, GlutathioneCycle, Piezo1System, Piezo1Diagnostics,
    initialize_redox_metabolites,
    // Phase 6b: Fully Integrated Solver (Glycolysis + PPP + Glutathione + O2 + Ions)
    FullyIntegratedSolver, FullyIntegratedConfig, FullyIntegratedDiagnostics,
    FullyIntegratedIndices, run_full_integration_diagnostics,
    // Phase 6b: Ion Homeostasis
    IonIndices, IonHomeostasisSystem, IonHomeostasisConfig, IonDiagnostics, NaKATPase,
    initialize_ion_metabolites,
    // Phase 7: Disease Models
    DiseaseModel, DiseaseDiagnostics, DiseaseRegistry,
    StorageLesionModel, StorageLesionConfig,
    DiabeticModel, DiabeticConfig,
    MalariaModel, MalariaConfig, ParasiteStage,
    SickleCellModel, SickleCellConfig,
};
pub use config::Parameters;
// Phase 8: Mechano-Metabolic Coupling
pub use coupling::{CoupledSolver, CoupledConfig, CoupledDiagnostics, TensionComputer, SpectrinModulator};
// Export module
pub use export::{CsvExporter, TimeSeriesRecord, export_state_json, save_screenshot};
pub use geometry::{Mesh, SpectrinNetwork};
pub use physics::{PhysicsConfig, PhysicsSolver};
pub use render::{Camera, ExportAction, HudColors, HudOverlay, HudState, HudTheme, RenderState};
pub use state::{CellState, DiseaseIndicator, MetaboliteStatus, SimulationMetrics, SimulationMode};
