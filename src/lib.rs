//! Cell Simulator X - GPU-accelerated human red blood cell simulation engine
//!
//! This library integrates mechanics, metabolism, and oxygen transport
//! to provide a comprehensive RBC simulation.

// Allow non-snake-case for unit suffixes in field names (mM, mmHg, μN, etc.)
// This follows the project convention of including units in names.
#![allow(non_snake_case)]

pub mod biochemistry;
pub mod config;
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
pub use geometry::{Mesh, SpectrinNetwork};
pub use physics::{PhysicsConfig, PhysicsSolver};
pub use render::{Camera, RenderState};
pub use state::CellState;
