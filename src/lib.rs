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
    ExtendedMetaboliteIndices,
    HemoglobinSolver, HemoglobinState, OxygenDiagnostics,
    STANDARD_PH, STANDARD_TEMPERATURE_K, STANDARD_DPG_MM, STANDARD_PCO2_MMHG,
};
pub use config::Parameters;
pub use geometry::{Mesh, SpectrinNetwork};
pub use physics::{PhysicsConfig, PhysicsSolver};
pub use render::{Camera, RenderState};
pub use state::CellState;
