//! State management for the cell simulation.
//!
//! Contains data structures representing the complete state of a red blood cell,
//! including geometry, biochemistry, physics, and environment.

mod biochemistry;
mod cell;
mod environment;
mod physics;

pub use biochemistry::{BiochemistryState, HemoglobinState, IonState, MetaboliteState};
pub use cell::{CellState, GeometryState};
pub use environment::EnvironmentState;
pub use physics::PhysicsState;
