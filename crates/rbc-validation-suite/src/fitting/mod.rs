//! Parameter sensitivity, identifiability, and refit machinery.
//!
//! Phase 10's parameter-identifiability analysis lives here. Two pieces:
//!
//! - [`sensitivity`] — Latin-hypercube parameter sweep over named ranges,
//!   with response captured per validation curve.
//! - [`identifiability`] — local Fisher-information / parameter-correlation
//!   analysis to flag parameter clusters that cannot be jointly identified
//!   from available data.

pub mod identifiability;
pub mod sensitivity;

pub use identifiability::{IdentifiabilityReport, ParameterCorrelation};
pub use sensitivity::{ParameterRange, SensitivitySample, SensitivitySweep};
