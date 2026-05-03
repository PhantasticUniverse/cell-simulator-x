//! Storage lesion at physiological timescale (Phase 14).
//!
//! This is the headline scientific deliverable: a 42-day blood-bank
//! storage simulation that produces a metabolite trajectory matching
//! published storage metabolomics (Hess 2010, D'Alessandro et al.).
//!
//! ## Multi-rate scheme
//!
//! Three timescales run together:
//!
//! | Loop | Cadence | Subsystem |
//! |------|---------|-----------|
//! | Inner | 1 ms | `FullyIntegratedSolver` (38-species ODE) |
//! | Middle | seconds | Biochemistry equilibration window per storage-day step |
//! | Outer | days | `StorageLesionModel` advances; pump efficiency / leak / oxidative stress recomputed |
//!
//! At each "storage-day step" the simulator (a) overrides ATP and
//! 2,3-DPG to the empirical Hess-2010 targets for that day, (b)
//! refreshes the Na⁺/K⁺-ATPase, leak conductances, and oxidative-stress
//! multiplier per `StorageLesionModel`, then (c) runs biochemistry for
//! `seconds_of_bio_per_step` to let the rest of the metabolite pool
//! equilibrate to the new storage state.
//!
//! ## Why analytic envelopes for the slowest processes
//!
//! Real storage lesion is driven by minute-to-hour processes (membrane
//! lipid peroxidation, gradual Hb autoxidation, ATPase glycation) that
//! are not in the 38-species kinetic model. The empirical envelope
//! (`StorageLesionModel`) absorbs them — ATP follows
//! `2 mM × exp(-ln 2 × t / 21 d)` per Hess 2010, 2,3-DPG drops linearly,
//! pump activity decays exponentially. The biochemistry's job at each
//! storage day is to find the *consistent* steady state of the rest of
//! the pool given those forced metabolites and modified parameters.
//!
//! ## Validation
//!
//! `tests/storage_curve.rs` verifies the curve at days 0, 14, 21, 42
//! matches Hess 2010's reported ranges:
//!
//! | Day | ATP (mM) | 2,3-DPG (mM) | Na+ (mM) | K+ (mM) |
//! |-----|----------|--------------|----------|---------|
//! | 0   | ≈ 2.0    | ≈ 5.0        | ≈ 10     | ≈ 140   |
//! | 14  | ≈ 1.5    | ≈ 0.5        | 20–30    | 110–130 |
//! | 42  | ≈ 0.5    | 0            | 50–70    | 80–100  |

pub mod additive;
pub mod sensitivity;
pub mod simulator;

pub use additive::AdditiveSolution;
pub use sensitivity::{
    run_oat_sensitivity, top_sensitivities_by_deformability, write_csv as write_sensitivity_csv,
    SensitivityParameter, SensitivityRow, ALL_PARAMETERS,
};
pub use simulator::{StorageCurveSimulator, StorageSample, StorageSimConfig};
