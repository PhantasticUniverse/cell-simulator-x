//! External fluid for single-cell-in-vessel scenarios (Phase 12).
//!
//! Phase 12.A (this module): analytic Poiseuille flow in a cylindrical
//! channel + per-vertex drag coupling. No actual Navier-Stokes solver
//! yet — the channel geometry and the velocity field are imposed
//! analytically. This lets us drive a cell with a known shear rate /
//! velocity profile and observe the membrane deformation response.
//!
//! Phase 12.B / 12.C (deferred): wiring drag into `PhysicsBackend.step()`
//! and quantitative validation (parachute-shape under Poiseuille,
//! tank-treading frequency vs shear rate per Fischer 2007).
//!
//! ## References
//! - Poiseuille flow: any fluid mechanics textbook; we use the
//!   axisymmetric form `v_z(r) = v_max (1 - (r/R)²)`.
//! - Drag coefficient: Stokes drag on a sphere `F = 6π μ a (v - u)`,
//!   adapted per vertex with an effective hydrodynamic radius.
//! - Fischer 2007: tank-treading frequency benchmark
//!   (Fischer, T. M. _Biophys. J._ 93, 2553–2561).

pub mod poiseuille;

pub use poiseuille::{
    apply_drag_to_external_forces, drag_force_uN, CylindricalChannel, Poiseuille,
};
