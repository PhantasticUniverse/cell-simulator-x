//! Mechano-metabolic coupling between physics and biochemistry.
//!
//! This module implements bidirectional coupling between the physical mechanics
//! of the RBC membrane and its biochemistry (metabolism):
//!
//! ## Forward Coupling (Mechanics → Biochemistry)
//! - Membrane deformation → Tension → Piezo1 activation → Ca²⁺ influx → ATP release
//! - Computed via TensionComputer from Skalak strain invariants
//!
//! ## Reverse Coupling (Biochemistry → Mechanics)
//! - ATP depletion → Reduced spectrin phosphorylation → Increased stiffness
//! - Implemented via SpectrinModulator
//!
//! ## Data Flow
//! ```text
//! Physics (1μs steps)                    Biochemistry (1ms steps)
//!        │                                       │
//!        ▼                                       ▼
//! ┌──────────────┐    TensionComputer    ┌──────────────┐
//! │ Skalak Model │ ─────────────────────►│  Piezo1      │
//! │ I₁, I₂       │    tension (pN/nm)    │  Ca²⁺→ATP    │
//! └──────────────┘                       └──────────────┘
//!        ▲                                       │
//!        │     SpectrinModulator                 │
//!        └────────────────────────────────────── │
//!               stiffness modifier        ATP level
//! ```
//!
//! ## References
//! - Piezo1: Wan et al., PNAS 2008
//! - Spectrin phosphorylation: Manno et al., PNAS 2002
//! - RBC mechanics: Evans & Waugh, Biophys J 1977

pub mod coupled_solver;
pub mod spectrin_modulator;
pub mod tension_computer;

pub use coupled_solver::{CoupledConfig, CoupledDiagnostics, CoupledSolver};
pub use spectrin_modulator::SpectrinModulator;
pub use tension_computer::TensionComputer;
