//! A single simulated cell, owning all per-cell state.
//!
//! `Cell` is the unit of parallelism inside a [`super::World`]. It bundles
//! mesh + spectrin network + physics state + biochemistry solver +
//! metabolite pool into one self-contained object so that
//! `World::step` can iterate cells with `par_iter_mut` safely.
//!
//! The dynamics inside `Cell::step` mirror
//! [`crate::coupling::CoupledSolver::step`] — physics substeps,
//! tension accumulation, biochemistry tick, ATP-driven stiffness update.
//! Refactoring CoupledSolver to internally delegate to `Cell` is a
//! Phase 10.5 follow-on; for the moment both code paths exist and stay
//! in sync.
//!
//! Phase 11.5 marks `CoupledSolver` `#[deprecated]`; this module still
//! references `CoupledConfig` for parity, so internal uses are silenced
//! here. External consumers continue to see the deprecation warning.

#![allow(deprecated)]

use glam::Vec3;

use crate::biochemistry::{
    FullyIntegratedConfig, FullyIntegratedSolver, MetabolitePool,
};
use crate::config::Parameters;
use crate::coupling::{CoupledConfig, CoupledDiagnostics, SpectrinModulator, TensionComputer};
use crate::geometry::{Mesh, SpectrinNetwork};
use crate::physics::{PhysicsConfig, PhysicsSolver};
use crate::state::PhysicsState;

/// One simulated cell. Owns mesh, spectrin network, physics state,
/// physics solver, biochemistry solver, metabolite pool, and tension
/// history.
pub struct Cell {
    pub mesh: Mesh,
    pub spectrin: SpectrinNetwork,
    pub physics_state: PhysicsState,
    pub physics: PhysicsSolver,
    pub biochemistry: FullyIntegratedSolver,
    pub metabolites: MetabolitePool,
    pub tension_computer: TensionComputer,
    /// Last computed averaged membrane tension (pN/nm).
    pub current_tension_pN_per_nm: f64,
    /// Last computed ATP-dependent stiffness modifier.
    pub current_stiffness_modifier: f32,
    /// Simulated time (seconds) for this cell.
    pub time_sec: f64,
    /// Base shear modulus before any ATP modulation (μN/m).
    pub base_shear_modulus: f32,
    /// Base spectrin persistence length before any ATP modulation (μm).
    pub base_persistence_length: f32,
}

impl Cell {
    /// Build a fresh cell from geometry parameters and the global coupled
    /// config. The cell is positioned at the origin; multi-cell flow
    /// geometry is a Phase 12 concern.
    pub fn new(params: &Parameters, config: &CoupledConfig) -> Self {
        let mesh = Mesh::generate_rbc(&params.geometry);
        let spectrin = SpectrinNetwork::generate(&mesh, &params.geometry);

        // Reference positions are needed for any future strain-vs-rest
        // diagnostics (the existing TensionComputer derives strain from
        // the SkalakSolver's element rest configuration, so this is for
        // downstream consumers, not the per-step path).
        let mut physics_state = PhysicsState::new(mesh.vertices.len());
        let reference_positions: Vec<Vec3> = mesh
            .vertices
            .iter()
            .map(|v| v.position_vec3())
            .collect();
        physics_state.init_reference_positions(&reference_positions);

        let physics_config = PhysicsConfig {
            dt_sec: config.physics_dt_sec as f32,
            temperature_K: 310.0,
            enable_thermal_noise: true,
            membrane_damping: 5.0,
        };
        let physics = PhysicsSolver::new(&mesh, physics_config);
        let base_shear_modulus = physics.skalak_solver.material.shear_modulus_uN_per_m;
        let base_persistence_length = physics.wlc_solver.params.persistence_length_um;

        let biochem_config = FullyIntegratedConfig {
            dt_sec: config.biochem_dt_sec,
            external_glucose_mM: config.external_glucose_mM,
            po2_mmHg: config.po2_mmHg,
            oxidative_stress_multiplier: config.oxidative_stress_multiplier,
            enable_piezo1: config.enable_tension_coupling,
            membrane_tension_pN_per_nm: 0.0,
            ..Default::default()
        };
        let biochemistry = FullyIntegratedSolver::new(biochem_config);
        let metabolites = MetabolitePool::default_fully_integrated();

        Self {
            mesh,
            spectrin,
            physics_state,
            physics,
            biochemistry,
            metabolites,
            tension_computer: TensionComputer::default(),
            current_tension_pN_per_nm: 0.0,
            current_stiffness_modifier: 1.0,
            time_sec: 0.0,
            base_shear_modulus,
            base_persistence_length,
        }
    }

    /// One coupled biochemistry step
    /// (`config.physics_substeps` physics substeps + 1 biochem step).
    ///
    /// Mirrors [`crate::coupling::CoupledSolver::step`] but operates on
    /// the cell's owned state instead of external buffers.
    pub fn step(&mut self, config: &CoupledConfig, modulator: &SpectrinModulator) {
        // Phase 1: physics substeps + tension accumulation.
        for _ in 0..config.physics_substeps {
            self.physics
                .step(&mut self.mesh, &self.spectrin, &mut self.physics_state);
            if config.enable_tension_coupling {
                self.tension_computer
                    .compute_global_tension_pN_per_nm(&self.physics.skalak_solver, &self.mesh);
            }
        }

        // Phase 2: write averaged tension into biochemistry config.
        if config.enable_tension_coupling {
            self.current_tension_pN_per_nm = self.tension_computer.averaged_tension();
            self.biochemistry.config.membrane_tension_pN_per_nm = self.current_tension_pN_per_nm;
        }

        // Phase 3: biochemistry tick.
        self.biochemistry.step(&mut self.metabolites, 0.0);

        // Phase 4: update stiffness modifier from ATP.
        if config.enable_atp_stiffness_coupling {
            let atp_mM = self.metabolites.get(self.biochemistry.indices.glycolysis.atp);
            self.current_stiffness_modifier = modulator.stiffness_modifier(atp_mM);
        }

        self.time_sec += config.biochem_dt_sec;
    }

    /// Force a specific tension value (bypasses the tension computer).
    /// Used for testing without running physics.
    pub fn set_tension_override(&mut self, tension_pN_per_nm: f64) {
        self.current_tension_pN_per_nm = tension_pN_per_nm;
        self.biochemistry.config.membrane_tension_pN_per_nm = tension_pN_per_nm;
    }

    /// Reset cell state: clears tension history, resets biochem time, but
    /// does not regenerate the mesh or reset metabolite concentrations.
    pub fn reset(&mut self) {
        self.biochemistry.reset();
        self.tension_computer.reset();
        self.time_sec = 0.0;
        self.current_tension_pN_per_nm = 0.0;
        self.current_stiffness_modifier = 1.0;
    }

    /// Build a coupled-solver-style diagnostic snapshot.
    pub fn diagnostics(&self) -> CoupledDiagnostics {
        let biochem_diag = self.biochemistry.diagnostics(&self.metabolites);

        CoupledDiagnostics {
            time_sec: self.time_sec,
            membrane_tension_pN_per_nm: self.current_tension_pN_per_nm,
            stiffness_modifier: self.current_stiffness_modifier,
            effective_shear_modulus_uN_per_m:
                self.base_shear_modulus * self.current_stiffness_modifier,
            effective_persistence_length_um:
                self.base_persistence_length / self.current_stiffness_modifier,
            atp_mM: biochem_diag.atp_mM,
            adp_mM: biochem_diag.adp_mM,
            ca_cytosolic_uM: biochem_diag.ca_cytosolic_uM,
            gsh_mM: biochem_diag.gsh_mM,
            h2o2_uM: biochem_diag.h2o2_uM,
            piezo1_open_probability: piezo1_open_probability(self.current_tension_pN_per_nm),
            stiffness_status: SpectrinModulator::default()
                .status_description(biochem_diag.atp_mM),
        }
    }
}

/// Hill-style approximation of Piezo1 open probability for diagnostics.
/// Mirrors the formula in [`crate::coupling::CoupledSolver::diagnostics`].
fn piezo1_open_probability(tension_pN_per_nm: f64) -> f64 {
    if tension_pN_per_nm <= 0.0 {
        return 0.0;
    }
    let half_activation: f64 = 2.0;
    let hill: f64 = 4.0;
    let t = tension_pN_per_nm;
    let p = t.powf(hill) / (half_activation.powf(hill) + t.powf(hill));
    p.min(1.0)
}
