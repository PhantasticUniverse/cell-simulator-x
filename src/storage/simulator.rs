//! Multi-rate storage curve simulator (Phase 14.A).

use std::path::Path;

use anyhow::{Context as _, Result};

use crate::biochemistry::disease::StorageLesionModel;
use crate::biochemistry::{
    FullyIntegratedConfig, FullyIntegratedIndices, FullyIntegratedSolver, MetabolitePool,
};

/// Per-day sample of the cell state.
#[derive(Debug, Clone, Copy)]
pub struct StorageSample {
    pub day: f64,
    pub atp_mM: f64,
    pub dpg23_mM: f64,
    pub gsh_mM: f64,
    pub gssg_mM: f64,
    pub h2o2_mM: f64,
    pub nadph_mM: f64,
    pub na_cyt_mM: f64,
    pub k_cyt_mM: f64,
    pub lactate_mM: f64,
    pub pump_efficiency: f64,
    pub leak_multiplier: f64,
    pub oxidative_stress: f64,
}

/// Configuration for the storage curve simulator.
#[derive(Debug, Clone)]
pub struct StorageSimConfig {
    /// Storage day at which the simulation ends (inclusive).
    pub end_day: f64,
    /// Days advanced per outer step. 1.0 = one sample per simulated day.
    pub days_per_step: f64,
    /// Seconds of biochemistry equilibration per outer step. Larger
    /// values give more equilibration; 1.0 s = 1000 RK4 ms-steps and
    /// reaches ~95% of the target steady state per
    /// `FullyIntegratedSolver`'s default rates.
    pub seconds_of_bio_per_step: f64,
    /// Inner ODE timestep.
    pub bio_dt_sec: f64,
    /// Whether to apply Hess-2010 ATP and 2,3-DPG envelopes directly
    /// (overrides the metabolite pool to the targets at each step). When
    /// `false`, the biochemistry's own kinetics drive the trajectory.
    pub force_atp_dpg_targets: bool,
}

impl Default for StorageSimConfig {
    fn default() -> Self {
        Self {
            end_day: 42.0,
            days_per_step: 1.0,
            seconds_of_bio_per_step: 1.0,
            bio_dt_sec: 1e-3,
            force_atp_dpg_targets: true,
        }
    }
}

/// Top-level orchestrator for a 42-day storage curve simulation.
pub struct StorageCurveSimulator {
    pub solver: FullyIntegratedSolver,
    pub pool: MetabolitePool,
    pub storage: StorageLesionModel,
    pub config: StorageSimConfig,
    samples: Vec<StorageSample>,
    indices: FullyIntegratedIndices,
}

impl StorageCurveSimulator {
    /// Build a fresh simulator at storage day 0 with default physiological
    /// metabolite levels.
    pub fn new(config: StorageSimConfig) -> Self {
        Self::with_solver_config(config, FullyIntegratedConfig::default())
    }

    /// Build a simulator with a custom `FullyIntegratedConfig` for the
    /// inner biochemistry (e.g., to set `temperature_K` to 4 °C / 277 K
    /// for refrigerated storage).
    pub fn with_solver_config(
        config: StorageSimConfig,
        mut solver_config: FullyIntegratedConfig,
    ) -> Self {
        solver_config.dt_sec = config.bio_dt_sec;
        let solver = FullyIntegratedSolver::new(solver_config);
        let pool = MetabolitePool::default_fully_integrated();
        let storage = StorageLesionModel::new(0.0);
        let indices = FullyIntegratedIndices::new();
        Self {
            solver,
            pool,
            storage,
            config,
            samples: Vec::new(),
            indices,
        }
    }

    /// Snapshot of the current state.
    fn sample(&self) -> StorageSample {
        let i = &self.indices;
        StorageSample {
            day: self.storage.storage_days,
            atp_mM: self.pool.get(i.glycolysis.atp),
            dpg23_mM: self.pool.get(i.bisphosphoglycerate_2_3),
            gsh_mM: self.pool.get(i.redox.gsh),
            gssg_mM: self.pool.get(i.redox.gssg),
            h2o2_mM: self.pool.get(i.redox.h2o2),
            nadph_mM: self.pool.get(i.redox.nadph),
            na_cyt_mM: self.pool.get(i.ions.na_plus_cytosolic),
            k_cyt_mM: self.pool.get(i.ions.k_plus_cytosolic),
            lactate_mM: self.pool.get(i.glycolysis.lactate),
            pump_efficiency: pump_efficiency(&self.storage),
            leak_multiplier: leak_multiplier(&self.storage),
            oxidative_stress: oxidative_stress(&self.storage),
        }
    }

    /// Advance to a target storage day, run the equilibration window,
    /// then push a sample.
    fn step_to_day(&mut self, target_day: f64) {
        // === 1. Update storage envelope ===
        self.storage = StorageLesionModel::new(target_day);

        // === 2. Optionally force ATP and 2,3-DPG to envelope targets ===
        if self.config.force_atp_dpg_targets {
            let target_atp = self.storage.expected_atp_mM();
            let target_dpg = self.storage.expected_dpg_mM();
            self.pool.set(self.indices.glycolysis.atp, target_atp);
            self.pool.set(self.indices.bisphosphoglycerate_2_3, target_dpg);
        }

        // === 3. Re-apply storage parameters to the solver ===
        // Oxidative stress (multiplier on H2O2 production).
        let stress = oxidative_stress(&self.storage);
        self.solver.glutathione.set_oxidative_stress(stress);

        // Pump efficiency (multiplier on Na/K-ATPase Vmax).
        let pump = pump_efficiency(&self.storage);
        self.solver.ion_homeostasis.na_k_pump.vmax_mM_per_sec = 0.055 * pump;

        // Leak conductances.
        let leak = leak_multiplier(&self.storage);
        self.solver.ion_homeostasis.config.g_na_per_sec = 0.00024 * leak;
        self.solver.ion_homeostasis.config.g_k_per_sec = 0.00015 * leak;

        // === 4. Run biochemistry equilibration ===
        let n_steps =
            (self.config.seconds_of_bio_per_step / self.config.bio_dt_sec).ceil() as usize;
        for _ in 0..n_steps {
            // No external ATP demand beyond the basal consumption — the
            // pumps are now the dominant ATP sink.
            self.solver.step(&mut self.pool, 0.0);
        }

        // === 5. Sample ===
        self.samples.push(self.sample());
    }

    /// Run the full storage curve from day 0 to `config.end_day`.
    pub fn run(&mut self) {
        // First sample at day 0 with no equilibration so the user sees
        // the actual physiological initial conditions.
        self.samples.push(self.sample());
        let mut day = self.config.days_per_step;
        // Use a small epsilon to avoid floating-point off-by-one at the
        // end day.
        while day <= self.config.end_day + 1e-9 {
            self.step_to_day(day);
            day += self.config.days_per_step;
        }
    }

    /// Borrow the recorded samples.
    pub fn samples(&self) -> &[StorageSample] {
        &self.samples
    }

    /// Find the sample whose `day` is closest to `target_day`.
    pub fn sample_at_day(&self, target_day: f64) -> Option<&StorageSample> {
        self.samples
            .iter()
            .min_by(|a, b| {
                let da = (a.day - target_day).abs();
                let db = (b.day - target_day).abs();
                da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
            })
    }

    /// Write the recorded samples as CSV. Header includes one column per
    /// `StorageSample` field.
    pub fn write_csv(&self, path: &Path) -> Result<()> {
        let mut wtr = csv::Writer::from_path(path)
            .with_context(|| format!("creating csv writer at {}", path.display()))?;
        wtr.write_record(&[
            "day", "atp_mM", "dpg23_mM", "gsh_mM", "gssg_mM", "h2o2_mM", "nadph_mM",
            "na_cyt_mM", "k_cyt_mM", "lactate_mM", "pump_efficiency",
            "leak_multiplier", "oxidative_stress",
        ])?;
        for s in &self.samples {
            wtr.write_record(&[
                format!("{:.4}", s.day),
                format!("{:.6}", s.atp_mM),
                format!("{:.6}", s.dpg23_mM),
                format!("{:.6}", s.gsh_mM),
                format!("{:.6}", s.gssg_mM),
                format!("{:.6}", s.h2o2_mM),
                format!("{:.6}", s.nadph_mM),
                format!("{:.6}", s.na_cyt_mM),
                format!("{:.6}", s.k_cyt_mM),
                format!("{:.6}", s.lactate_mM),
                format!("{:.6}", s.pump_efficiency),
                format!("{:.6}", s.leak_multiplier),
                format!("{:.6}", s.oxidative_stress),
            ])?;
        }
        wtr.flush()?;
        Ok(())
    }
}

// === Storage envelope inspection helpers =============================
//
// `StorageLesionModel`'s `pump_efficiency`, `leak_multiplier`, and
// `oxidative_stress` fields are private. Re-derive them from the
// public-config-driven formula so the simulator can read them without
// touching the model's internals.

fn pump_efficiency(s: &StorageLesionModel) -> f64 {
    (1.0 - s.config.pump_efficiency_decay_per_day * s.storage_days).max(0.1)
}

fn leak_multiplier(s: &StorageLesionModel) -> f64 {
    1.0 + s.config.leak_increase_per_day * s.storage_days
}

fn oxidative_stress(s: &StorageLesionModel) -> f64 {
    1.0 + s.config.oxidative_stress_increase_per_day * s.storage_days
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn day_0_sample_matches_physiological_defaults() {
        let mut sim = StorageCurveSimulator::new(StorageSimConfig {
            end_day: 0.0,
            days_per_step: 1.0,
            seconds_of_bio_per_step: 0.0,
            bio_dt_sec: 1e-3,
            force_atp_dpg_targets: true,
        });
        sim.run();
        let s = sim.samples().first().expect("at least one sample");
        // Day 0 should reflect default physiological values.
        assert!((s.atp_mM - 2.0).abs() < 0.05, "ATP day-0: {}", s.atp_mM);
        assert!((s.dpg23_mM - 5.0).abs() < 0.05, "DPG day-0: {}", s.dpg23_mM);
        assert!((s.na_cyt_mM - 10.0).abs() < 0.5, "Na day-0: {}", s.na_cyt_mM);
        assert!((s.k_cyt_mM - 140.0).abs() < 1.0, "K day-0: {}", s.k_cyt_mM);
    }

    #[test]
    fn pump_efficiency_decays_monotonically() {
        let s0 = StorageLesionModel::new(0.0);
        let s14 = StorageLesionModel::new(14.0);
        let s42 = StorageLesionModel::new(42.0);
        assert!((pump_efficiency(&s0) - 1.0).abs() < 1e-6);
        assert!(pump_efficiency(&s14) < pump_efficiency(&s0));
        assert!(pump_efficiency(&s42) < pump_efficiency(&s14));
        assert!(pump_efficiency(&s42) >= 0.1, "pump efficiency floor");
    }

    #[test]
    fn sample_lookup_finds_closest_day() {
        let mut sim = StorageCurveSimulator::new(StorageSimConfig {
            end_day: 5.0,
            days_per_step: 1.0,
            seconds_of_bio_per_step: 0.05,
            bio_dt_sec: 1e-3,
            force_atp_dpg_targets: true,
        });
        sim.run();
        let s = sim.sample_at_day(3.5).expect("found");
        // Closest to 3.5 is 3.0 or 4.0 — either is fine.
        assert!((s.day - 3.0).abs() < 1e-6 || (s.day - 4.0).abs() < 1e-6);
    }
}
