//! Multi-rate storage curve simulator (Phase 14.A).

use std::path::Path;

use anyhow::{Context as _, Result};

use crate::biochemistry::disease::StorageLesionModel;
use crate::biochemistry::{
    FullyIntegratedConfig, FullyIntegratedIndices, FullyIntegratedSolver, MetabolitePool,
};
use crate::coupling::SpectrinModulator;

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
    /// Phase 14.C — deformability index relative to day 0. 1.0 = fresh,
    /// < 1.0 = less deformable. Currently driven by `SpectrinModulator`'s
    /// ATP→stiffness coupling (Manno et al., PNAS 2002): lower ATP →
    /// higher spectrin stiffness → lower membrane deformability. This is
    /// a coarse proxy; refining to ektacytometry-style deformation
    /// indices is part of Phase 15 (disease & drug-screening
    /// extensions).
    pub deformability_relative: f64,
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
    /// Phase 14.B: solve the analytic quasi-steady-state for Na+/K+ at
    /// the start of each outer step and set the metabolite pool
    /// accordingly. This bridges the ~hour ion-equilibration timescale
    /// to the ~second bio-equilibration window, letting the curve match
    /// Hess 2010's quantitative ion targets without a ~10-minute
    /// wall-clock per simulation.
    pub force_ion_qss: bool,
}

impl Default for StorageSimConfig {
    fn default() -> Self {
        Self {
            end_day: 42.0,
            days_per_step: 1.0,
            seconds_of_bio_per_step: 1.0,
            bio_dt_sec: 1e-3,
            force_atp_dpg_targets: true,
            force_ion_qss: true,
        }
    }
}

/// Top-level orchestrator for a 42-day storage curve simulation.
pub struct StorageCurveSimulator {
    pub solver: FullyIntegratedSolver,
    pub pool: MetabolitePool,
    pub storage: StorageLesionModel,
    pub config: StorageSimConfig,
    /// ATP→stiffness coupling (Phase 8). Used to project the membrane
    /// deformability index from current ATP at each storage day.
    pub spectrin_modulator: SpectrinModulator,
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
        // Default Phase 8 coupling: ATP_ref = 2.0 mM, max stiffening
        // factor = 0.5 (50% stiffer at zero ATP). This sets the ATP →
        // deformability mapping used in `sample()`.
        let spectrin_modulator = SpectrinModulator::new(2.0, 0.5);
        Self {
            solver,
            pool,
            storage,
            config,
            spectrin_modulator,
            samples: Vec::new(),
            indices,
        }
    }

    /// Snapshot of the current state.
    fn sample(&self) -> StorageSample {
        let i = &self.indices;
        let atp = self.pool.get(i.glycolysis.atp);
        // Deformability index relative to fresh (day 0). Phase 14.C uses
        // the inverse of `SpectrinModulator`'s ATP→stiffness modifier so
        // a fresh ATP=2.0 mM gives 1.0 and ATP→0 gives 1/(1 + max_stiff).
        // This is a coarse coupling — adding GSH-driven oxidative
        // stiffening is a clean future extension once D'Alessandro et
        // al. spectrin oxidation kinetics are fit.
        let stiff = self.spectrin_modulator.stiffness_modifier(atp) as f64;
        let deformability_relative = if stiff > 1e-9 { 1.0 / stiff } else { 0.0 };

        StorageSample {
            day: self.storage.storage_days,
            atp_mM: atp,
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
            deformability_relative,
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

        // === 3.5. Phase 14.B: solve ion QSS ===
        // Sets Na/K to their quasi-steady-state values for the modified
        // pump+leak parameters. The biochemistry equilibration that
        // follows is then a small correction, not a multi-hour drive.
        if self.config.force_ion_qss {
            let atp = self.pool.get(self.indices.glycolysis.atp);
            let qss = solve_ion_qss(pump, leak, atp);
            self.pool
                .set(self.indices.ions.na_plus_cytosolic, qss.na_cyt_mM);
            self.pool
                .set(self.indices.ions.k_plus_cytosolic, qss.k_cyt_mM);
        }

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
    ///
    /// Day 0 is sampled the same way as every later day (envelope +
    /// optional QSS + equilibration) so all samples are produced by the
    /// same procedure. With `force_ion_qss = true`, day-0 K is the QSS
    /// value (~145 mM, slightly above the bare-pool default of 140 mM)
    /// — that's the K value at which pump and leak balance under
    /// physiological pump/leak parameters.
    pub fn run(&mut self) {
        let mut day = 0.0;
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
            "leak_multiplier", "oxidative_stress", "deformability_relative",
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
                format!("{:.6}", s.deformability_relative),
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

// === Ion quasi-steady-state solver (Phase 14.B) ======================
//
// At the modified pump/leak parameters of a storage day `d`, the Na+
// quasi-steady-state satisfies
//
//   3 · pump_rate(Na, ATP) = leak_Na(Na)
//
// where
//
//   pump_rate(Na, ATP) = vmax · pump_eff
//                       · (Na³ / (Km_Na³ + Na³))
//                       · (K_ext² / (Km_K² + K_ext²))
//                       · (ATP / (Km_ATP + ATP))
//   leak_Na(Na)        = g_na · leak_mult · (Na_ext - Na)
//
// Bisection on Na ∈ [0.001, Na_ext) finds the unique root. K+ at QSS
// follows from `2 · pump_rate = g_k · leak_mult · (K_cyt - K_ext)`.

/// Result of the ion QSS solver.
#[derive(Debug, Clone, Copy)]
pub struct IonQss {
    pub na_cyt_mM: f64,
    pub k_cyt_mM: f64,
}

const NA_EXT_MM: f64 = 140.0;
const K_EXT_MM: f64 = 5.0;
const KM_NA_MM: f64 = 15.0;
const KM_K_MM: f64 = 1.5;
const KM_ATP_MM: f64 = 0.2;
const VMAX_PUMP: f64 = 0.055;
const G_NA_BASE: f64 = 0.00024;
const G_K_BASE: f64 = 0.00015;

fn pump_rate(na_cyt: f64, atp_mM: f64, pump_efficiency: f64) -> f64 {
    let na3 = na_cyt * na_cyt * na_cyt;
    let na_term = na3 / (KM_NA_MM.powi(3) + na3);
    let k_term = K_EXT_MM * K_EXT_MM / (KM_K_MM * KM_K_MM + K_EXT_MM * K_EXT_MM);
    let atp_term = atp_mM / (KM_ATP_MM + atp_mM);
    VMAX_PUMP * pump_efficiency * na_term * k_term * atp_term
}

fn leak_na(na_cyt: f64, leak_multiplier: f64) -> f64 {
    G_NA_BASE * leak_multiplier * (NA_EXT_MM - na_cyt)
}

/// Solve the Na+ QSS via bisection, then derive K+ from the same pump
/// rate. Returns a sane default when ATP or pump efficiency is so low
/// that no equilibrium exists in (0, 140) — in that limit, leak
/// dominates and Na approaches Na_ext.
pub fn solve_ion_qss(
    pump_efficiency: f64,
    leak_multiplier: f64,
    atp_mM: f64,
) -> IonQss {
    // Function whose zero we seek: f(Na) = 3·pump - leak_Na.
    // f(0) is small negative (no pump, leak full inward → leak > 0
    // means f = -leak < 0). f(Na_ext) is positive (no leak, pump > 0).
    // Strictly increasing in Na on (0, Na_ext) since pump rises with
    // Na and leak falls.
    let f = |na: f64| 3.0 * pump_rate(na, atp_mM, pump_efficiency)
        - leak_na(na, leak_multiplier);

    let mut lo = 0.001_f64;
    let mut hi = NA_EXT_MM - 1e-3;
    let f_lo = f(lo);
    let f_hi = f(hi);
    let na_cyt = if f_lo >= 0.0 {
        // Pump dominates even at Na ≈ 0; QSS is at lo.
        lo
    } else if f_hi <= 0.0 {
        // Leak dominates even near Na_ext; QSS is near the wall.
        hi
    } else {
        for _ in 0..80 {
            let mid = 0.5 * (lo + hi);
            let f_mid = f(mid);
            if f_mid.abs() < 1e-10 {
                lo = mid;
                hi = mid;
                break;
            }
            if f_mid > 0.0 { hi = mid; } else { lo = mid; }
        }
        0.5 * (lo + hi)
    };

    // K+ QSS: 2·pump_rate = g_k·leak_mult·(K_cyt - K_ext).
    let pump = pump_rate(na_cyt, atp_mM, pump_efficiency);
    let g_k = G_K_BASE * leak_multiplier;
    let k_cyt = if g_k > 1e-12 {
        K_EXT_MM + 2.0 * pump / g_k
    } else {
        140.0 // unreachable in practice; sentinel
    };

    IonQss {
        na_cyt_mM: na_cyt,
        k_cyt_mM: k_cyt,
    }
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
            force_ion_qss: true,
        });
        sim.run();
        let s = sim.samples().first().expect("at least one sample");
        // Day 0 reflects default physiological values modified by ion
        // QSS — Na ≈ 10 mM and K near 140 mM (QSS slightly above bare
        // initial K=140).
        assert!((s.atp_mM - 2.0).abs() < 0.05, "ATP day-0: {}", s.atp_mM);
        assert!((s.dpg23_mM - 5.0).abs() < 0.05, "DPG day-0: {}", s.dpg23_mM);
        assert!((s.na_cyt_mM - 10.0).abs() < 1.0, "Na day-0: {}", s.na_cyt_mM);
        assert!(s.k_cyt_mM > 138.0, "K day-0: {}", s.k_cyt_mM);
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
    fn ion_qss_day_0_matches_physiological() {
        // pump_eff = 1.0, leak_mult = 1.0, ATP = 2.0 mM should land at
        // Na ~ 10 mM, K ~ 140 mM.
        let qss = solve_ion_qss(1.0, 1.0, 2.0);
        assert!(
            (qss.na_cyt_mM - 10.0).abs() < 2.0,
            "day-0 QSS Na: {} (target ~10)",
            qss.na_cyt_mM
        );
        assert!(
            (qss.k_cyt_mM - 140.0).abs() < 5.0,
            "day-0 QSS K: {} (target ~140)",
            qss.k_cyt_mM
        );
    }

    #[test]
    fn ion_qss_day_42_in_pathological_range() {
        // pump_eff = 1 - 0.02*42 = 0.16; leak_mult = 1 + 0.015*42 = 1.63;
        // ATP at day 42 ≈ 0.5. With those envelope parameters, the QSS
        // sits at Na ≈ 95 mM — pathologically high (Hess 2010 reports
        // ~60 mM). The discrepancy reflects the linear envelope
        // parameters, not a solver bug; closing it is a Phase 14.B
        // tuning step (re-fit `pump_efficiency_decay_per_day` and
        // `leak_increase_per_day` to Hess 2010 day-42 / day-14 anchors).
        let pump_eff = 1.0 - 0.02 * 42.0;
        let leak_mult = 1.0 + 0.015 * 42.0;
        let qss = solve_ion_qss(pump_eff, leak_mult, 0.5);
        println!("day-42 QSS Na={:.1} K={:.1}", qss.na_cyt_mM, qss.k_cyt_mM);
        // Sanity bounds: pathological Na (>40 mM) and reduced K. Refining
        // to Hess 2010's exact Na=60 is part of envelope re-fitting in a
        // follow-on; for now just verify the solver moves Na firmly into
        // the "membrane-failing" regime.
        assert!(
            qss.na_cyt_mM > 40.0,
            "day-42 QSS Na: {} (>40 mM expected)",
            qss.na_cyt_mM
        );
        assert!(
            qss.k_cyt_mM < 130.0,
            "day-42 QSS K: {} (<130 mM expected)",
            qss.k_cyt_mM
        );
    }

    #[test]
    fn ion_qss_monotone_in_pump_efficiency() {
        // Lower pump efficiency → higher Na (less efflux).
        let na_full = solve_ion_qss(1.0, 1.0, 2.0).na_cyt_mM;
        let na_half = solve_ion_qss(0.5, 1.0, 2.0).na_cyt_mM;
        let na_low = solve_ion_qss(0.1, 1.0, 2.0).na_cyt_mM;
        assert!(na_full < na_half);
        assert!(na_half < na_low);
    }

    #[test]
    fn sample_lookup_finds_closest_day() {
        let mut sim = StorageCurveSimulator::new(StorageSimConfig {
            end_day: 5.0,
            days_per_step: 1.0,
            seconds_of_bio_per_step: 0.05,
            bio_dt_sec: 1e-3,
            force_atp_dpg_targets: true,
            force_ion_qss: true,
        });
        sim.run();
        let s = sim.sample_at_day(3.5).expect("found");
        // Closest to 3.5 is 3.0 or 4.0 — either is fine.
        assert!((s.day - 3.0).abs() < 1e-6 || (s.day - 4.0).abs() < 1e-6);
    }
}
