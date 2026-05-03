//! Phase 14.D.3: storage envelope sensitivity analysis.
//!
//! One-at-a-time (OAT) sensitivity sweep over named envelope parameters.
//! For each parameter, ±20% perturbations are applied around the
//! baseline; the simulator re-runs and the resulting fractional change
//! in day-42 deformability and ATP is recorded.
//!
//! The output is a table that surfaces which biological levers most
//! strongly drive the day-42 storage outcome — useful for clinical
//! interpretation (which preservation strategy is highest leverage?)
//! and for parameter-fitting prioritization (which parameter most needs
//! tighter literature anchoring?).

use std::path::Path;

use anyhow::{Context as _, Result};

use crate::biochemistry::disease::StorageLesionConfig;
use crate::storage::simulator::{StorageCurveSimulator, StorageSimConfig};

/// Named storage envelope parameter selectable for sensitivity analysis.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SensitivityParameter {
    AtpHalfLifeDays,
    PumpEfficiencyDecayPerDay,
    LeakIncreasePerDay,
    DpgLossRateMmPerDay,
    OxidativeStressIncreasePerDay,
}

impl SensitivityParameter {
    pub fn name(self) -> &'static str {
        match self {
            Self::AtpHalfLifeDays => "atp_half_life_days",
            Self::PumpEfficiencyDecayPerDay => "pump_efficiency_decay_per_day",
            Self::LeakIncreasePerDay => "leak_increase_per_day",
            Self::DpgLossRateMmPerDay => "dpg_loss_rate_mM_per_day",
            Self::OxidativeStressIncreasePerDay => "oxidative_stress_increase_per_day",
        }
    }

    pub fn baseline_value(self, cfg: &StorageLesionConfig) -> f64 {
        match self {
            Self::AtpHalfLifeDays => cfg.atp_decay_half_life_days,
            Self::PumpEfficiencyDecayPerDay => cfg.pump_efficiency_decay_per_day,
            Self::LeakIncreasePerDay => cfg.leak_increase_per_day,
            Self::DpgLossRateMmPerDay => cfg.dpg_loss_rate_mM_per_day,
            Self::OxidativeStressIncreasePerDay => cfg.oxidative_stress_increase_per_day,
        }
    }

    pub fn set(self, cfg: &mut StorageLesionConfig, value: f64) {
        match self {
            Self::AtpHalfLifeDays => cfg.atp_decay_half_life_days = value,
            Self::PumpEfficiencyDecayPerDay => cfg.pump_efficiency_decay_per_day = value,
            Self::LeakIncreasePerDay => cfg.leak_increase_per_day = value,
            Self::DpgLossRateMmPerDay => cfg.dpg_loss_rate_mM_per_day = value,
            Self::OxidativeStressIncreasePerDay => {
                cfg.oxidative_stress_increase_per_day = value
            }
        }
    }
}

/// One row of sensitivity output: one parameter perturbed by ±20%.
#[derive(Debug, Clone, Copy)]
pub struct SensitivityRow {
    pub parameter: SensitivityParameter,
    pub baseline_value: f64,
    pub perturbed_value: f64,
    pub perturbation_pct: f64,
    pub day42_atp_baseline: f64,
    pub day42_atp_perturbed: f64,
    pub day42_atp_relative_change: f64,
    pub day42_deformability_baseline: f64,
    pub day42_deformability_perturbed: f64,
    pub day42_deformability_relative_change: f64,
}

/// Run a one-at-a-time sensitivity sweep over the supplied parameters
/// at ±perturbation_fraction (e.g. 0.20 for ±20%) around each parameter's
/// baseline value in `base_config.lesion_config`.
///
/// Returns 2 rows per parameter (one positive, one negative perturbation).
pub fn run_oat_sensitivity(
    base_config: &StorageSimConfig,
    parameters: &[SensitivityParameter],
    perturbation_fraction: f64,
) -> Vec<SensitivityRow> {
    let baseline = run_and_extract_day42(base_config);

    let mut rows = Vec::with_capacity(parameters.len() * 2);
    for &param in parameters {
        for &sign in &[-1.0_f64, 1.0_f64] {
            let factor = 1.0 + sign * perturbation_fraction;
            let mut perturbed_lesion = base_config.lesion_config.clone();
            let baseline_val = param.baseline_value(&base_config.lesion_config);
            let perturbed_val = baseline_val * factor;
            param.set(&mut perturbed_lesion, perturbed_val);

            let perturbed_config = StorageSimConfig {
                lesion_config: perturbed_lesion,
                ..base_config.clone()
            };
            let perturbed = run_and_extract_day42(&perturbed_config);

            rows.push(SensitivityRow {
                parameter: param,
                baseline_value: baseline_val,
                perturbed_value: perturbed_val,
                perturbation_pct: sign * perturbation_fraction * 100.0,
                day42_atp_baseline: baseline.atp_mM,
                day42_atp_perturbed: perturbed.atp_mM,
                day42_atp_relative_change: rel_change(baseline.atp_mM, perturbed.atp_mM),
                day42_deformability_baseline: baseline.deformability_relative,
                day42_deformability_perturbed: perturbed.deformability_relative,
                day42_deformability_relative_change: rel_change(
                    baseline.deformability_relative,
                    perturbed.deformability_relative,
                ),
            });
        }
    }
    rows
}

#[derive(Debug, Clone, Copy)]
struct Day42 {
    atp_mM: f64,
    deformability_relative: f64,
}

fn run_and_extract_day42(config: &StorageSimConfig) -> Day42 {
    let mut sim = StorageCurveSimulator::new(config.clone());
    sim.run();
    let s = sim.sample_at_day(42.0).expect("day-42 sample");
    Day42 {
        atp_mM: s.atp_mM,
        deformability_relative: s.deformability_relative,
    }
}

fn rel_change(baseline: f64, perturbed: f64) -> f64 {
    if baseline.abs() < 1e-12 {
        0.0
    } else {
        (perturbed - baseline) / baseline
    }
}

/// Write sensitivity rows to a CSV file with one header line and one
/// row per perturbation.
pub fn write_csv(rows: &[SensitivityRow], path: &Path) -> Result<()> {
    let mut wtr = csv::Writer::from_path(path)
        .with_context(|| format!("creating csv writer at {}", path.display()))?;
    wtr.write_record(&[
        "parameter",
        "baseline_value",
        "perturbed_value",
        "perturbation_pct",
        "day42_atp_baseline",
        "day42_atp_perturbed",
        "day42_atp_relative_change",
        "day42_deformability_baseline",
        "day42_deformability_perturbed",
        "day42_deformability_relative_change",
    ])?;
    for r in rows {
        wtr.write_record(&[
            r.parameter.name().to_string(),
            format!("{:.6}", r.baseline_value),
            format!("{:.6}", r.perturbed_value),
            format!("{:.2}", r.perturbation_pct),
            format!("{:.6}", r.day42_atp_baseline),
            format!("{:.6}", r.day42_atp_perturbed),
            format!("{:.6}", r.day42_atp_relative_change),
            format!("{:.6}", r.day42_deformability_baseline),
            format!("{:.6}", r.day42_deformability_perturbed),
            format!("{:.6}", r.day42_deformability_relative_change),
        ])?;
    }
    wtr.flush()?;
    Ok(())
}

/// All parameters in scope of Phase 14.D.3.
pub const ALL_PARAMETERS: &[SensitivityParameter] = &[
    SensitivityParameter::AtpHalfLifeDays,
    SensitivityParameter::PumpEfficiencyDecayPerDay,
    SensitivityParameter::LeakIncreasePerDay,
    SensitivityParameter::DpgLossRateMmPerDay,
    SensitivityParameter::OxidativeStressIncreasePerDay,
];

/// Find the top-N most-sensitive parameters by absolute fractional
/// change in day-42 deformability. Considers the larger of the +20% and
/// -20% perturbations per parameter.
pub fn top_sensitivities_by_deformability(
    rows: &[SensitivityRow],
    n: usize,
) -> Vec<SensitivityParameter> {
    use std::collections::HashMap;
    let mut max_per_param: HashMap<SensitivityParameter, f64> = HashMap::new();
    for r in rows {
        let abs_change = r.day42_deformability_relative_change.abs();
        let entry = max_per_param.entry(r.parameter).or_insert(0.0);
        if abs_change > *entry {
            *entry = abs_change;
        }
    }
    let mut entries: Vec<_> = max_per_param.into_iter().collect();
    entries.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
    entries.into_iter().take(n).map(|(p, _)| p).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn run_oat_sweep_returns_two_rows_per_parameter() {
        let cfg = StorageSimConfig {
            seconds_of_bio_per_step: 0.5,
            ..StorageSimConfig::default()
        };
        let rows = run_oat_sensitivity(&cfg, ALL_PARAMETERS, 0.20);
        assert_eq!(rows.len(), 2 * ALL_PARAMETERS.len());
        // Each parameter must have one +20% and one -20% perturbation.
        for &param in ALL_PARAMETERS {
            let pos = rows.iter().filter(|r| r.parameter == param && r.perturbation_pct > 0.0).count();
            let neg = rows.iter().filter(|r| r.parameter == param && r.perturbation_pct < 0.0).count();
            assert_eq!(pos, 1, "param {:?} missing +20% perturbation", param);
            assert_eq!(neg, 1, "param {:?} missing -20% perturbation", param);
        }
    }

    #[test]
    fn atp_half_life_dominates_atp_sensitivity() {
        let cfg = StorageSimConfig {
            seconds_of_bio_per_step: 0.5,
            ..StorageSimConfig::default()
        };
        let rows = run_oat_sensitivity(&cfg, ALL_PARAMETERS, 0.20);
        // The ATP half-life parameter should produce by far the largest
        // change in day-42 ATP — it's literally the parameter that
        // controls ATP decay.
        let atp_t_half_change = rows.iter()
            .filter(|r| r.parameter == SensitivityParameter::AtpHalfLifeDays)
            .map(|r| r.day42_atp_relative_change.abs())
            .fold(0.0_f64, f64::max);
        let other_changes_max = rows.iter()
            .filter(|r| r.parameter != SensitivityParameter::AtpHalfLifeDays)
            .map(|r| r.day42_atp_relative_change.abs())
            .fold(0.0_f64, f64::max);
        println!(
            "ATP sensitivity: t_half effect {:.3} vs others max {:.3}",
            atp_t_half_change, other_changes_max
        );
        assert!(
            atp_t_half_change > other_changes_max,
            "ATP half-life should dominate ATP sensitivity ({} vs {})",
            atp_t_half_change, other_changes_max
        );
    }

    #[test]
    fn write_csv_round_trips() {
        let cfg = StorageSimConfig {
            seconds_of_bio_per_step: 0.5,
            end_day: 42.0,
            ..StorageSimConfig::default()
        };
        let rows = run_oat_sensitivity(&cfg, &[SensitivityParameter::AtpHalfLifeDays], 0.20);
        let path = std::env::temp_dir().join("storage_sensitivity_test.csv");
        write_csv(&rows, &path).unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.starts_with("parameter,baseline_value"));
        // Header + 2 data rows.
        assert_eq!(content.lines().count(), 3);
        let _ = std::fs::remove_file(&path);
    }
}
