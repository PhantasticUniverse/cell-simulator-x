//! Storage additive solution presets (Phase 14.D.1).
//!
//! Each additive solution prolongs RBC viability through different
//! biochemical strategies (adenine for ATP precursor, citrate for pH
//! buffering, mannitol for membrane stabilization, guanosine for purine
//! salvage). Day-42 ATP retention is the headline literature anchor.
//!
//! | Additive | Day-42 ATP retention | T_half (days) | Reference |
//! |----------|----------------------|---------------|-----------|
//! | CPD baseline | ~25% | 21 | Hess 2010 (current calibration) |
//! | AS-3 (Nutricel) | ~70% | 81.6 | Hess 2010 Table 2 |
//! | SAGM | ~50% | 42 | Hess 2010 |
//! | PAGGSM | ~85% | 179 | Burger 2008, de Korte 2008 |
//!
//! Pump efficiency, leak conductance, and oxidative stress decay rates
//! are scaled inversely with retention as a first-pass approximation —
//! additives that preserve ATP also tend to preserve membrane integrity.
//! This is a simplifying assumption, refinable from per-additive
//! membrane-protein literature in a follow-on phase.

use crate::biochemistry::disease::StorageLesionConfig;

/// Storage additive solution preset.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AdditiveSolution {
    /// CPD baseline (no adenine). Day-42 ATP retention ~25%.
    /// Current default — matches the Phase 14.B'' headline calibration
    /// against Hess 2010.
    Cpd,
    /// AS-3 (Nutricel): adenine + citrate + dextrose + NaH₂PO₄.
    /// Day-42 ATP retention ~70% per Hess 2010 Table 2.
    As3,
    /// SAGM: saline + adenine + glucose + mannitol.
    /// Day-42 ATP retention ~50% per Hess 2010.
    Sagm,
    /// PAGGSM: phosphate + adenine + guanosine + glucose + saline +
    /// mannitol. Day-42 ATP retention ~85% per Burger 2008,
    /// de Korte 2008.
    Paggsm,
}

impl AdditiveSolution {
    /// Day-42 ATP retention fraction (0–1).
    pub fn day42_atp_retention(self) -> f64 {
        match self {
            Self::Cpd => 0.25,
            Self::As3 => 0.70,
            Self::Sagm => 0.50,
            Self::Paggsm => 0.85,
        }
    }

    /// Day-42 ATP target in mM (assuming ATP_0 = 2.0 mM).
    pub fn day42_atp_target_mM(self) -> f64 {
        2.0 * self.day42_atp_retention()
    }

    /// ATP exponential-decay half-life in days, derived from the day-42
    /// retention target via `T_half = -42·ln 2 / ln(retention)`.
    pub fn atp_half_life_days(self) -> f64 {
        let r = self.day42_atp_retention();
        -42.0 * std::f64::consts::LN_2 / r.ln()
    }

    /// Lesion envelope config calibrated for this additive.
    ///
    /// CPD returns the Phase 14.B'-calibrated baseline used throughout
    /// Phase 14. Other additives scale pump-decay and leak-increase rates
    /// inversely with retention — additives that preserve ATP also tend
    /// to preserve membrane-protein function.
    pub fn lesion_config(self) -> StorageLesionConfig {
        let t_half = self.atp_half_life_days();
        let retention = self.day42_atp_retention();
        // CPD baseline preservation factor = 1.0 (uses 0.0168/day pump
        // decay). Other additives reduce pump+leak decay by their
        // retention-ratio improvement over CPD.
        let preservation_scale = 0.25 / retention;

        StorageLesionConfig {
            atp_decay_half_life_days: t_half,
            // 2,3-DPG depletion is comparable across additives (Hess
            // 2010, Zimrin 2009 — all reach ~0 by day 14 regardless of
            // additive). Keep at baseline.
            dpg_loss_rate_mM_per_day: 0.4,
            pump_efficiency_decay_per_day: 0.0168 * preservation_scale,
            leak_increase_per_day: 0.015 * preservation_scale,
            // Oxidative stress mostly tracks plasma residuals + iron
            // chemistry; first-pass approximation keeps it at baseline.
            oxidative_stress_increase_per_day: 0.03,
            initial_atp_mM: 2.0,
            initial_dpg_mM: 5.0,
        }
    }

    /// Whether the Phase 14.B'' Hess-calibrated exponential pump envelope
    /// applies. Only the CPD baseline is calibrated; other additives use
    /// the linear envelope from `lesion_config`.
    ///
    /// Per-additive exponential calibrations are a clean follow-on once
    /// per-additive ion-gradient time-courses become available.
    pub fn uses_exponential_pump_envelope(self) -> bool {
        matches!(self, Self::Cpd)
    }

    /// Display name for diagnostics and CSV columns.
    pub fn name(self) -> &'static str {
        match self {
            Self::Cpd => "CPD",
            Self::As3 => "AS-3",
            Self::Sagm => "SAGM",
            Self::Paggsm => "PAGGSM",
        }
    }
}

impl Default for AdditiveSolution {
    fn default() -> Self {
        Self::Cpd
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cpd_matches_existing_calibration() {
        let cfg = AdditiveSolution::Cpd.lesion_config();
        assert!((cfg.atp_decay_half_life_days - 21.0).abs() < 0.1);
        assert!((cfg.pump_efficiency_decay_per_day - 0.0168).abs() < 1e-6);
        assert!((cfg.leak_increase_per_day - 0.015).abs() < 1e-6);
    }

    #[test]
    fn as3_t_half_matches_70pct_retention() {
        let cfg = AdditiveSolution::As3.lesion_config();
        // T_half = -42·ln2 / ln(0.70) ≈ 81.6 days.
        assert!(
            (cfg.atp_decay_half_life_days - 81.6).abs() < 1.0,
            "AS-3 T_half: {}",
            cfg.atp_decay_half_life_days
        );
        // Pump preservation: 0.0168 * (0.25/0.70) ≈ 0.006/day.
        assert!(
            (cfg.pump_efficiency_decay_per_day - 0.006).abs() < 1e-3,
            "AS-3 pump decay: {}",
            cfg.pump_efficiency_decay_per_day
        );
    }

    #[test]
    fn sagm_t_half_matches_50pct_retention() {
        let cfg = AdditiveSolution::Sagm.lesion_config();
        // T_half = 42 days exactly (50% retention at day 42).
        assert!(
            (cfg.atp_decay_half_life_days - 42.0).abs() < 0.1,
            "SAGM T_half: {}",
            cfg.atp_decay_half_life_days
        );
    }

    #[test]
    fn paggsm_t_half_matches_85pct_retention() {
        let cfg = AdditiveSolution::Paggsm.lesion_config();
        // T_half = -42·ln2 / ln(0.85) ≈ 179 days.
        assert!(
            (cfg.atp_decay_half_life_days - 179.0).abs() < 5.0,
            "PAGGSM T_half: {}",
            cfg.atp_decay_half_life_days
        );
    }

    #[test]
    fn day42_atp_targets_match_literature() {
        assert!((AdditiveSolution::Cpd.day42_atp_target_mM() - 0.5).abs() < 0.01);
        assert!((AdditiveSolution::As3.day42_atp_target_mM() - 1.4).abs() < 0.01);
        assert!((AdditiveSolution::Sagm.day42_atp_target_mM() - 1.0).abs() < 0.01);
        assert!((AdditiveSolution::Paggsm.day42_atp_target_mM() - 1.7).abs() < 0.01);
    }

    #[test]
    fn only_cpd_uses_exponential_envelope() {
        assert!(AdditiveSolution::Cpd.uses_exponential_pump_envelope());
        assert!(!AdditiveSolution::As3.uses_exponential_pump_envelope());
        assert!(!AdditiveSolution::Sagm.uses_exponential_pump_envelope());
        assert!(!AdditiveSolution::Paggsm.uses_exponential_pump_envelope());
    }
}
