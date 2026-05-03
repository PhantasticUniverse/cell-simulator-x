//! Phase 14.D.1 validation: additive solution comparator presets.
//!
//! Each additive solution (CPD baseline, AS-3, SAGM, PAGGSM) produces a
//! different day-42 ATP retention. The simulator's `AdditiveSolution`
//! enum drives `StorageLesionConfig` parameters (T_half, pump-decay,
//! leak-increase) per Hess 2010 / Burger 2008 / de Korte 2008 anchors.
//!
//! | Additive | Day-42 retention | Day-42 ATP target | Source |
//! |----------|------------------|-------------------|--------|
//! | CPD | 25% | 0.5 mM | Hess 2010 (current calibration) |
//! | AS-3 | 70% | 1.4 mM | Hess 2010 Table 2 |
//! | SAGM | 50% | 1.0 mM | Hess 2010 |
//! | PAGGSM | 85% | 1.7 mM | Burger 2008, de Korte 2008 |

use cell_simulator_x::storage::{
    AdditiveSolution, StorageCurveSimulator, StorageSimConfig,
};

const SECONDS_BIO_PER_STEP: f64 = 1.0;

fn run_additive(additive: AdditiveSolution) -> StorageCurveSimulator {
    let config = StorageSimConfig {
        seconds_of_bio_per_step: SECONDS_BIO_PER_STEP,
        ..StorageSimConfig::with_additive(additive)
    };
    let mut sim = StorageCurveSimulator::new(config);
    sim.run();
    sim
}

#[test]
fn cpd_baseline_day42_atp_matches_hess_2010() {
    // The CPD baseline must continue to match the headline result
    // (~0.5 mM at day 42 per Hess 2010 21-day half-life).
    let sim = run_additive(AdditiveSolution::Cpd);
    let s = sim.sample_at_day(42.0).unwrap();
    let target = AdditiveSolution::Cpd.day42_atp_target_mM();
    println!("CPD day-42 ATP: {:.3} (target {:.3})", s.atp_mM, target);
    let rel_err = (s.atp_mM - target).abs() / target;
    assert!(
        rel_err < 0.15,
        "CPD day-42 ATP {} mM vs target {} mM (rel err {:.1}%)",
        s.atp_mM, target, rel_err * 100.0
    );
}

#[test]
fn as3_day42_atp_matches_hess_2010_table_2() {
    let sim = run_additive(AdditiveSolution::As3);
    let s = sim.sample_at_day(42.0).unwrap();
    let target = AdditiveSolution::As3.day42_atp_target_mM(); // 1.4 mM
    println!("AS-3 day-42 ATP: {:.3} (target {:.3})", s.atp_mM, target);
    let rel_err = (s.atp_mM - target).abs() / target;
    assert!(
        rel_err < 0.15,
        "AS-3 day-42 ATP {} mM vs target {} mM (rel err {:.1}%)",
        s.atp_mM, target, rel_err * 100.0
    );
}

#[test]
fn sagm_day42_atp_matches_hess_2010() {
    let sim = run_additive(AdditiveSolution::Sagm);
    let s = sim.sample_at_day(42.0).unwrap();
    let target = AdditiveSolution::Sagm.day42_atp_target_mM(); // 1.0 mM
    println!("SAGM day-42 ATP: {:.3} (target {:.3})", s.atp_mM, target);
    let rel_err = (s.atp_mM - target).abs() / target;
    assert!(
        rel_err < 0.15,
        "SAGM day-42 ATP {} mM vs target {} mM (rel err {:.1}%)",
        s.atp_mM, target, rel_err * 100.0
    );
}

#[test]
fn paggsm_day42_atp_matches_burger_2008() {
    let sim = run_additive(AdditiveSolution::Paggsm);
    let s = sim.sample_at_day(42.0).unwrap();
    let target = AdditiveSolution::Paggsm.day42_atp_target_mM(); // 1.7 mM
    println!("PAGGSM day-42 ATP: {:.3} (target {:.3})", s.atp_mM, target);
    let rel_err = (s.atp_mM - target).abs() / target;
    assert!(
        rel_err < 0.15,
        "PAGGSM day-42 ATP {} mM vs target {} mM (rel err {:.1}%)",
        s.atp_mM, target, rel_err * 100.0
    );
}

#[test]
fn additive_atp_retention_is_monotonic() {
    // Better additives must produce monotonically higher day-42 ATP.
    // CPD < SAGM < AS-3 < PAGGSM per literature.
    let cpd = run_additive(AdditiveSolution::Cpd)
        .sample_at_day(42.0).unwrap().atp_mM;
    let sagm = run_additive(AdditiveSolution::Sagm)
        .sample_at_day(42.0).unwrap().atp_mM;
    let as3 = run_additive(AdditiveSolution::As3)
        .sample_at_day(42.0).unwrap().atp_mM;
    let paggsm = run_additive(AdditiveSolution::Paggsm)
        .sample_at_day(42.0).unwrap().atp_mM;
    println!(
        "Day-42 ATP: CPD={:.3} SAGM={:.3} AS-3={:.3} PAGGSM={:.3}",
        cpd, sagm, as3, paggsm
    );
    assert!(cpd < sagm, "CPD ATP {} should be < SAGM ATP {}", cpd, sagm);
    assert!(sagm < as3, "SAGM ATP {} should be < AS-3 ATP {}", sagm, as3);
    assert!(as3 < paggsm, "AS-3 ATP {} should be < PAGGSM ATP {}", as3, paggsm);
}

#[test]
fn additive_deformability_tracks_atp_retention() {
    // Higher day-42 ATP → lower spectrin stiffness → higher deformability.
    // This validates that the SpectrinModulator coupling reflects the
    // additive's metabolic preservation.
    let cpd = run_additive(AdditiveSolution::Cpd)
        .sample_at_day(42.0).unwrap().deformability_relative;
    let sagm = run_additive(AdditiveSolution::Sagm)
        .sample_at_day(42.0).unwrap().deformability_relative;
    let as3 = run_additive(AdditiveSolution::As3)
        .sample_at_day(42.0).unwrap().deformability_relative;
    let paggsm = run_additive(AdditiveSolution::Paggsm)
        .sample_at_day(42.0).unwrap().deformability_relative;
    println!(
        "Day-42 deformability: CPD={:.3} SAGM={:.3} AS-3={:.3} PAGGSM={:.3}",
        cpd, sagm, as3, paggsm
    );
    assert!(cpd < sagm, "CPD def {} should be < SAGM def {}", cpd, sagm);
    assert!(sagm < as3, "SAGM def {} should be < AS-3 def {}", sagm, as3);
    assert!(as3 < paggsm, "AS-3 def {} should be < PAGGSM def {}", as3, paggsm);
    // PAGGSM should preserve deformability close to fresh (>=0.9)
    // while CPD ends near 0.7 per the Phase 14.C calibration.
    assert!(
        paggsm >= 0.85,
        "PAGGSM day-42 deformability {} should be >= 0.85 (best preservation)",
        paggsm
    );
}

#[test]
fn additive_day0_states_match() {
    // All additives must agree at day 0 (no time evolution yet).
    let cpd = run_additive(AdditiveSolution::Cpd)
        .sample_at_day(0.0).unwrap().atp_mM;
    let as3 = run_additive(AdditiveSolution::As3)
        .sample_at_day(0.0).unwrap().atp_mM;
    let sagm = run_additive(AdditiveSolution::Sagm)
        .sample_at_day(0.0).unwrap().atp_mM;
    let paggsm = run_additive(AdditiveSolution::Paggsm)
        .sample_at_day(0.0).unwrap().atp_mM;
    println!(
        "Day-0 ATP: CPD={:.3} SAGM={:.3} AS-3={:.3} PAGGSM={:.3}",
        cpd, sagm, as3, paggsm
    );
    // All within 5% of 2.0 mM at day 0.
    for (name, atp) in [("CPD", cpd), ("AS-3", as3), ("SAGM", sagm), ("PAGGSM", paggsm)] {
        assert!(
            (atp - 2.0).abs() < 0.1,
            "{} day-0 ATP {} mM differs from baseline 2.0",
            name, atp
        );
    }
}

#[test]
fn additive_dpg_depletion_uniform_by_day_14() {
    // 2,3-DPG depletes by ~day 14 across all additives (Hess 2010,
    // Zimrin 2009 — additive solution does not strongly modulate
    // BPGM/BPGP kinetics).
    for additive in [
        AdditiveSolution::Cpd,
        AdditiveSolution::As3,
        AdditiveSolution::Sagm,
        AdditiveSolution::Paggsm,
    ] {
        let sim = run_additive(additive);
        let d14 = sim.sample_at_day(14.0).unwrap();
        assert!(
            d14.dpg23_mM < 1.0,
            "{} day-14 DPG: {} (should be < 1 mM)",
            additive.name(), d14.dpg23_mM
        );
    }
}
