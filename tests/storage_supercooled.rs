//! Phase 14.D.2 validation: supercooled (-4 °C) storage with Q10 = 2.5.
//!
//! Per Bruinsma 2019 / Berendsen 2014 supercooled storage protocols,
//! lowering RBC storage temperature from 4 °C to -4 °C extends viable
//! storage to ~100+ days. The simulator should reproduce this through
//! Q10 temperature scaling alone — no additional parameter retuning.
//!
//! Validation criterion (per Phase 14 plan):
//!   At -4 °C / day 100, deformability remains within 10% of the
//!   standard 4 °C / day 42 value.

use cell_simulator_x::storage::{
    AdditiveSolution, StorageCurveSimulator, StorageSimConfig,
};
use cell_simulator_x::storage::simulator::q10_rate_scale;

const SECONDS_BIO_PER_STEP: f64 = 1.0;

#[test]
fn q10_scale_matches_canonical_formula() {
    // Q10 = 2.5, T_ref = 4°C. At T = -4°C, scale = 2.5^(-0.8) ≈ 0.4774.
    let scale = q10_rate_scale(-4.0);
    assert!(
        (scale - 0.4774).abs() < 0.01,
        "Q10 scale at -4°C: {} (expected ≈ 0.4774)",
        scale
    );
    // Reference temperature: scale = 1.0.
    assert!((q10_rate_scale(4.0) - 1.0).abs() < 1e-9);
    // Higher T: scale > 1 (faster decay).
    assert!(q10_rate_scale(14.0) > 1.0);
}

#[test]
fn supercooled_extends_atp_lifetime() {
    let mut standard = StorageCurveSimulator::new(StorageSimConfig {
        seconds_of_bio_per_step: SECONDS_BIO_PER_STEP,
        ..StorageSimConfig::default()
    });
    standard.run();
    let mut sub = StorageCurveSimulator::new(StorageSimConfig {
        seconds_of_bio_per_step: SECONDS_BIO_PER_STEP,
        ..StorageSimConfig::supercooled(AdditiveSolution::Cpd)
    });
    sub.run();

    let standard_d42 = standard.sample_at_day(42.0).unwrap().atp_mM;
    let sub_d42 = sub.sample_at_day(42.0).unwrap().atp_mM;
    println!(
        "ATP @ d42: standard 4°C = {:.3} mM, supercooled -4°C = {:.3} mM",
        standard_d42, sub_d42
    );
    // Supercooled must preserve more ATP at the same calendar day.
    assert!(
        sub_d42 > standard_d42,
        "supercooled ATP at d42 ({:.3}) should exceed standard ({:.3})",
        sub_d42, standard_d42
    );
}

#[test]
fn supercooled_day100_deformability_matches_standard_day42() {
    // Phase 14.D.2 headline criterion: at -4°C / day 100, deformability
    // is within 10% of standard 4°C / day 42.
    let mut standard = StorageCurveSimulator::new(StorageSimConfig {
        seconds_of_bio_per_step: SECONDS_BIO_PER_STEP,
        ..StorageSimConfig::default()
    });
    standard.run();
    let standard_d42 = standard.sample_at_day(42.0).unwrap().deformability_relative;

    let mut sub = StorageCurveSimulator::new(StorageSimConfig {
        seconds_of_bio_per_step: SECONDS_BIO_PER_STEP,
        ..StorageSimConfig::supercooled(AdditiveSolution::Cpd)
    });
    sub.run();
    let sub_d100 = sub.sample_at_day(100.0).unwrap().deformability_relative;

    println!(
        "deformability: standard 4°C d42 = {:.3}, supercooled -4°C d100 = {:.3}",
        standard_d42, sub_d100
    );
    let rel_err = (sub_d100 - standard_d42).abs() / standard_d42;
    assert!(
        rel_err < 0.10,
        "supercooled d100 deformability {} differs from standard d42 {} by {:.1}% (>10%)",
        sub_d100, standard_d42, rel_err * 100.0
    );
}

#[test]
fn supercooled_pump_efficiency_higher_at_same_day() {
    // At any given storage day, supercooled storage should produce
    // higher pump efficiency than standard storage.
    let mut standard = StorageCurveSimulator::new(StorageSimConfig {
        seconds_of_bio_per_step: SECONDS_BIO_PER_STEP,
        ..StorageSimConfig::default()
    });
    standard.run();
    let std_d21 = standard.sample_at_day(21.0).unwrap().pump_efficiency;

    let mut sub = StorageCurveSimulator::new(StorageSimConfig {
        seconds_of_bio_per_step: SECONDS_BIO_PER_STEP,
        ..StorageSimConfig::supercooled(AdditiveSolution::Cpd)
    });
    sub.run();
    let sub_d21 = sub.sample_at_day(21.0).unwrap().pump_efficiency;

    println!(
        "pump efficiency d21: standard 4°C = {:.3}, supercooled -4°C = {:.3}",
        std_d21, sub_d21
    );
    assert!(
        sub_d21 > std_d21,
        "supercooled pump eff at d21 ({:.3}) should exceed standard ({:.3})",
        sub_d21, std_d21
    );
}

#[test]
fn supercooled_paggsm_outperforms_supercooled_cpd() {
    // Best-of-both-worlds: PAGGSM additive + supercooled storage should
    // preserve ATP better than CPD + supercooled.
    for additive in [AdditiveSolution::Cpd, AdditiveSolution::Paggsm] {
        let mut sim = StorageCurveSimulator::new(StorageSimConfig {
            seconds_of_bio_per_step: SECONDS_BIO_PER_STEP,
            ..StorageSimConfig::supercooled(additive)
        });
        sim.run();
        let d100 = sim.sample_at_day(100.0).unwrap();
        println!(
            "supercooled {} day-100: ATP {:.3} mM, def {:.3}",
            additive.name(), d100.atp_mM, d100.deformability_relative
        );
    }

    let cpd_atp = {
        let mut sim = StorageCurveSimulator::new(StorageSimConfig {
            seconds_of_bio_per_step: SECONDS_BIO_PER_STEP,
            ..StorageSimConfig::supercooled(AdditiveSolution::Cpd)
        });
        sim.run();
        sim.sample_at_day(100.0).unwrap().atp_mM
    };
    let paggsm_atp = {
        let mut sim = StorageCurveSimulator::new(StorageSimConfig {
            seconds_of_bio_per_step: SECONDS_BIO_PER_STEP,
            ..StorageSimConfig::supercooled(AdditiveSolution::Paggsm)
        });
        sim.run();
        sim.sample_at_day(100.0).unwrap().atp_mM
    };
    assert!(
        paggsm_atp > cpd_atp,
        "supercooled PAGGSM ATP at d100 ({:.3}) should exceed supercooled CPD ({:.3})",
        paggsm_atp, cpd_atp
    );
}

#[test]
fn standard_storage_unchanged_by_q10_at_4c() {
    // Sanity: at the reference temperature (4°C), Q10 scaling is 1.0
    // and the standard storage curve is identical to the Phase 14.B''
    // headline result.
    let mut sim = StorageCurveSimulator::new(StorageSimConfig {
        seconds_of_bio_per_step: SECONDS_BIO_PER_STEP,
        // Default temperature_celsius = 4.0
        ..StorageSimConfig::default()
    });
    sim.run();
    let d42 = sim.sample_at_day(42.0).unwrap();
    println!(
        "Standard 4°C day-42: ATP {:.3} mM, Na {:.1}, K {:.1}, def {:.3}",
        d42.atp_mM, d42.na_cyt_mM, d42.k_cyt_mM, d42.deformability_relative
    );
    // Hess 2010 anchors (must continue to match Phase 14.B'' headline).
    assert!(
        (d42.atp_mM - 0.5).abs() < 0.1,
        "day-42 ATP {} (Hess 2010 ≈ 0.5)",
        d42.atp_mM
    );
    assert!(
        (d42.na_cyt_mM - 60.0).abs() < 5.0,
        "day-42 Na {} (Hess 2010 ≈ 60)",
        d42.na_cyt_mM
    );
}
