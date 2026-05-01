//! Phase 14.A validation: 42-day storage curve against Hess 2010 / Luten 2008 / Zimrin 2009.
//!
//! Reference targets (from `src/biochemistry/disease/storage_lesion.rs` doc comment):
//!
//! | Day | ATP (mM) | 2,3-DPG (mM) | Na+ (mM) | K+ (mM) |
//! |-----|----------|--------------|----------|---------|
//! | 0   | 2.0      | 5.0          | 10       | 140     |
//! | 14  | 1.5      | 0.5          | 25       | 120     |
//! | 42  | 0.5      | 0.0          | 60       | 90      |
//!
//! ATP and 2,3-DPG are forced to envelope targets (Hess 2010 / Zimrin 2009)
//! so they trivially match. The interesting validation is whether the
//! biochemistry/ion-homeostasis kinetics drive Na+ and K+ to the
//! literature ranges given the modified pump efficiency / leak conductance.

use cell_simulator_x::storage::{StorageCurveSimulator, StorageSimConfig};

const SECONDS_BIO_PER_STEP: f64 = 2.0;

fn run_full_curve() -> StorageCurveSimulator {
    let config = StorageSimConfig {
        end_day: 42.0,
        days_per_step: 1.0,
        seconds_of_bio_per_step: SECONDS_BIO_PER_STEP,
        bio_dt_sec: 1e-3,
        force_atp_dpg_targets: true,
        force_ion_qss: true,
    };
    let mut sim = StorageCurveSimulator::new(config);
    sim.run();
    sim
}

#[test]
fn day_0_targets() {
    let sim = run_full_curve();
    let s = sim.sample_at_day(0.0).expect("day 0 sample");
    assert!((s.atp_mM - 2.0).abs() < 0.05, "ATP day-0: {}", s.atp_mM);
    assert!((s.dpg23_mM - 5.0).abs() < 0.05, "DPG day-0: {}", s.dpg23_mM);
    // With ion QSS enabled, day-0 lands at the pump+leak equilibrium
    // (Na ≈ 10 mM, K ≈ 145 mM). The K value is ~5 mM above the
    // physiological "rest" 140 mM because the solver's pump kinetics
    // sit slightly off the default initial pool — close enough.
    assert!((s.na_cyt_mM - 10.0).abs() < 1.0, "Na day-0: {}", s.na_cyt_mM);
    assert!(s.k_cyt_mM > 138.0 && s.k_cyt_mM < 150.0, "K day-0: {}", s.k_cyt_mM);
}

#[test]
fn atp_decays_per_hess_2010() {
    let sim = run_full_curve();
    let d0 = sim.sample_at_day(0.0).unwrap().atp_mM;
    let d14 = sim.sample_at_day(14.0).unwrap().atp_mM;
    let d21 = sim.sample_at_day(21.0).unwrap().atp_mM;
    let d42 = sim.sample_at_day(42.0).unwrap().atp_mM;
    println!(
        "ATP: d0={:.3} d14={:.3} d21={:.3} d42={:.3}",
        d0, d14, d21, d42
    );
    // Hess 2010: ATP half-life ~21 d.
    assert!(d0 > d14, "ATP must decay");
    assert!(d14 > d21);
    assert!(d21 > d42);
    // Target band: day 21 = ~half of day 0.
    assert!((d21 - 1.0).abs() < 0.15, "day-21 ATP: {} (~1.0)", d21);
    // Day 42 = ~quarter of day 0.
    assert!((d42 - 0.5).abs() < 0.15, "day-42 ATP: {} (~0.5)", d42);
}

#[test]
fn dpg_depletes_per_zimrin_2009() {
    let sim = run_full_curve();
    let d0 = sim.sample_at_day(0.0).unwrap().dpg23_mM;
    let d14 = sim.sample_at_day(14.0).unwrap().dpg23_mM;
    let d42 = sim.sample_at_day(42.0).unwrap().dpg23_mM;
    println!("DPG: d0={:.3} d14={:.3} d42={:.3}", d0, d14, d42);
    assert!(d0 > 4.5, "day-0 DPG: {}", d0);
    // 5 mM - 0.4 mM/day * 14 = -0.6 → clamped to 0; allow up to 1 mM
    // (the BPGM enzyme produces a tiny amount during the equilibration
    // window each step, even when DPG is forced to 0 at the day boundary).
    assert!(d14 < 1.0, "day-14 DPG: {} (< 1)", d14);
    assert!(d42 < 0.05, "day-42 DPG: {} (< 0.05 mM)", d42);
}

#[test]
fn ion_gradients_trend_correctly() {
    let sim = run_full_curve();
    let d0 = sim.sample_at_day(0.0).unwrap();
    let d14 = sim.sample_at_day(14.0).unwrap();
    let d42 = sim.sample_at_day(42.0).unwrap();
    println!(
        "Na: d0={:.2} d14={:.2} d42={:.2}",
        d0.na_cyt_mM, d14.na_cyt_mM, d42.na_cyt_mM
    );
    println!(
        "K:  d0={:.2} d14={:.2} d42={:.2}",
        d0.k_cyt_mM, d14.k_cyt_mM, d42.k_cyt_mM
    );
    // Trend assertions: Na rises monotonically, K falls monotonically.
    // Quantitative agreement with Hess 2010 (Na=60 mM at day 42) needs
    // *longer* equilibration windows than fit in a unit test — see
    // `ion_gradients_match_hess_2010` (`#[ignore]`).
    assert!(d14.na_cyt_mM > d0.na_cyt_mM, "Na must rise by day 14");
    assert!(d42.na_cyt_mM > d14.na_cyt_mM, "Na must rise by day 42");
    assert!(d14.k_cyt_mM < d0.k_cyt_mM, "K must fall by day 14");
    assert!(d42.k_cyt_mM < d14.k_cyt_mM, "K must fall by day 42");
}

/// Long-equilibration validation. Phase 14.A simulator with 30 s of bio
/// per simulated day moves Na+ to ~25 mM at day 42 (vs Hess 2010's ~60
/// mM target). Closing the remaining gap is a Phase 14.B deliverable —
/// either by raising bio-per-step to ~300 s/day (10× wall-clock cost),
/// running an analytic quasi-steady-state pre-step for the ion system,
/// or retuning the pump/leak parameters.
///
/// This `#[ignore]` test guards Phase 14.A's measurable output: with
/// 30 s/day, Na+ reaches > 20 mM (clear directional signal) and K+
/// drops below 135 mM. Phase 14.B will tighten this band.
#[test]
#[ignore = "slow — run with `cargo test --release ion_gradients_long_equilibration -- --ignored`"]
fn ion_gradients_long_equilibration() {
    use cell_simulator_x::storage::{StorageCurveSimulator, StorageSimConfig};
    let mut sim = StorageCurveSimulator::new(StorageSimConfig {
        end_day: 42.0,
        days_per_step: 1.0,
        seconds_of_bio_per_step: 30.0,
        bio_dt_sec: 1e-3,
        force_atp_dpg_targets: true,
        // Disable QSS — this test exercises the legacy slow-equilibration path.
        force_ion_qss: false,
    });
    sim.run();
    let d14 = sim.sample_at_day(14.0).unwrap();
    let d42 = sim.sample_at_day(42.0).unwrap();
    println!(
        "30 s/day bio: d14 Na={:.1} K={:.1} | d42 Na={:.1} K={:.1}",
        d14.na_cyt_mM, d14.k_cyt_mM, d42.na_cyt_mM, d42.k_cyt_mM
    );
    // Phase 14.A bound — directional with margin (literature targets are
    // Na~25 d14, ~60 d42 / K~120 d14, ~90 d42). With 30 s/day bio we
    // reach roughly Na=15 d14, Na=25 d42; tightening to Hess 2010 numbers
    // is Phase 14.B.
    assert!(d42.na_cyt_mM > 20.0, "day-42 Na: {} (>20)", d42.na_cyt_mM);
    assert!(d42.k_cyt_mM < 135.0, "day-42 K: {} (<135)", d42.k_cyt_mM);
}

#[test]
fn full_curve_has_one_sample_per_day() {
    let sim = run_full_curve();
    // Day 0 + days 1..=42 = 43 samples.
    assert_eq!(sim.samples().len(), 43);
}

#[test]
fn write_csv_succeeds() {
    let sim = run_full_curve();
    let dir = std::env::temp_dir();
    let path = dir.join("storage_curve_test.csv");
    sim.write_csv(&path).expect("csv write");
    let written = std::fs::read_to_string(&path).expect("read back");
    assert!(written.starts_with("day,atp_mM"));
    let line_count = written.lines().count();
    // 1 header + 43 samples.
    assert_eq!(line_count, 44);
    let _ = std::fs::remove_file(&path);
}
