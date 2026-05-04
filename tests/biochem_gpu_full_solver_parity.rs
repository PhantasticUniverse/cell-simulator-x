//! Phase 11.2.E parity test: GPU vs the CPU `FullyIntegratedSolver` after
//! 1 s of simulated time. This is the strictest possible parity check —
//! the CPU side calls the same solver used in production
//! (`FullyIntegratedSolver::step`), and the GPU side enables every flag
//! (Piezo1 + ions + Hb + pH + the three inline homeostasis corrections).
//!
//! Asserts per-species relative error < 1% (or absolute < 1e-3 mM for
//! sub-millimolar species) AND hemoglobin saturation matches within 1%
//! absolute. Skipped if no GPU adapter is available.

use cell_simulator_x::biochemistry::{
    FullyIntegratedConfig, FullyIntegratedSolver, MetabolitePool,
};
use cell_simulator_x::compute::{
    run_full_biochem_batch_with_hb, ComputeContext, FullBiochemBatchConfig,
};

const DURATION_SEC: f64 = 1.0;
const DT_SEC: f64 = 0.001;
const MEMBRANE_TENSION: f64 = 0.5;
const PO2_MMHG: f64 = 100.0;
const INITIAL_HB_SAT: f64 = 0.75;

fn species_label(idx: usize) -> &'static str {
    match idx {
        0  => "Glucose",  1  => "G6P",      2  => "F6P",      3  => "F1,6BP",
        4  => "DHAP",     5  => "GAP",      6  => "1,3-BPG",  7  => "3-PG",
        8  => "2-PG",     9  => "PEP",      10 => "Pyruvate", 11 => "Lactate",
        12 => "ATP",      13 => "ADP",      14 => "NAD+",     15 => "NADH",
        16 => "Pi",       17 => "2,3-BPG",  18 => "6-PGL",    19 => "6-PG",
        20 => "Ru5P",     21 => "R5P",      22 => "Xu5P",     23 => "S7P",
        24 => "E4P",      25 => "NADPH",    26 => "NADP+",    27 => "GSH",
        28 => "GSSG",     29 => "H2O2",     30 => "Ca2+ uM",  31 => "Glu",
        32 => "Cys",      33 => "Gly",      34 => "γ-Glu-Cys",35 => "Na+",
        36 => "K+",       37 => "Cl-",      _  => "(other)",
    }
}

#[test]
fn gpu_matches_full_integrated_solver_one_second() {
    let ctx = match ComputeContext::new_headless_blocking() {
        Ok(c) => c,
        Err(e) => {
            eprintln!("skipping GPU parity test: no adapter: {e}");
            return;
        }
    };

    // === CPU run: the production FullyIntegratedSolver ===
    let mut cpu_pool = MetabolitePool::default_fully_integrated();
    let mut cpu_config = FullyIntegratedConfig::default();
    cpu_config.dt_sec = DT_SEC;
    cpu_config.po2_mmHg = PO2_MMHG;
    cpu_config.membrane_tension_pN_per_nm = MEMBRANE_TENSION;
    let mut cpu_solver = FullyIntegratedSolver::new(cpu_config);
    cpu_solver.hb_state =
        cell_simulator_x::biochemistry::HemoglobinState::at_saturation(INITIAL_HB_SAT, 5.0);
    let n_steps = (DURATION_SEC / DT_SEC).ceil() as usize;
    for _ in 0..n_steps {
        // basal_atp_consumption_mM_per_sec is built into the solver config;
        // pass it through `step` to match the production loop exactly.
        cpu_solver.step(&mut cpu_pool, cpu_solver.config.basal_atp_consumption_mM_per_sec);
    }
    let cpu_hb_sat = cpu_solver.hb_state.saturation;

    // === GPU run: the same parameter set, with all corrections enabled ===
    let mut pools = vec![MetabolitePool::default_fully_integrated()];
    let mut hb_sats = vec![INITIAL_HB_SAT];
    let config = FullBiochemBatchConfig {
        dt_sec: DT_SEC,
        n_steps: n_steps as u32,
        membrane_tension_pN_per_nm: MEMBRANE_TENSION,
        po2_mmHg: PO2_MMHG,
        // Match the FullyIntegratedConfig defaults verbatim.
        atp_consumption_mM_per_sec: cpu_solver.config.basal_atp_consumption_mM_per_sec,
        enable_homeostasis_corrections: true,
        ..FullBiochemBatchConfig::default()
    };
    run_full_biochem_batch_with_hb(&ctx, &mut pools, &mut hb_sats, &config)
        .expect("GPU dispatch");
    let gpu_pool = pools.into_iter().next().unwrap();
    let gpu_hb_sat = hb_sats[0];

    // === Compare ===
    let mut max_rel_err = 0.0f64;
    let mut max_rel_err_species: i32 = -1;
    let mut failures = Vec::new();
    println!("\n  {:<11} {:>14} {:>14} {:>10}", "species", "cpu", "gpu", "rel-err");
    for idx in 0..38 {
        let cpu_v = cpu_pool.get(idx);
        let gpu_v = gpu_pool.get(idx);
        let rel_err = if cpu_v.abs() > 1e-6 {
            ((cpu_v - gpu_v) / cpu_v).abs()
        } else {
            (cpu_v - gpu_v).abs()
        };
        let label = species_label(idx);
        println!("  {:<11} {:>14.6} {:>14.6} {:>10.4}", label, cpu_v, gpu_v, rel_err);
        if rel_err > max_rel_err {
            max_rel_err = rel_err;
            max_rel_err_species = idx as i32;
        }
        let pass = if cpu_v.abs() < 1e-3 {
            (cpu_v - gpu_v).abs() < 1e-3
        } else {
            rel_err < 0.01
        };
        if !pass {
            failures.push((label, cpu_v, gpu_v, rel_err));
        }
    }

    let hb_diff = (cpu_hb_sat - gpu_hb_sat).abs();
    println!(
        "\n  Hb sat      cpu = {:.6}  gpu = {:.6}  abs diff = {:.6}",
        cpu_hb_sat, gpu_hb_sat, hb_diff
    );

    let species_msg = if max_rel_err_species < 0 {
        "(none)".to_string()
    } else {
        format!(
            "{} (rel-err = {:.4}%)",
            species_label(max_rel_err_species as usize),
            max_rel_err * 100.0
        )
    };
    println!("\nWorst-fit metabolite: {}", species_msg);

    assert!(
        failures.is_empty(),
        "GPU vs FullyIntegratedSolver parity failures (>{}% rel-err): {:?}",
        1.0,
        failures
    );
    assert!(
        hb_diff < 0.01,
        "Hb saturation diverged: cpu = {:.6}, gpu = {:.6}, abs diff = {:.6}",
        cpu_hb_sat,
        gpu_hb_sat,
        hb_diff
    );
}
