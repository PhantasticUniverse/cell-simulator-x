//! Phase 11.2.C parity test: GPU vs CPU full 38-species biochem
//! after 1 s of simulated time.
//!
//! Builds an identical pool, runs the full kernel on both backends for 1000
//! RK4 ms-steps, and asserts per-species relative error < 1% (or absolute
//! < 1e-3 mM for sub-millimolar species). The CPU reference mirrors the GPU
//! `derivatives` function exactly: glycolysis + 2,3-BPG shunt + PPP +
//! glutathione + Piezo1 + Na/K-ATPase + leaks + GLUT1 + MCT1 + external ATP
//! demand. NO inline homeostasis hacks (deferred to 11.2.E), NO Hb / pH
//! (deferred to 11.2.D).
//!
//! Skipped if no GPU adapter is available.

use cell_simulator_x::biochemistry::{
    GlutathioneCycle, GlycolysisSolver, IonHomeostasisSystem, MetaboliteIndices, MetabolitePool,
    PentosePhosphatePathway, Piezo1System, RapoportLueberingSolver,
};
use cell_simulator_x::compute::{run_full_biochem_batch, ComputeContext, FullBiochemBatchConfig};
use cell_simulator_x::FullyIntegratedIndices;

const DURATION_SEC: f64 = 1.0;
const DT_SEC: f64 = 0.001;
const MEMBRANE_TENSION: f64 = 0.5;

/// CPU reference: full biochem RK4 mirroring the GPU kernel scope exactly.
/// Glycolysis + 2,3-BPG shunt + PPP + glutathione + Piezo1 + ion homeostasis +
/// GLUT1 + MCT1 + external ATP demand. No Hb, no pH, no inline homeostasis hacks.
#[allow(non_snake_case)]
fn run_cpu_full_biochem(pool: &mut MetabolitePool, duration_sec: f64) {
    let full = FullyIntegratedIndices::new();
    let glyco = GlycolysisSolver::new();
    let shunt = RapoportLueberingSolver::new(&full.glycolysis, full.bisphosphoglycerate_2_3);
    let ppp = PentosePhosphatePathway::new(&full.redox);
    let mut glutathione = GlutathioneCycle::new(&full.redox);
    glutathione.h2o2_production_rate_mM_per_sec = 0.005;
    let piezo1 = Piezo1System::new(&full.redox);
    let ion_homeo = IonHomeostasisSystem::new(&full.glycolysis, &full.redox, 35);

    let external_glucose_mM = 5.0;
    let atp_consumption_mM_per_sec = 0.001;
    let max_change_mM = 0.1;
    let min_concentration_mM = 1e-9;
    let dt = DT_SEC;
    let n = pool.len();
    let indices = MetaboliteIndices::default();

    let mut k1 = vec![0f64; n];
    let mut k2 = vec![0f64; n];
    let mut k3 = vec![0f64; n];
    let mut k4 = vec![0f64; n];
    let mut y_temp = vec![0f64; n];

    let n_steps = (duration_sec / dt).ceil() as usize;

    for _ in 0..n_steps {
        let derivatives = |state: &[f64], dydt: &mut [f64]| {
            for d in dydt.iter_mut() {
                *d = 0.0;
            }
            let temp = MetabolitePool {
                concentrations_mM: state.to_vec(),
            };
            glyco.compute_derivatives(&temp, dydt);
            shunt.compute_derivatives(&temp, dydt);
            ppp.compute_derivatives(&temp, dydt);
            glutathione.compute_derivatives(&temp, dydt);
            piezo1.compute_derivatives(&temp, MEMBRANE_TENSION, dydt);
            ion_homeo.compute_derivatives(&temp, dydt);
            // External ATP consumption.
            dydt[indices.atp] -= atp_consumption_mM_per_sec;
            dydt[indices.adp] += atp_consumption_mM_per_sec;
            // GLUT1 transport.
            let glc = state[indices.glucose];
            dydt[indices.glucose] += 0.5 * (external_glucose_mM - glc);
            // MCT1 lactate export.
            let lac = state[indices.lactate];
            if lac > 1.0 {
                dydt[indices.lactate] -= 0.1 * (lac - 1.0);
            }
        };

        let y = pool.as_mut_slice();
        derivatives(y, &mut k1);
        for i in 0..n { y_temp[i] = y[i] + 0.5 * dt * k1[i]; }
        derivatives(&y_temp, &mut k2);
        for i in 0..n { y_temp[i] = y[i] + 0.5 * dt * k2[i]; }
        derivatives(&y_temp, &mut k3);
        for i in 0..n { y_temp[i] = y[i] + dt * k3[i]; }
        derivatives(&y_temp, &mut k4);
        let dt6 = dt / 6.0;
        for i in 0..n {
            let dy = dt6 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
            let dy_clamped = dy.clamp(-max_change_mM, max_change_mM);
            y[i] += dy_clamped;
            if y[i] < min_concentration_mM {
                y[i] = min_concentration_mM;
            }
        }
    }
}

fn build_pool() -> MetabolitePool {
    MetabolitePool::default_fully_integrated()
}

fn species_label(idx: usize) -> &'static str {
    match idx {
        0  => "Glucose",
        1  => "G6P",
        2  => "F6P",
        3  => "F1,6BP",
        4  => "DHAP",
        5  => "GAP",
        6  => "1,3-BPG",
        7  => "3-PG",
        8  => "2-PG",
        9  => "PEP",
        10 => "Pyruvate",
        11 => "Lactate",
        12 => "ATP",
        13 => "ADP",
        14 => "NAD+",
        15 => "NADH",
        16 => "Pi",
        17 => "2,3-BPG",
        18 => "6-PGL",
        19 => "6-PG",
        20 => "Ru5P",
        21 => "R5P",
        22 => "Xu5P",
        23 => "S7P",
        24 => "E4P",
        25 => "NADPH",
        26 => "NADP+",
        27 => "GSH",
        28 => "GSSG",
        29 => "H2O2",
        30 => "Ca2+ (uM)",
        31 => "Glu",
        32 => "Cys",
        33 => "Gly",
        34 => "γ-Glu-Cys",
        35 => "Na+",
        36 => "K+",
        37 => "Cl-",
        _  => "(other)",
    }
}

#[test]
fn full_biochem_cpu_gpu_parity_one_second() {
    let ctx = match ComputeContext::new_headless_blocking() {
        Ok(c) => c,
        Err(e) => {
            eprintln!("skipping GPU parity test: no adapter: {e}");
            return;
        }
    };

    // === CPU run ===
    let mut cpu_pool = build_pool();
    run_cpu_full_biochem(&mut cpu_pool, DURATION_SEC);

    // === GPU run ===
    let gpu_pool = build_pool();
    let mut pools = vec![gpu_pool.clone()];
    let config = FullBiochemBatchConfig {
        dt_sec: DT_SEC,
        n_steps: (DURATION_SEC / DT_SEC).ceil() as u32,
        membrane_tension_pN_per_nm: MEMBRANE_TENSION,
        enable_hb: false,
        ..FullBiochemBatchConfig::default()
    };
    run_full_biochem_batch(&ctx, &mut pools, &config).expect("GPU dispatch");
    let gpu_pool = pools.into_iter().next().unwrap();

    // === Compare ===
    let mut max_rel_err = 0.0f64;
    let mut max_rel_err_species = 0usize;
    let mut failures = Vec::new();
    println!(
        "\n  {:<11} {:>14} {:>14} {:>10}",
        "species", "cpu", "gpu", "rel-err"
    );
    for idx in 0..38 {
        let cpu_v = cpu_pool.get(idx);
        let gpu_v = gpu_pool.get(idx);
        let rel_err = if cpu_v.abs() > 1e-6 {
            ((cpu_v - gpu_v) / cpu_v).abs()
        } else {
            (cpu_v - gpu_v).abs()
        };
        let label = species_label(idx);
        println!(
            "  {:<11} {:>14.6} {:>14.6} {:>10.4}",
            label, cpu_v, gpu_v, rel_err
        );
        if rel_err > max_rel_err {
            max_rel_err = rel_err;
            max_rel_err_species = idx;
        }
        // Tolerance: 1% relative for species above 1e-3, 1e-3 absolute below.
        let pass = if cpu_v.abs() < 1e-3 {
            (cpu_v - gpu_v).abs() < 1e-3
        } else {
            rel_err < 0.01
        };
        if !pass {
            failures.push((label, cpu_v, gpu_v, rel_err));
        }
    }

    println!(
        "\nWorst-fit: {} (rel-err = {:.4}%)",
        species_label(max_rel_err_species),
        max_rel_err * 100.0
    );

    assert!(
        failures.is_empty(),
        "GPU/CPU full biochem parity failures (>{}% rel-err): {:?}",
        1.0,
        failures
    );
}
