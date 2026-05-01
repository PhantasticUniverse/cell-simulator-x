//! Phase 11.2.A/B parity test: GPU vs CPU glycolysis + 2,3-BPG shunt
//! after 1 s of simulated time.
//!
//! Builds an identical pool, runs the 18-species kernel on both backends
//! for 1000 RK4 ms-steps, and asserts per-species relative error < 1%
//! (or absolute < 1e-3 mM for sub-millimolar species). The CPU reference
//! mirrors the GPU `derivatives` function exactly so the comparison is
//! apples-to-apples (glycolysis + shunt + GLUT1 + MCT1 + ATP demand,
//! nothing more, nothing less).
//!
//! Skipped if no GPU adapter is available.

use cell_simulator_x::biochemistry::{
    GlycolysisSolver, MetaboliteIndices, MetabolitePool, RapoportLueberingSolver,
};
use cell_simulator_x::compute::{run_glycolysis_batch, ComputeContext, GlycolysisBatchConfig};
use cell_simulator_x::ExtendedMetaboliteIndices;

const DURATION_SEC: f64 = 1.0;
const DT_SEC: f64 = 0.001;

/// CPU reference: glycolysis + 2,3-BPG shunt RK4. Matches the GPU
/// kernel scope after Phase 11.2.B (shunt added).
#[allow(non_snake_case)]
fn run_cpu_glycolysis_only(pool: &mut MetabolitePool, duration_sec: f64) {
    let extended = ExtendedMetaboliteIndices::default();
    let glyco = GlycolysisSolver::new();
    let shunt = RapoportLueberingSolver::new(
        &extended.glycolysis,
        extended.bisphosphoglycerate_2_3,
    );
    // Match GPU defaults.
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
        // Closure that mirrors the GPU `derivatives` function.
        let derivatives = |state: &[f64], dydt: &mut [f64]| {
            for d in dydt.iter_mut() {
                *d = 0.0;
            }
            let temp = MetabolitePool {
                concentrations_mM: state.to_vec(),
            };
            glyco.compute_derivatives(&temp, dydt);
            shunt.compute_derivatives(&temp, dydt);
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

        // RK4 step (matches `RK4Integrator::step`).
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
    MetabolitePool::default_physiological()
}

fn extended() -> ExtendedMetaboliteIndices {
    ExtendedMetaboliteIndices::default()
}

fn species_label(idx: usize, indices: &MetaboliteIndices) -> &'static str {
    if idx == indices.glucose                  { "Glucose" }
    else if idx == indices.glucose_6_phosphate { "G6P" }
    else if idx == indices.fructose_6_phosphate { "F6P" }
    else if idx == indices.fructose_1_6_bisphosphate { "F1,6BP" }
    else if idx == indices.dihydroxyacetone_phosphate { "DHAP" }
    else if idx == indices.glyceraldehyde_3_phosphate { "GAP" }
    else if idx == indices.bisphosphoglycerate_1_3 { "1,3-BPG" }
    else if idx == indices.phosphoglycerate_3 { "3-PG" }
    else if idx == indices.phosphoglycerate_2 { "2-PG" }
    else if idx == indices.phosphoenolpyruvate { "PEP" }
    else if idx == indices.pyruvate { "Pyruvate" }
    else if idx == indices.lactate { "Lactate" }
    else if idx == indices.atp { "ATP" }
    else if idx == indices.adp { "ADP" }
    else if idx == indices.nad { "NAD+" }
    else if idx == indices.nadh { "NADH" }
    else if idx == indices.pi { "Pi" }
    else if idx == 17 { "2,3-BPG" }
    else { "(other)" }
}

#[test]
fn glycolysis_cpu_gpu_parity_one_second() {
    let ctx = match ComputeContext::new_headless_blocking() {
        Ok(c) => c,
        Err(e) => {
            eprintln!("skipping GPU parity test: no adapter: {e}");
            return;
        }
    };

    let indices = MetaboliteIndices::default();

    // === CPU run ===
    let mut cpu_pool = build_pool();
    run_cpu_glycolysis_only(&mut cpu_pool, DURATION_SEC);

    // === GPU run ===
    let mut gpu_pool = build_pool();
    let mut pools = vec![gpu_pool.clone()];
    let config = GlycolysisBatchConfig {
        dt_sec: DT_SEC,
        n_steps: (DURATION_SEC / DT_SEC).ceil() as u32,
        external_glucose_mM: 5.0,
        atp_consumption_mM_per_sec: 0.001,
        enable_glucose_transport: true,
        enable_lactate_export: true,
        max_change_mM: 0.1,
        min_concentration_mM: 1e-9,
    };
    run_glycolysis_batch(&ctx, &mut pools, &config).expect("GPU dispatch");
    gpu_pool = pools.into_iter().next().unwrap();

    // === Compare ===
    let ext = extended();
    let species: Vec<usize> = vec![
        indices.glucose,
        indices.glucose_6_phosphate,
        indices.fructose_6_phosphate,
        indices.fructose_1_6_bisphosphate,
        indices.dihydroxyacetone_phosphate,
        indices.glyceraldehyde_3_phosphate,
        indices.bisphosphoglycerate_1_3,
        indices.phosphoglycerate_3,
        indices.phosphoglycerate_2,
        indices.phosphoenolpyruvate,
        indices.pyruvate,
        indices.lactate,
        indices.atp,
        indices.adp,
        indices.nad,
        indices.nadh,
        indices.pi,
        ext.bisphosphoglycerate_2_3,
    ];

    let mut max_rel_err = 0.0f64;
    let mut max_rel_err_species = 0usize;
    let mut failures = Vec::new();
    println!("\n  {:<10} {:>14} {:>14} {:>10}", "species", "cpu (mM)", "gpu (mM)", "rel-err");
    for &idx in &species {
        let cpu_v = cpu_pool.get(idx);
        let gpu_v = gpu_pool.get(idx);
        let rel_err = if cpu_v.abs() > 1e-6 {
            ((cpu_v - gpu_v) / cpu_v).abs()
        } else {
            (cpu_v - gpu_v).abs()
        };
        let label = species_label(idx, &indices);
        println!("  {:<10} {:>14.6} {:>14.6} {:>10.4}", label, cpu_v, gpu_v, rel_err);
        if rel_err > max_rel_err {
            max_rel_err = rel_err;
            max_rel_err_species = idx;
        }
        // Tolerance: 1% relative for species above 1e-3 mM, 1e-3 absolute below.
        let pass = if cpu_v.abs() < 1e-3 {
            (cpu_v - gpu_v).abs() < 1e-3
        } else {
            rel_err < 0.01
        };
        if !pass {
            failures.push((label, cpu_v, gpu_v, rel_err));
        }
    }

    println!("\nWorst-fit: {} (rel-err = {:.4}%)",
        species_label(max_rel_err_species, &indices),
        max_rel_err * 100.0);

    assert!(
        failures.is_empty(),
        "GPU/CPU glycolysis parity failures (>{}% rel-err): {:?}",
        1.0,
        failures
    );
}
