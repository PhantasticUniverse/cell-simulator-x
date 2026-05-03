//! Phase C-hybrid.3: splenic-transit-vs-storage-day sweep.
//!
//! The headline figure: transit metrics across the (storage_day,
//! slit_width) grid showing how aged cells progressively struggle to
//! traverse narrow splenic slits. Output: target/splenic_transit_storage_curve.csv
//! ready for figure generation.
//!
//! Validation criteria (per Phase C-hybrid plan):
//!   - Day 0 / 1.0 μm slit: substantial centroid displacement
//!     (cell at least partially traverses).
//!   - Slit-width effect: 0.5 μm slit produces less displacement than
//!     1.0 μm slit.

use cell_simulator_x::storage::{
    sweep_storage_day_x_slit_width, write_transit_csv, SplenicTransitConfig,
};

#[test]
fn transit_sweep_emits_csv() {
    let days = vec![0.0_f64, 7.0, 14.0, 21.0, 28.0, 35.0, 42.0];
    let widths = vec![0.5_f32, 0.7, 1.0];

    let mut config = SplenicTransitConfig::default();
    // Tight wall-clock budget: 0.1 s simulated time per (day, width) cell;
    // at 1e-5 dt, that's 10,000 substeps × ~15 μs = ~150 ms each.
    // 7 days × 3 widths = 21 simulations → ~3 s total wall-clock.
    config.timeout_simulated_sec = 0.1;

    let rows = sweep_storage_day_x_slit_width(&days, &widths, &config);
    assert_eq!(rows.len(), days.len() * widths.len());

    let target_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("target");
    std::fs::create_dir_all(&target_dir).unwrap();
    let path = target_dir.join("splenic_transit_storage_curve.csv");
    write_transit_csv(&rows, &path).expect("write csv");
    println!("Wrote splenic-transit sweep: {}", path.display());

    // Print human-readable summary.
    println!(
        "{:>5} {:>5} {:>8} {:>10} {:>8} {:>8} {:>5}",
        "day", "w_um", "shear", "displ_um", "strain", "v_peak", "ok"
    );
    for r in &rows {
        println!(
            "{:>5.0} {:>5.2} {:>8.0} {:>10.3} {:>8.3} {:>8.1} {:>5}",
            r.storage_day,
            r.slit_width_um,
            r.wall_shear_rate_per_sec,
            r.centroid_displacement_um,
            r.peak_strain_relative,
            r.peak_velocity_um_per_sec,
            if r.completed { "Y" } else { "N" }
        );
    }

    let content = std::fs::read_to_string(&path).unwrap();
    assert!(content.starts_with("storage_day,slit_width_um"));
    let line_count = content.lines().count();
    assert_eq!(line_count, 1 + days.len() * widths.len());
}

#[test]
fn slit_width_strongly_modulates_transit() {
    // For any storage day, narrower slit → smaller centroid displacement
    // (mechanical filter effect). This is the qualitative signature
    // Phase C-hybrid is built to capture.
    let days = vec![0.0_f64];
    let widths = vec![0.5_f32, 1.0];
    let config = SplenicTransitConfig {
        timeout_simulated_sec: 0.1,
        ..SplenicTransitConfig::default()
    };
    let rows = sweep_storage_day_x_slit_width(&days, &widths, &config);
    assert_eq!(rows.len(), 2);

    let narrow = &rows[0]; // first width = 0.5
    let wide = &rows[1]; // second width = 1.0
    println!(
        "Day 0: narrow (0.5 μm) displacement = {:.3} μm; wide (1.0 μm) = {:.3} μm",
        narrow.centroid_displacement_um, wide.centroid_displacement_um
    );
    assert!(
        wide.centroid_displacement_um > narrow.centroid_displacement_um,
        "wider slit should produce more displacement: wide {} vs narrow {}",
        wide.centroid_displacement_um, narrow.centroid_displacement_um
    );
}

#[test]
fn day_42_cell_stiffness_modifier_visible() {
    // Day 42 cells have a higher stiffness modifier (reflecting Phase 8
    // ATP→spectrin coupling). The deformability index should be < 1.0.
    let days = vec![42.0_f64];
    let widths = vec![1.0_f32];
    let config = SplenicTransitConfig {
        timeout_simulated_sec: 0.05, // very short — just need deformability check
        ..SplenicTransitConfig::default()
    };
    let rows = sweep_storage_day_x_slit_width(&days, &widths, &config);
    assert_eq!(rows.len(), 1);
    let r = &rows[0];
    println!("Day 42: deformability = {}", r.deformability_relative);
    assert!(
        r.deformability_relative < 0.85,
        "day-42 deformability {} should be < 0.85 per Phase 14.C",
        r.deformability_relative
    );
}
