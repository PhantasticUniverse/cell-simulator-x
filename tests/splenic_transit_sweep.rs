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

use cell_simulator_x::flow::SplenicSlit;
use cell_simulator_x::storage::{
    run_splenic_transit, sweep_storage_day_x_slit_width, write_transit_csv,
    SplenicTransitConfig, StorageCurveSimulator, StorageSimConfig,
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

/// Stream A drag-regime diagnostic: 3D sweep over (day, width, drag_coeff).
///
/// **Goal:** prove or refute the drag-saturation hypothesis. The
/// existing 2D sweep in `transit_sweep_emits_csv` shows only ~0.001 μm
/// difference between day-0 and day-42 at constant slit width. The
/// candidate explanations are:
///   (a) **drag-saturated** — the drag force at default coeff=5.0 dominates
///       membrane elasticity so much that the 37.5% stiffness delta
///       (modifier 1.0 → 1.375) between day-0 and day-42 is masked.
///   (b) **stiffness-saturated** — the stiffness change is genuinely too
///       small to matter regardless of drag.
///
/// Sweeping `drag_coeff ∈ {0.5, 1.0, 2.0, 5.0}` answers which: if lower
/// drag exposes a visible day-0 → day-42 displacement gap, hypothesis
/// (a) is confirmed and the 2D sweep was just operating in saturation.
///
/// Output: `target/splenic_transit_drag_sensitivity.csv` (36 rows = 3 days
/// × 3 widths × 4 drags).
#[test]
fn transit_drag_sensitivity_sweep_emits_csv() {
    let days = vec![0.0_f64, 21.0, 42.0];
    let widths = vec![0.5_f32, 0.7, 1.0];
    let drag_coeffs = vec![0.5_f32, 1.0, 2.0, 5.0];

    // Build the day-snapshots once, reuse across (width, drag) cells.
    // 3 × 3 × 4 = 36 simulations × ~150 ms each ≈ 5–6 s wall-clock.
    let mut rows: Vec<DragRow> = Vec::with_capacity(days.len() * widths.len() * drag_coeffs.len());

    for &day in &days {
        let cfg = StorageSimConfig {
            seconds_of_bio_per_step: 1.0,
            end_day: day,
            ..StorageSimConfig::default()
        };
        let sim = StorageCurveSimulator::new(cfg);
        let snap = sim.cell_state_at_day(day);

        for &w in &widths {
            for &drag in &drag_coeffs {
                let tc = SplenicTransitConfig {
                    timeout_simulated_sec: 0.1,
                    slit: SplenicSlit::with_width(w),
                    drag_coeff: drag,
                    ..SplenicTransitConfig::default()
                };
                let r = run_splenic_transit(&snap, &tc);
                rows.push(DragRow {
                    storage_day: r.storage_day,
                    slit_width_um: r.slit_width_um,
                    drag_coeff: drag,
                    wall_shear_rate_per_sec: r.wall_shear_rate_per_sec,
                    transit_time_sec: r.transit_time_sec,
                    centroid_displacement_um: r.centroid_displacement_um,
                    completed: r.completed,
                    peak_strain_relative: r.peak_strain_relative,
                    peak_velocity_um_per_sec: r.peak_velocity_um_per_sec,
                    deformability_relative: r.deformability_relative,
                });
            }
        }
    }

    assert_eq!(rows.len(), 36, "expected 36 sweep rows");

    let target_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("target");
    std::fs::create_dir_all(&target_dir).unwrap();
    let path = target_dir.join("splenic_transit_drag_sensitivity.csv");
    write_drag_csv(&rows, &path).expect("write drag-sensitivity csv");
    println!("Wrote drag-sensitivity sweep: {}", path.display());

    // Print human-readable summary table.
    println!(
        "{:>5} {:>5} {:>5} {:>8} {:>10} {:>8} {:>8} {:>5}",
        "day", "w_um", "drag", "shear", "displ_um", "strain", "v_peak", "ok"
    );
    for r in &rows {
        println!(
            "{:>5.0} {:>5.2} {:>5.2} {:>8.0} {:>10.4} {:>8.3} {:>8.1} {:>5}",
            r.storage_day,
            r.slit_width_um,
            r.drag_coeff,
            r.wall_shear_rate_per_sec,
            r.centroid_displacement_um,
            r.peak_strain_relative,
            r.peak_velocity_um_per_sec,
            if r.completed { "Y" } else { "N" }
        );
    }

    // Drag-vs-stiffness saturation diagnostic.
    //
    // For each (width, drag) pair compute the relative storage-day spread
    // (|Δdisplacement(d42 − d0)| / displacement(d0)). The hypothesis:
    //
    //   - If the system were drag-saturated, reducing drag_coeff from 5.0
    //     toward 0.5 would unmask the stiffness signal — we'd expect at
    //     least one cell with > 5% spread at the lowest drag.
    //   - If no (width, drag) cell exceeds 5%, the system is
    //     stiffness-saturated: the 37.5% stiffness modifier delta is
    //     intrinsically too small to translate into measurable centroid
    //     displacement, regardless of drag regime.
    println!("\n=== Drag-vs-stiffness saturation diagnostic ===");
    println!(
        "{:>5} {:>5} {:>10} {:>10} {:>12} {:>10}",
        "w_um", "drag", "day0_disp", "day42_disp", "rel_delta", "regime"
    );
    let mut max_rel_delta = 0.0_f32;
    let mut max_rel_delta_at: (f32, f32) = (0.0, 0.0);
    for &w in &widths {
        for &drag in &drag_coeffs {
            let day0 = rows
                .iter()
                .find(|r| {
                    r.storage_day == 0.0 && r.slit_width_um == w && r.drag_coeff == drag
                })
                .expect("day-0 row");
            let day42 = rows
                .iter()
                .find(|r| {
                    r.storage_day == 42.0 && r.slit_width_um == w && r.drag_coeff == drag
                })
                .expect("day-42 row");
            let d0 = day0.centroid_displacement_um;
            let d42 = day42.centroid_displacement_um;
            let rel_delta = if d0.abs() > 1e-9 { (d42 - d0) / d0 } else { 0.0 };
            // Per-cell label: > 5% spread = stiffness-sensitive at this
            // drag; otherwise the stiffness signal is below detection here.
            let regime = if rel_delta.abs() > 0.05 {
                "STIFF-SENS"
            } else {
                "stiff-sat"
            };
            println!(
                "{:>5.2} {:>5.2} {:>10.4} {:>10.4} {:>11.2}% {:>10}",
                w,
                drag,
                d0,
                d42,
                rel_delta * 100.0,
                regime
            );
            if rel_delta.abs() > max_rel_delta {
                max_rel_delta = rel_delta.abs();
                max_rel_delta_at = (w, drag);
            }
        }
    }
    println!(
        "\nMax |Δ(day42 − day0)/day0| = {:.2}% at (width={:.2} μm, drag={:.2})",
        max_rel_delta * 100.0,
        max_rel_delta_at.0,
        max_rel_delta_at.1
    );
    if max_rel_delta > 0.05 {
        println!(
            "Finding: DRAG-SATURATION REFUTED in favor of stiffness-sensitive regime. \
             The day-0/day-42 gap exceeds 5% at drag_coeff = {:.2} — lowering drag \
             unmasks the storage signal, confirming the original 7×3 sweep was \
             operating in drag saturation.",
            max_rel_delta_at.1
        );
    } else {
        println!(
            "Finding: STIFFNESS-SATURATED across the entire sweep. Reducing drag_coeff \
             10× (5.0 → 0.5) does NOT expose a > 5% storage-day spread anywhere in the \
             3×3×4 grid. The hypothesis that default drag=5.0 was masking stiffness \
             signal is REFUTED. The 37.5% spectrin-modifier delta (Phase 8 \
             max_stiffening_factor=0.5) is intrinsically below the centroid-displacement \
             detection floor at all tested drag regimes — amplifying the stiffness \
             coupling itself, not the drag, is the bottleneck."
        );
    }

    // Sanity: the CSV should contain 36 data rows + 1 header.
    let content = std::fs::read_to_string(&path).unwrap();
    assert!(content.starts_with("storage_day,slit_width_um,drag_coeff"));
    let line_count = content.lines().count();
    assert_eq!(line_count, 1 + 36);
}

/// Helper struct for the drag-sensitivity row (adds `drag_coeff` column).
struct DragRow {
    storage_day: f64,
    slit_width_um: f32,
    drag_coeff: f32,
    wall_shear_rate_per_sec: f32,
    transit_time_sec: f32,
    centroid_displacement_um: f32,
    completed: bool,
    peak_strain_relative: f32,
    peak_velocity_um_per_sec: f32,
    deformability_relative: f64,
}

fn write_drag_csv(rows: &[DragRow], path: &std::path::Path) -> std::io::Result<()> {
    use std::io::Write as _;
    let mut f = std::fs::File::create(path)?;
    writeln!(
        f,
        "storage_day,slit_width_um,drag_coeff,wall_shear_rate_per_sec,\
         transit_time_sec,centroid_displacement_um,completed,\
         peak_strain_relative,peak_velocity_um_per_sec,deformability_relative"
    )?;
    for r in rows {
        writeln!(
            f,
            "{:.2},{:.3},{:.3},{:.0},{:.6},{:.4},{},{:.4},{:.2},{:.4}",
            r.storage_day,
            r.slit_width_um,
            r.drag_coeff,
            r.wall_shear_rate_per_sec,
            r.transit_time_sec,
            r.centroid_displacement_um,
            if r.completed { 1 } else { 0 },
            r.peak_strain_relative,
            r.peak_velocity_um_per_sec,
            r.deformability_relative,
        )?;
    }
    Ok(())
}

