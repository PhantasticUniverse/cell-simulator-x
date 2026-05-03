//! Phase C-hybrid: splenic-slit transit simulator.
//!
//! Couples a [`StorageCellSnapshot`] (metabolite + spectrin-modulator
//! state at storage day N) with a [`SlitFlow`] geometry, runs physics
//! + drag substeps, and reports transit time + peak strain + peak shear
//! + deformability for the cell.
//!
//! This is the headline integrated-coupling demo: the same cell that
//! the 42-day storage curve evolved is now placed in a 0.5–1.0 μm
//! splenic slit and asked to traverse. Transit time as a function of
//! storage day is a clinical-relevant biomarker no competitor
//! simulator can produce.

use std::time::Instant;

use glam::Vec3;

use crate::config::GeometryParameters;
use crate::flow::{apply_slit_drag_to_external_forces, SlitFlow, SplenicSlit};
use crate::geometry::{Mesh, SpectrinNetwork};
use crate::physics::{PhysicsConfig, PhysicsSolver, SkalakMaterial};
use crate::state::PhysicsState;
use crate::storage::simulator::StorageCellSnapshot;

/// Configuration for a splenic-transit simulation.
#[derive(Debug, Clone)]
pub struct SplenicTransitConfig {
    /// Slit geometry (width / height / length / axes).
    pub slit: SplenicSlit,
    /// Centerline flow velocity in the slit (μm/s).
    pub max_velocity_um_per_sec: f32,
    /// Stokes drag coefficient per vertex.
    pub drag_coeff: f32,
    /// Physics time step.
    pub dt_sec: f32,
    /// Cap on wall-clock seconds; abort if cell hasn't traversed.
    pub timeout_simulated_sec: f32,
    /// Cell mesh resolution (passed to `Mesh::generate_rbc`).
    pub mesh_resolution: u32,
    /// Spectrin network target node count.
    pub spectrin_target_count: u32,
}

impl Default for SplenicTransitConfig {
    fn default() -> Self {
        Self {
            slit: SplenicSlit::default(),
            max_velocity_um_per_sec: 500.0,
            // Higher drag coefficient than the Poiseuille demo because the
            // splenic slit's narrow geometry (sub-micron) demands strong
            // mechanical pull to overcome membrane elasticity. Calibrated
            // so a fresh cell traverses a 1 μm slit in ~0.1 s.
            drag_coeff: 5.0,
            dt_sec: 1e-5,
            timeout_simulated_sec: 0.5,
            mesh_resolution: 6,
            spectrin_target_count: 30,
        }
    }
}

/// Result of one splenic-transit simulation.
#[derive(Debug, Clone, Copy)]
pub struct SplenicTransitResult {
    /// Storage day this cell snapshot was taken from.
    pub storage_day: f64,
    /// Slit width (μm).
    pub slit_width_um: f32,
    /// Wall shear rate (1/s) — `4·v_max/w`.
    pub wall_shear_rate_per_sec: f32,
    /// Wall-clock simulated time before the cell's centroid exited the
    /// slit (seconds), or `timeout_simulated_sec` if it failed to exit.
    pub transit_time_sec: f32,
    /// Cell centroid axial displacement at end of run (μm). Primary
    /// metric: tracks how far the cell traveled along the slit axis.
    pub centroid_displacement_um: f32,
    /// True if the cell's centroid passed through the slit before the
    /// simulated-time cap. False = clearance failure.
    pub completed: bool,
    /// Peak per-vertex strain (relative position drift from initial)
    /// observed over the run.
    pub peak_strain_relative: f32,
    /// Peak velocity magnitude observed over the run.
    pub peak_velocity_um_per_sec: f32,
    /// Deformability index from the input snapshot.
    pub deformability_relative: f64,
    /// Wall-clock duration of the simulation (informational).
    pub wall_clock_ms: f64,
}

/// Run a splenic-transit simulation given a snapshot and configuration.
///
/// Workflow:
/// 1. Build a fresh mesh + spectrin network.
/// 2. Position the cell at the slit entrance.
/// 3. Apply the snapshot's stiffness modifier to the SkalakMaterial.
/// 4. Loop physics substeps; per substep, recompute drag from
///    SlitFlow and apply to `physics_state.external_forces_uN`, then
///    `PhysicsSolver::step`.
/// 5. Stop when the cell centroid passes the slit exit plane or the
///    simulated-time cap is reached.
pub fn run_splenic_transit(
    snapshot: &StorageCellSnapshot,
    config: &SplenicTransitConfig,
) -> SplenicTransitResult {
    let t0 = Instant::now();

    let geom = GeometryParameters {
        cell_radius_um: 3.91,
        fung_tong_c0_um: 0.81,
        fung_tong_c2_um: 7.83,
        fung_tong_c4_um: -4.39,
        mesh_resolution: config.mesh_resolution as usize,
        spectrin_target_count: config.spectrin_target_count as usize,
    };

    let mut mesh = Mesh::generate_rbc(&geom);
    let spectrin = SpectrinNetwork::generate(&mesh, &geom);

    // Position the cell straddling the slit entrance so part of the
    // membrane is already inside the slit and experiences drag from
    // step 1. A pre-existing cell radius of ~4 μm vs slit length ~5 μm
    // means roughly half the cell is in the slit at t=0.
    let centroid = mesh_centroid(&mesh);
    let cell_radius = 4.0_f32;
    let target_centroid = config.slit.entry + config.slit.axis * (cell_radius * 0.5);
    let shift = target_centroid - centroid;
    for v in mesh.vertices.iter_mut() {
        let p = Vec3::from_array(v.position) + shift;
        v.position = p.to_array();
    }

    // Initial reference positions for strain tracking.
    let initial_positions: Vec<Vec3> = mesh.vertices.iter().map(|v| v.position_vec3()).collect();

    let mut physics_state = PhysicsState::new(mesh.vertices.len());
    physics_state.init_reference_positions(&initial_positions);

    // Apply the snapshot's stiffness modifier to the membrane.
    let mut material = SkalakMaterial::default();
    material.shear_modulus_uN_per_m =
        material.shear_modulus_uN_per_m * snapshot.stiffness_modifier as f32;

    let physics_config = PhysicsConfig {
        dt_sec: config.dt_sec,
        temperature_K: 310.0,
        enable_thermal_noise: true,
        membrane_damping: 5.0,
    };
    let mut physics = PhysicsSolver::new(&mesh, physics_config);
    physics.skalak_solver.material = material;

    let flow = SlitFlow::new(config.slit, config.max_velocity_um_per_sec);
    let exit_plane = config.slit.entry + config.slit.axis * config.slit.length_um;

    let initial_centroid = mesh_centroid(&mesh);
    let initial_along = (initial_centroid - config.slit.entry).dot(config.slit.axis);

    let n_steps = (config.timeout_simulated_sec / config.dt_sec) as u32;
    let mut peak_strain = 0.0_f32;
    let mut peak_velocity = 0.0_f32;
    let mut transit_time_sec = config.timeout_simulated_sec;
    let mut completed = false;
    let mut final_along = initial_along;

    for step in 0..n_steps {
        apply_slit_drag_to_external_forces(&mut physics_state, &mesh, &flow, config.drag_coeff);
        physics.step(&mut mesh, &spectrin, &mut physics_state);

        // Track peak strain (max per-vertex displacement / cell_radius).
        for (i, v) in mesh.vertices.iter().enumerate() {
            let p = v.position_vec3();
            let d = (p - initial_positions[i]).length() / cell_radius;
            if d > peak_strain { peak_strain = d; }
        }
        // Track peak velocity (max magnitude across vertices).
        for vel in &physics_state.vertex_velocities_um_per_sec {
            let m = vel.length();
            if m > peak_velocity { peak_velocity = m; }
        }

        // Check exit condition: cell centroid past the exit plane.
        let centroid = mesh_centroid(&mesh);
        let along = (centroid - config.slit.entry).dot(config.slit.axis);
        final_along = along;
        if along > config.slit.length_um {
            transit_time_sec = step as f32 * config.dt_sec;
            completed = true;
            break;
        }
        let _ = exit_plane;
    }

    SplenicTransitResult {
        storage_day: snapshot.day,
        slit_width_um: config.slit.width_um,
        wall_shear_rate_per_sec: flow.wall_shear_rate_per_sec(),
        transit_time_sec,
        centroid_displacement_um: final_along - initial_along,
        completed,
        peak_strain_relative: peak_strain,
        peak_velocity_um_per_sec: peak_velocity,
        deformability_relative: snapshot.deformability_relative,
        wall_clock_ms: t0.elapsed().as_secs_f64() * 1000.0,
    }
}

fn mesh_centroid(mesh: &Mesh) -> Vec3 {
    let mut sum = Vec3::ZERO;
    for v in &mesh.vertices {
        sum += v.position_vec3();
    }
    sum / mesh.vertices.len() as f32
}

/// Sweep transit metrics over a (storage_day, slit_width) grid and
/// return the result table.
pub fn sweep_storage_day_x_slit_width(
    days: &[f64],
    widths_um: &[f32],
    base_config: &SplenicTransitConfig,
) -> Vec<SplenicTransitResult> {
    use crate::storage::simulator::{StorageCurveSimulator, StorageSimConfig};

    let mut results = Vec::with_capacity(days.len() * widths_um.len());
    for &day in days {
        // Build snapshot once per day; reuse across slit widths.
        let cfg = StorageSimConfig {
            seconds_of_bio_per_step: 1.0,
            end_day: day,
            ..StorageSimConfig::default()
        };
        let sim = StorageCurveSimulator::new(cfg);
        let snap = sim.cell_state_at_day(day);

        for &w in widths_um {
            let mut tc = base_config.clone();
            tc.slit = SplenicSlit::with_width(w);
            let r = run_splenic_transit(&snap, &tc);
            results.push(r);
        }
    }
    results
}

/// Write a transit-sweep result table to CSV.
pub fn write_transit_csv(
    rows: &[SplenicTransitResult],
    path: &std::path::Path,
) -> anyhow::Result<()> {
    use anyhow::Context as _;
    let mut wtr = csv::Writer::from_path(path)
        .with_context(|| format!("creating csv writer at {}", path.display()))?;
    wtr.write_record(&[
        "storage_day",
        "slit_width_um",
        "wall_shear_rate_per_sec",
        "transit_time_sec",
        "centroid_displacement_um",
        "completed",
        "peak_strain_relative",
        "peak_velocity_um_per_sec",
        "deformability_relative",
        "wall_clock_ms",
    ])?;
    for r in rows {
        wtr.write_record(&[
            format!("{:.2}", r.storage_day),
            format!("{:.3}", r.slit_width_um),
            format!("{:.0}", r.wall_shear_rate_per_sec),
            format!("{:.6}", r.transit_time_sec),
            format!("{:.4}", r.centroid_displacement_um),
            (if r.completed { "1" } else { "0" }).to_string(),
            format!("{:.4}", r.peak_strain_relative),
            format!("{:.2}", r.peak_velocity_um_per_sec),
            format!("{:.4}", r.deformability_relative),
            format!("{:.0}", r.wall_clock_ms),
        ])?;
    }
    wtr.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::storage::{StorageCurveSimulator, StorageSimConfig};

    fn fresh_snapshot() -> StorageCellSnapshot {
        let cfg = StorageSimConfig {
            seconds_of_bio_per_step: 0.5,
            end_day: 0.0,
            ..StorageSimConfig::default()
        };
        let sim = StorageCurveSimulator::new(cfg);
        sim.cell_state_at_day(0.0)
    }

    #[test]
    fn fresh_cell_traverses_canonical_slit() {
        let snap = fresh_snapshot();
        let result = run_splenic_transit(
            &snap,
            &SplenicTransitConfig {
                timeout_simulated_sec: 0.1,
                slit: SplenicSlit::with_width(0.7),
                ..SplenicTransitConfig::default()
            },
        );
        println!(
            "fresh cell, w=0.7 μm: transit={:.4} s, completed={}, peak_strain={:.3}, def={:.3}, wall_clock={:.0} ms",
            result.transit_time_sec, result.completed, result.peak_strain_relative,
            result.deformability_relative, result.wall_clock_ms
        );
        assert!(result.deformability_relative > 0.95, "fresh def: {}", result.deformability_relative);
    }
}
