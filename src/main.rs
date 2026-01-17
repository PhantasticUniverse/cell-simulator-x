// Allow non-snake-case for unit suffixes in field names (mM, mmHg, K, etc.)
#![allow(non_snake_case)]

//! Cell Simulator X - Entry point
//!
//! GPU-accelerated human red blood cell simulation engine.
//!
//! CLI Usage:
//!   cargo run                           # Run interactive simulation
//!   cargo run -- --diagnose             # Run physics diagnostics (no GUI)
//!   cargo run -- --diagnose -n 1000 -f 10.0  # Custom steps and force
//!   cargo run -- --diagnose-metabolism  # Run metabolism diagnostics (no GUI)
//!   cargo run -- --diagnose-metabolism -d 10.0  # 10 second simulation
//!   cargo run -- --diagnose-oxygen      # Run oxygen transport diagnostics (no GUI)
//!   cargo run -- --diagnose-oxygen --ph 7.0  # Test at different pH

use std::sync::Arc;
use std::time::Instant;

use anyhow::Result;
use cell_simulator_x::{
    config::Parameters,
    physics::{PhysicsConfig, PhysicsSolver},
    render::RenderState,
    state::CellState,
    MetabolismSolver, MetabolismConfig, MetabolitePool,
    HemoglobinSolver, HemoglobinState,
    STANDARD_PH, STANDARD_DPG_MM, STANDARD_PCO2_MMHG,
    run_integrated_diagnostics,
};
use glam::Vec3;
use winit::{
    event::{DeviceEvent, ElementState, Event, KeyEvent, MouseButton, WindowEvent},
    event_loop::{ControlFlow, EventLoop},
    keyboard::{KeyCode, PhysicalKey},
    window::WindowBuilder,
};

/// Run physics diagnostics without GUI
fn run_diagnostics(steps: usize, force_magnitude: f32) -> Result<()> {
    println!("=== Cell Simulator X - Physics Diagnostics ===\n");

    // Load parameters
    let params = Parameters::load_or_default();
    println!("Cell radius: {:.2} μm", params.geometry.cell_radius_um);

    // Create cell state
    let mut cell_state = CellState::new(&params);
    let n_vertices = cell_state.geometry.mesh.vertices.len();
    println!("Mesh vertices: {}", n_vertices);
    println!("Spectrin nodes: {}", cell_state.geometry.spectrin_network.nodes.len());

    // Store initial positions
    let initial_positions: Vec<Vec3> = cell_state
        .geometry
        .mesh
        .vertices
        .iter()
        .map(|v| v.position_vec3())
        .collect();

    // Find vertices on upper surface within a radius (simulating micropipette contact)
    let pipette_radius = 1.0; // 1 μm radius pipette
    let center = Vec3::new(0.0, 0.0, initial_positions.iter()
        .filter(|p| p.z > 0.0)
        .map(|p| p.z)
        .fold(0.0f32, f32::max)); // Top of cell

    let mut force_vertices: Vec<usize> = Vec::new();
    let mut center_idx = 0;
    let mut min_dist = f32::MAX;

    for (i, pos) in initial_positions.iter().enumerate() {
        if pos.z > 0.0 {
            let dist_xy = (pos.truncate() - center.truncate()).length();
            if dist_xy < pipette_radius {
                force_vertices.push(i);
            }
            if dist_xy < min_dist {
                min_dist = dist_xy;
                center_idx = i;
            }
        }
    }

    println!("Force target: {} vertices within {:.1} μm radius of center",
             force_vertices.len(), pipette_radius);
    println!("Center vertex {} at {:?}", center_idx, initial_positions[center_idx]);

    // Initialize physics solver
    let physics_config = PhysicsConfig {
        dt_sec: 1e-5,
        temperature_K: 310.0,
        enable_thermal_noise: false, // Disable noise for reproducible diagnostics
        membrane_damping: 5.0, // Damping for stable response
    };
    let mut physics_solver = PhysicsSolver::new(&cell_state.geometry.mesh, physics_config);

    // Apply force to single center vertex for clear behavior
    println!("Total force: {:.2} μN applied to center vertex", force_magnitude);
    println!("\n--- Running {} physics steps ---\n", steps);

    let force = Vec3::new(0.0, 0.0, -force_magnitude);  // Downward
    cell_state.physics.set_external_force(center_idx, force);

    // Run physics steps
    let start_time = Instant::now();
    for step in 0..steps {
        // Step physics
        physics_solver.step(
            &mut cell_state.geometry.mesh,
            &cell_state.geometry.spectrin_network,
            &mut cell_state.physics,
        );

        // Report progress every 10%
        if steps >= 10 && step % (steps / 10) == 0 {
            let progress = (step as f32 / steps as f32) * 100.0;
            let target_pos = cell_state.geometry.mesh.vertices[center_idx].position_vec3();
            let displacement = target_pos - initial_positions[center_idx];
            println!(
                "  {:3.0}%: step={}, target_z={:.4} μm, Δz={:.4} μm",
                progress, step, target_pos.z, displacement.z
            );
        }
    }
    let elapsed = start_time.elapsed();

    // Compute final statistics
    let final_positions: Vec<Vec3> = cell_state
        .geometry
        .mesh
        .vertices
        .iter()
        .map(|v| v.position_vec3())
        .collect();

    let displacements: Vec<f32> = initial_positions
        .iter()
        .zip(final_positions.iter())
        .map(|(i, f)| (*f - *i).length())
        .collect();

    let max_displacement = displacements.iter().cloned().fold(0.0f32, f32::max);
    let avg_displacement = displacements.iter().sum::<f32>() / n_vertices as f32;
    let target_displacement = final_positions[center_idx] - initial_positions[center_idx];

    let max_velocity = cell_state.physics.max_velocity_um_per_sec();
    let total_energy = cell_state.physics.total_energy_pJ();

    println!("\n=== Results ===");
    println!("Elapsed time: {:.2?}", elapsed);
    println!("Steps per second: {:.0}", steps as f32 / elapsed.as_secs_f32());
    println!("Simulation time: {:.4} ms", cell_state.physics.simulation_time_sec * 1000.0);
    println!();
    println!("Target vertex displacement: ({:.4}, {:.4}, {:.4}) μm",
        target_displacement.x, target_displacement.y, target_displacement.z);
    println!("Target vertex |displacement|: {:.4} μm", target_displacement.length());
    println!();
    println!("Max displacement (any vertex): {:.4} μm", max_displacement);
    println!("Avg displacement (all vertices): {:.6} μm", avg_displacement);
    println!("Max velocity: {:.2} μm/s", max_velocity);
    println!("Total energy: {:.4} pJ", total_energy);

    // Diagnostic checks
    println!("\n=== Diagnostic Checks ===");
    if max_displacement < 1e-6 {
        println!("⚠️  WARNING: No significant displacement detected!");
        println!("   Possible causes:");
        println!("   - Force too small (try -f 100)");
        println!("   - Damping too high");
        println!("   - Not enough steps (try -n 10000)");
    } else if max_displacement > 1.0 {
        println!("⚠️  WARNING: Very large displacement - simulation may be unstable");
    } else {
        println!("✓ Displacement looks reasonable");
    }

    if max_velocity < 1e-6 {
        println!("⚠️  WARNING: Velocities near zero - system may be over-damped");
    } else if max_velocity > 1e6 {
        println!("⚠️  WARNING: Very high velocities - simulation may be unstable");
    } else {
        println!("✓ Velocities look reasonable");
    }

    Ok(())
}

/// CLI arguments for the simulator
struct CliArgs {
    diagnose_physics: bool,
    diagnose_metabolism: bool,
    diagnose_oxygen: bool,
    diagnose_integrated: bool,
    steps: usize,
    force: f32,
    duration_sec: f64,
    glucose_step: Option<f64>,
    // Oxygen diagnostics parameters
    ph: f64,
    dpg_mM: f64,
    temperature_C: f64,
    pco2_mmHg: f64,
    // Integrated diagnostics parameters
    po2_mmHg: f64,
    atp_stress: f64,
}

impl Default for CliArgs {
    fn default() -> Self {
        Self {
            diagnose_physics: false,
            diagnose_metabolism: false,
            diagnose_oxygen: false,
            diagnose_integrated: false,
            steps: 1000,
            force: 5.0,
            duration_sec: 60.0,  // 60 seconds default for metabolism
            glucose_step: None,
            ph: STANDARD_PH,
            dpg_mM: STANDARD_DPG_MM,
            temperature_C: 37.0,
            pco2_mmHg: STANDARD_PCO2_MMHG,
            po2_mmHg: 100.0,     // Default arterial pO2
            atp_stress: 1.0,     // Default normal ATP consumption
        }
    }
}

/// Parse CLI arguments
fn parse_args() -> CliArgs {
    let args: Vec<String> = std::env::args().collect();
    let mut cli = CliArgs::default();

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--diagnose" => cli.diagnose_physics = true,
            "--diagnose-metabolism" => cli.diagnose_metabolism = true,
            "--diagnose-oxygen" => cli.diagnose_oxygen = true,
            "--diagnose-integrated" => cli.diagnose_integrated = true,
            "-n" | "--steps" => {
                i += 1;
                if i < args.len() {
                    cli.steps = args[i].parse().unwrap_or(1000);
                }
            }
            "-f" | "--force" => {
                i += 1;
                if i < args.len() {
                    cli.force = args[i].parse().unwrap_or(5.0);
                }
            }
            "-d" | "--duration" => {
                i += 1;
                if i < args.len() {
                    cli.duration_sec = args[i].parse().unwrap_or(60.0);
                }
            }
            "--glucose-step" => {
                i += 1;
                if i < args.len() {
                    cli.glucose_step = Some(args[i].parse().unwrap_or(10.0));
                }
            }
            "--ph" => {
                i += 1;
                if i < args.len() {
                    cli.ph = args[i].parse().unwrap_or(STANDARD_PH);
                }
            }
            "--dpg" => {
                i += 1;
                if i < args.len() {
                    cli.dpg_mM = args[i].parse().unwrap_or(STANDARD_DPG_MM);
                }
            }
            "--temp" => {
                i += 1;
                if i < args.len() {
                    cli.temperature_C = args[i].parse().unwrap_or(37.0);
                }
            }
            "--pco2" => {
                i += 1;
                if i < args.len() {
                    cli.pco2_mmHg = args[i].parse().unwrap_or(STANDARD_PCO2_MMHG);
                }
            }
            "--po2" => {
                i += 1;
                if i < args.len() {
                    cli.po2_mmHg = args[i].parse().unwrap_or(100.0);
                }
            }
            "--stress" => {
                i += 1;
                if i < args.len() {
                    cli.atp_stress = args[i].parse().unwrap_or(1.0);
                }
            }
            "--help" | "-h" => {
                println!("Cell Simulator X");
                println!();
                println!("Usage: cell-simulator-x [OPTIONS]");
                println!();
                println!("Options:");
                println!("  --diagnose             Run physics diagnostics (no GUI)");
                println!("  --diagnose-metabolism  Run metabolism diagnostics (no GUI)");
                println!("  --diagnose-oxygen      Run oxygen transport diagnostics (no GUI)");
                println!("  --diagnose-integrated  Run integrated metabolism-oxygen diagnostics (no GUI)");
                println!("  -n, --steps N          Number of physics steps (default: 1000)");
                println!("  -f, --force F          Force magnitude in μN (default: 5.0)");
                println!("  -d, --duration D       Duration in seconds for metabolism/integrated (default: 60.0)");
                println!("  --glucose-step G       Apply glucose step change at 50% duration");
                println!();
                println!("Oxygen diagnostics options:");
                println!("  --ph PH                Intracellular pH (default: 7.4)");
                println!("  --dpg DPG              2,3-DPG concentration in mM (default: 5.0)");
                println!("  --temp TEMP            Temperature in °C (default: 37)");
                println!("  --pco2 PCO2            CO2 partial pressure in mmHg (default: 40)");
                println!();
                println!("Integrated diagnostics options:");
                println!("  --po2 PO2              O2 partial pressure in mmHg (default: 100)");
                println!("  --stress S             ATP consumption multiplier (default: 1.0)");
                println!();
                println!("  --help, -h             Show this help");
                std::process::exit(0);
            }
            _ => {}
        }
        i += 1;
    }

    cli
}

/// Run metabolism diagnostics without GUI
fn run_metabolism_diagnostics(duration_sec: f64, glucose_step: Option<f64>) -> Result<()> {
    println!("=== Cell Simulator X - Metabolism Diagnostics ===\n");

    // Create metabolism solver
    let config = MetabolismConfig::default();
    let mut solver = MetabolismSolver::new(config);
    let mut metabolites = MetabolitePool::default_physiological();

    let indices = solver.indices;

    println!("Initial State:");
    println!("  ATP:      {:.3} mM", metabolites.get(indices.glycolysis.atp));
    println!("  ADP:      {:.3} mM", metabolites.get(indices.glycolysis.adp));
    println!("  2,3-DPG:  {:.3} mM", metabolites.get(indices.bisphosphoglycerate_2_3));
    println!("  Glucose:  {:.3} mM", metabolites.get(indices.glycolysis.glucose));
    println!("  Lactate:  {:.3} mM", metabolites.get(indices.glycolysis.lactate));
    println!();

    // ATP consumption rate (membrane pumps, etc.)
    // Reference: ~0.001 mM/s baseline consumption
    let atp_consumption = 0.001;

    // Time points for reporting
    let report_interval = duration_sec / 10.0;
    let mut next_report = report_interval;
    let step_time = if glucose_step.is_some() { duration_sec / 2.0 } else { f64::MAX };
    let mut step_applied = false;

    println!("--- Running {:.1} second simulation ---\n", duration_sec);
    println!("{:>8} {:>8} {:>8} {:>8} {:>8} {:>8}",
        "Time(s)", "ATP", "ADP", "2,3-DPG", "Glucose", "Lactate");
    println!("{:>8} {:>8} {:>8} {:>8} {:>8} {:>8}",
        "", "(mM)", "(mM)", "(mM)", "(mM)", "(mM)");
    println!("{}", "-".repeat(56));

    let start = Instant::now();

    // Integration loop
    let dt = solver.config.dt_sec;
    let n_steps = (duration_sec / dt).ceil() as usize;

    for step in 0..n_steps {
        let t = step as f64 * dt;

        // Apply glucose step if requested
        if let Some(glucose_change) = glucose_step {
            if t >= step_time && !step_applied {
                let current_glc = metabolites.get(indices.glycolysis.glucose);
                metabolites.set(indices.glycolysis.glucose, current_glc + glucose_change);
                println!("\n*** Glucose step: +{:.1} mM at t={:.1}s ***\n", glucose_change, t);
                step_applied = true;
            }
        }

        // Integrate
        solver.step(&mut metabolites, atp_consumption);

        // Report at intervals
        if t >= next_report || step == n_steps - 1 {
            println!("{:8.2} {:8.3} {:8.3} {:8.3} {:8.3} {:8.3}",
                solver.time_sec,
                metabolites.get(indices.glycolysis.atp),
                metabolites.get(indices.glycolysis.adp),
                metabolites.get(indices.bisphosphoglycerate_2_3),
                metabolites.get(indices.glycolysis.glucose),
                metabolites.get(indices.glycolysis.lactate));
            next_report += report_interval;
        }
    }

    let elapsed = start.elapsed();

    println!();
    println!("=== Final Diagnostics ===\n");

    let diag = solver.diagnostics(&metabolites);
    diag.print_summary();
    println!();
    diag.print_rates();

    // Validation checks
    println!("\n=== Validation Checks ===\n");
    let warnings = solver.validate_state(&metabolites);
    if warnings.is_empty() {
        println!("✓ All metabolite concentrations within physiological range");
    } else {
        for warning in &warnings {
            println!("⚠️  {}", warning);
        }
    }

    // Target validation
    let atp = metabolites.get(indices.glycolysis.atp);
    let dpg = metabolites.get(indices.bisphosphoglycerate_2_3);

    println!();
    if atp >= 1.5 && atp <= 2.5 {
        println!("✓ ATP concentration in target range (1.5-2.5 mM): {:.3} mM", atp);
    } else {
        println!("✗ ATP concentration outside target range (1.5-2.5 mM): {:.3} mM", atp);
    }

    if dpg >= 4.5 && dpg <= 5.5 {
        println!("✓ 2,3-DPG concentration in target range (4.5-5.5 mM): {:.3} mM", dpg);
    } else if dpg >= 3.0 && dpg <= 7.0 {
        println!("~ 2,3-DPG concentration acceptable (3.0-7.0 mM): {:.3} mM", dpg);
    } else {
        println!("✗ 2,3-DPG concentration outside acceptable range: {:.3} mM", dpg);
    }

    // Performance stats
    println!();
    println!("=== Performance ===");
    println!("Wall clock time: {:.2?}", elapsed);
    println!("Simulation time: {:.2} s", solver.time_sec);
    println!("Steps: {}", n_steps);
    println!("Steps/second: {:.0}", n_steps as f32 / elapsed.as_secs_f32());

    Ok(())
}

/// Run oxygen transport diagnostics without GUI
fn run_oxygen_diagnostics(ph: f64, dpg_mM: f64, temperature_C: f64, pco2_mmHg: f64) -> Result<()> {
    println!("=== Cell Simulator X - Oxygen Transport Diagnostics ===\n");

    let temperature_K = temperature_C + 273.15;
    let solver = HemoglobinSolver::default();

    // Print conditions
    println!("Conditions:");
    println!("  pH:          {:.2}", ph);
    println!("  2,3-DPG:     {:.1} mM", dpg_mM);
    println!("  Temperature: {:.1}°C ({:.1} K)", temperature_C, temperature_K);
    println!("  pCO2:        {:.1} mmHg", pco2_mmHg);
    println!();

    // Calculate key parameters
    let p50 = solver.calculate_p50(ph, dpg_mM, temperature_K, pco2_mmHg);
    let hill_n = solver.calculate_hill_coefficient(ph, dpg_mM, temperature_K, pco2_mmHg);
    let bohr_coeff = solver.measured_bohr_coefficient(dpg_mM, temperature_K, pco2_mmHg);

    println!("=== Calculated OEC Parameters ===\n");
    println!("P50:              {:.1} mmHg", p50);
    println!("Hill coefficient: {:.2}", hill_n);
    println!("Bohr coefficient: {:.2}", bohr_coeff);
    println!();

    // Validation against targets
    println!("=== Validation Against Imai 1982 ===\n");

    let p50_target = 26.8;
    let p50_tolerance = 1.0;
    let p50_pass = (p50 - p50_target).abs() <= p50_tolerance;
    println!("P50: {:.1} mmHg {} (target: {:.1} ± {:.1} mmHg)",
        p50,
        if p50_pass { "✓" } else { "✗" },
        p50_target, p50_tolerance);

    let hill_target = 2.7;
    let hill_tolerance = 0.1;
    let hill_pass = (hill_n - hill_target).abs() <= hill_tolerance;
    println!("Hill coefficient: {:.2} {} (target: {:.1} ± {:.1})",
        hill_n,
        if hill_pass { "✓" } else { "✗" },
        hill_target, hill_tolerance);

    let bohr_target = -0.48;
    let bohr_tolerance = 0.05;
    let bohr_pass = (bohr_coeff - bohr_target).abs() <= bohr_tolerance;
    println!("Bohr coefficient: {:.2} {} (target: {:.2} ± {:.2})",
        bohr_coeff,
        if bohr_pass { "✓" } else { "✗" },
        bohr_target, bohr_tolerance);

    // Generate and display OEC
    println!("\n=== Oxygen Equilibrium Curve ===\n");
    println!("{:>8} {:>12}", "pO2", "Saturation");
    println!("{:>8} {:>12}", "(mmHg)", "(%)");
    println!("{}", "-".repeat(22));

    let conditions = (ph, dpg_mM, temperature_K, pco2_mmHg);
    let oec_points = [0.0, 5.0, 10.0, 15.0, 20.0, 26.8, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0, 150.0];

    for &po2 in &oec_points {
        let sat = solver.calculate_saturation(po2, ph, dpg_mM, temperature_K, pco2_mmHg);
        let marker = if (po2 - p50).abs() < 1.0 { " ← P50" } else { "" };
        println!("{:8.1} {:11.1}%{}", po2, sat * 100.0, marker);
    }

    // Compare OEC at different pH values (Bohr effect demonstration)
    println!("\n=== Bohr Effect Comparison ===\n");
    println!("{:>8} {:>12} {:>12} {:>12}", "pO2", "pH 7.2", "pH 7.4", "pH 7.6");
    println!("{:>8} {:>12} {:>12} {:>12}", "(mmHg)", "Sat %", "Sat %", "Sat %");
    println!("{}", "-".repeat(50));

    let test_po2s = [10.0, 20.0, 26.8, 40.0, 60.0, 100.0];
    for &po2 in &test_po2s {
        let sat_72 = solver.calculate_saturation(po2, 7.2, dpg_mM, temperature_K, pco2_mmHg);
        let sat_74 = solver.calculate_saturation(po2, 7.4, dpg_mM, temperature_K, pco2_mmHg);
        let sat_76 = solver.calculate_saturation(po2, 7.6, dpg_mM, temperature_K, pco2_mmHg);
        println!("{:8.1} {:11.1}% {:11.1}% {:11.1}%", po2, sat_72 * 100.0, sat_74 * 100.0, sat_76 * 100.0);
    }

    // P50 comparison
    let p50_72 = solver.calculate_p50(7.2, dpg_mM, temperature_K, pco2_mmHg);
    let p50_74 = solver.calculate_p50(7.4, dpg_mM, temperature_K, pco2_mmHg);
    let p50_76 = solver.calculate_p50(7.6, dpg_mM, temperature_K, pco2_mmHg);
    println!();
    println!("P50 comparison:");
    println!("  pH 7.2: {:.1} mmHg", p50_72);
    println!("  pH 7.4: {:.1} mmHg", p50_74);
    println!("  pH 7.6: {:.1} mmHg", p50_76);

    // 2,3-DPG effect demonstration
    println!("\n=== 2,3-DPG Effect ===\n");
    let dpg_levels = [3.0, 5.0, 7.0];
    for &dpg in &dpg_levels {
        let p50_dpg = solver.calculate_p50(STANDARD_PH, dpg, temperature_K, pco2_mmHg);
        println!("  {:.1} mM 2,3-DPG: P50 = {:.1} mmHg", dpg, p50_dpg);
    }

    // Calculate 2,3-DPG sensitivity
    let p50_low_dpg = solver.calculate_p50(STANDARD_PH, 3.0, temperature_K, pco2_mmHg);
    let p50_high_dpg = solver.calculate_p50(STANDARD_PH, 7.0, temperature_K, pco2_mmHg);
    let dpg_sensitivity = (p50_high_dpg - p50_low_dpg) / (7.0 - 3.0);
    println!();
    println!("2,3-DPG sensitivity: {:.1} mmHg/mM (target: ~2-3 mmHg/mM)", dpg_sensitivity);

    // Dynamic binding test
    println!("\n=== Dynamic Binding Kinetics ===\n");
    let mut state = HemoglobinState::at_saturation(0.5, 5.0);
    println!("Starting saturation: {:.1}%", state.saturation * 100.0);
    println!("Applied pO2: 100 mmHg (high O2 environment)");
    println!();

    println!("{:>8} {:>12} {:>12}", "Time", "Saturation", "Bound O2");
    println!("{:>8} {:>12} {:>12}", "(ms)", "(%)", "(mM)");
    println!("{}", "-".repeat(36));

    let dt = 0.001;  // 1 ms timestep
    for i in 0..=10 {
        let time_ms = i as f64 * 10.0;
        println!("{:8.1} {:11.1}% {:11.2}", time_ms, state.saturation * 100.0, state.bound_o2_mM);
        // 10 steps per print
        for _ in 0..10 {
            solver.step(&mut state, 100.0, conditions, dt);
        }
    }

    println!();
    println!("Final saturation: {:.1}%", state.saturation * 100.0);

    // Overall validation summary
    println!("\n=== Validation Summary ===\n");
    let all_pass = p50_pass && hill_pass && bohr_pass;
    if all_pass {
        println!("✓ All validation targets met - Phase 4 complete!");
    } else {
        println!("Some validation targets not met:");
        if !p50_pass { println!("  - P50 outside tolerance"); }
        if !hill_pass { println!("  - Hill coefficient outside tolerance"); }
        if !bohr_pass { println!("  - Bohr coefficient outside tolerance"); }
    }

    Ok(())
}

fn main() -> Result<()> {
    env_logger::init();

    // Parse CLI arguments
    let cli = parse_args();

    if cli.diagnose_physics {
        return run_diagnostics(cli.steps, cli.force);
    }

    if cli.diagnose_metabolism {
        return run_metabolism_diagnostics(cli.duration_sec, cli.glucose_step);
    }

    if cli.diagnose_oxygen {
        return run_oxygen_diagnostics(cli.ph, cli.dpg_mM, cli.temperature_C, cli.pco2_mmHg);
    }

    if cli.diagnose_integrated {
        run_integrated_diagnostics(cli.duration_sec, cli.po2_mmHg, cli.atp_stress);
        return Ok(());
    }

    log::info!("Cell Simulator X starting...");

    // Load parameters
    let params = Parameters::load_or_default();
    log::info!("Parameters loaded: {:?}", params.geometry.cell_radius_um);

    // Create cell state with geometry
    let mut cell_state = CellState::new(&params);
    log::info!(
        "Cell state created: {} mesh vertices, {} spectrin nodes",
        cell_state.geometry.mesh.vertices.len(),
        cell_state.geometry.spectrin_network.nodes.len()
    );

    // Initialize physics solver
    let physics_config = PhysicsConfig {
        dt_sec: 1e-5,           // 10 microsecond timestep
        temperature_K: 310.0,    // 37°C body temperature
        enable_thermal_noise: true,
        membrane_damping: 5.0,   // Damping for stable visualization
    };
    let mut physics_solver = PhysicsSolver::new(&cell_state.geometry.mesh, physics_config);
    log::info!("Physics solver initialized");

    // Create window and event loop
    let event_loop = EventLoop::new()?;
    let window = Arc::new(
        WindowBuilder::new()
            .with_title("Cell Simulator X - Red Blood Cell")
            .with_inner_size(winit::dpi::LogicalSize::new(1280, 720))
            .build(&event_loop)?,
    );

    // Initialize render state
    let mut render_state = pollster::block_on(RenderState::new(window.clone(), &cell_state))?;

    // Input state
    let mut mouse_pressed = false;
    let mut show_spectrin = true;
    let mut physics_running = false;  // Start with physics paused
    let mut last_physics_time = Instant::now();
    let mut physics_substeps = 1;     // Number of physics substeps per frame (keep low for performance)

    // Force application state
    let mut apply_force = false;

    log::info!("Controls:");
    log::info!("  Mouse drag: Orbit camera");
    log::info!("  S: Toggle spectrin network");
    log::info!("  R: Reset camera");
    log::info!("  P: Toggle physics simulation");
    log::info!("  F: Apply force to center vertex");
    log::info!("  +/-: Adjust physics substeps");
    log::info!("  Escape: Exit");

    event_loop.run(move |event, elwt| {
        elwt.set_control_flow(ControlFlow::Poll);

        match event {
            Event::WindowEvent { event, .. } => match event {
                WindowEvent::CloseRequested => {
                    elwt.exit();
                }
                WindowEvent::KeyboardInput {
                    event:
                        KeyEvent {
                            physical_key: PhysicalKey::Code(key_code),
                            state: ElementState::Pressed,
                            ..
                        },
                    ..
                } => match key_code {
                    KeyCode::Escape => elwt.exit(),
                    KeyCode::KeyS => {
                        show_spectrin = !show_spectrin;
                        render_state.set_show_spectrin(show_spectrin);
                        log::info!("Spectrin overlay: {}", show_spectrin);
                    }
                    KeyCode::KeyR => {
                        render_state.camera.reset();
                        log::info!("Camera reset");
                    }
                    KeyCode::KeyP => {
                        physics_running = !physics_running;
                        log::info!(
                            "Physics simulation: {}",
                            if physics_running { "RUNNING" } else { "PAUSED" }
                        );
                        if physics_running {
                            last_physics_time = Instant::now();
                        }
                    }
                    KeyCode::KeyF => {
                        apply_force = !apply_force;
                        if apply_force {
                            // Find center vertex (closest to origin on upper surface)
                            let mut min_dist = f32::MAX;
                            let mut center_idx = 0;
                            for (i, v) in cell_state.geometry.mesh.vertices.iter().enumerate() {
                                let pos = v.position_vec3();
                                // Only consider upper surface (z > 0)
                                if pos.z > 0.0 {
                                    let dist = pos.truncate().length();
                                    if dist < min_dist {
                                        min_dist = dist;
                                        center_idx = i;
                                    }
                                }
                            }
                            // Apply persistent force (large enough for visible deformation)
                            let force = Vec3::new(0.0, 0.0, -100.0); // μN downward
                            physics_solver.apply_external_force(&mut cell_state.physics, center_idx, force);
                            log::info!("Applying {:.1} μN force to vertex {}", force.length(), center_idx);
                        } else {
                            physics_solver.clear_external_forces(&mut cell_state.physics);
                            log::info!("Force cleared");
                        }
                    }
                    KeyCode::Equal | KeyCode::NumpadAdd => {
                        physics_substeps = (physics_substeps + 5).min(100);
                        log::info!("Physics substeps: {}", physics_substeps);
                    }
                    KeyCode::Minus | KeyCode::NumpadSubtract => {
                        physics_substeps = (physics_substeps - 5).max(1);
                        log::info!("Physics substeps: {}", physics_substeps);
                    }
                    _ => {}
                },
                WindowEvent::MouseInput { state, button, .. } => {
                    if button == MouseButton::Left {
                        mouse_pressed = state == ElementState::Pressed;
                    }
                }
                WindowEvent::Resized(new_size) => {
                    render_state.resize(new_size);
                }
                WindowEvent::RedrawRequested => {
                    // Physics update
                    if physics_running {
                        let now = Instant::now();
                        let _frame_time = (now - last_physics_time).as_secs_f32();
                        last_physics_time = now;

                        // External forces are now persistent (set via F key)

                        // Run physics substeps
                        for _ in 0..physics_substeps {
                            physics_solver.step(
                                &mut cell_state.geometry.mesh,
                                &cell_state.geometry.spectrin_network,
                                &mut cell_state.physics,
                            );
                        }

                        // Update render buffers with new mesh
                        render_state.update_mesh(&cell_state.geometry.mesh);

                        // Log physics stats occasionally
                        if cell_state.physics.step_count % 1000 == 0 && cell_state.physics.step_count > 0 {
                            log::info!(
                                "Physics step {}: time={:.3}ms, energy={:.4}pJ, max_vel={:.2}μm/s",
                                cell_state.physics.step_count,
                                cell_state.physics.simulation_time_sec * 1000.0,
                                cell_state.physics.total_energy_pJ(),
                                cell_state.physics.max_velocity_um_per_sec()
                            );
                        }
                    }

                    // Render update
                    render_state.update();
                    match render_state.render() {
                        Ok(_) => {}
                        Err(wgpu::SurfaceError::Lost) => render_state.resize(render_state.size),
                        Err(wgpu::SurfaceError::OutOfMemory) => elwt.exit(),
                        Err(e) => log::error!("Render error: {:?}", e),
                    }
                }
                _ => {}
            },
            Event::DeviceEvent {
                event: DeviceEvent::MouseMotion { delta },
                ..
            } => {
                if mouse_pressed {
                    render_state
                        .camera
                        .orbit(delta.0 as f32 * 0.01, delta.1 as f32 * 0.01);
                }
            }
            Event::AboutToWait => {
                window.request_redraw();
            }
            _ => {}
        }
    })?;

    Ok(())
}
