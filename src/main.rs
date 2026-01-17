//! Cell Simulator X - Entry point
//!
//! GPU-accelerated human red blood cell simulation engine.
//!
//! CLI Usage:
//!   cargo run                    # Run interactive simulation
//!   cargo run -- --diagnose      # Run physics diagnostics (no GUI)
//!   cargo run -- --diagnose -n 1000 -f 10.0  # Custom steps and force

use std::sync::Arc;
use std::time::Instant;

use anyhow::Result;
use cell_simulator_x::{
    config::Parameters,
    physics::{PhysicsConfig, PhysicsSolver},
    render::RenderState,
    state::CellState,
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

/// Parse CLI arguments
fn parse_args() -> (bool, usize, f32) {
    let args: Vec<String> = std::env::args().collect();
    let mut diagnose = false;
    let mut steps = 1000;
    let mut force = 5.0;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--diagnose" | "-d" => diagnose = true,
            "-n" | "--steps" => {
                i += 1;
                if i < args.len() {
                    steps = args[i].parse().unwrap_or(1000);
                }
            }
            "-f" | "--force" => {
                i += 1;
                if i < args.len() {
                    force = args[i].parse().unwrap_or(5.0);
                }
            }
            "--help" | "-h" => {
                println!("Cell Simulator X");
                println!();
                println!("Usage: cell-simulator-x [OPTIONS]");
                println!();
                println!("Options:");
                println!("  --diagnose, -d     Run physics diagnostics (no GUI)");
                println!("  -n, --steps N      Number of physics steps (default: 1000)");
                println!("  -f, --force F      Force magnitude in μN (default: 5.0)");
                println!("  --help, -h         Show this help");
                std::process::exit(0);
            }
            _ => {}
        }
        i += 1;
    }

    (diagnose, steps, force)
}

fn main() -> Result<()> {
    env_logger::init();

    // Parse CLI arguments
    let (diagnose, steps, force) = parse_args();

    if diagnose {
        return run_diagnostics(steps, force);
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
