//! Cell Simulator X - Entry point
//!
//! GPU-accelerated human red blood cell simulation engine.

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

fn main() -> Result<()> {
    env_logger::init();
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
        membrane_damping: 0.5,   // Increased damping for stability
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
    let mut force_vertex_idx: Option<usize> = None;

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
                            force_vertex_idx = Some(center_idx);
                            log::info!("Applying force to vertex {}", center_idx);
                        } else {
                            force_vertex_idx = None;
                            log::info!("Force application stopped");
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

                        // Apply external force if enabled
                        if let Some(idx) = force_vertex_idx {
                            // Apply downward force to simulate micropipette aspiration
                            let force = Vec3::new(0.0, 0.0, -0.001); // μN
                            physics_solver.apply_external_force(&mut cell_state.physics, idx, force);
                        }

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
