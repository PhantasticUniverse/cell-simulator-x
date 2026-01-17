//! Cell Simulator X - Entry point
//!
//! GPU-accelerated human red blood cell simulation engine.

use std::sync::Arc;

use anyhow::Result;
use cell_simulator_x::{config::Parameters, render::RenderState, state::CellState};
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
    let cell_state = CellState::new(&params);
    log::info!(
        "Cell state created: {} mesh vertices, {} spectrin nodes",
        cell_state.geometry.mesh.vertices.len(),
        cell_state.geometry.spectrin_network.nodes.len()
    );

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
