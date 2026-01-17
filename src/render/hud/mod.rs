//! HUD overlay system using egui.
//!
//! Provides a professional scientific instrument interface for displaying
//! simulation metrics and controls.

mod panels;
mod state;
mod theme;
mod widgets;

pub use state::HudState;
pub use theme::{HudColors, HudTheme};

use egui::Context;
use egui_wgpu::ScreenDescriptor;
use winit::event::WindowEvent;
use winit::window::Window;

use crate::state::SimulationMetrics;

/// Export action requested by the user through the HUD
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExportAction {
    Screenshot,
    JsonState,
    CsvTimeSeries,
}

/// HUD overlay manager integrating egui with wgpu
pub struct HudOverlay {
    /// Panel visibility state
    pub state: HudState,
    /// Theme configuration
    pub theme: HudTheme,
    /// egui context
    ctx: Context,
    /// egui-winit state
    egui_state: egui_winit::State,
    /// egui-wgpu renderer
    renderer: egui_wgpu::Renderer,
    /// Pending export action
    pending_export: Option<ExportAction>,
}

impl HudOverlay {
    /// Create a new HUD overlay
    pub fn new(
        window: &Window,
        device: &wgpu::Device,
        surface_format: wgpu::TextureFormat,
    ) -> Self {
        let ctx = Context::default();
        let theme = HudTheme::default();

        // Apply theme to context
        theme.apply(&ctx);

        // Create egui-winit state
        let viewport_id = ctx.viewport_id();
        let egui_state = egui_winit::State::new(
            ctx.clone(),
            viewport_id,
            window,
            Some(window.scale_factor() as f32),
            None,
        );

        // Create egui-wgpu renderer
        let renderer = egui_wgpu::Renderer::new(device, surface_format, None, 1);

        Self {
            state: HudState::default(),
            theme,
            ctx,
            egui_state,
            renderer,
            pending_export: None,
        }
    }

    /// Handle window events, returns true if egui consumed the event
    pub fn handle_event(&mut self, window: &Window, event: &WindowEvent) -> bool {
        let response = self.egui_state.on_window_event(window, event);
        response.consumed
    }

    /// Take any pending export action
    pub fn take_export_action(&mut self) -> Option<ExportAction> {
        self.pending_export.take()
    }

    /// Render the HUD overlay
    ///
    /// Returns paint jobs and texture delta for the renderer
    pub fn render(
        &mut self,
        window: &Window,
        metrics: &SimulationMetrics,
    ) -> (Vec<egui::ClippedPrimitive>, egui::TexturesDelta) {
        // Get input
        let raw_input = self.egui_state.take_egui_input(window);

        // Run the UI
        let output = self.ctx.run(raw_input, |ctx| {
            panels::render_panels(ctx, &self.state, metrics);
        });

        // Handle platform output (cursor changes, etc.)
        self.egui_state
            .handle_platform_output(window, output.platform_output);

        // Tessellate shapes
        let pixels_per_point = self.ctx.pixels_per_point();
        let primitives = self.ctx.tessellate(output.shapes, pixels_per_point);

        (primitives, output.textures_delta)
    }

    /// Paint the HUD to the screen
    pub fn paint(
        &mut self,
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        encoder: &mut wgpu::CommandEncoder,
        view: &wgpu::TextureView,
        screen_descriptor: ScreenDescriptor,
        paint_jobs: Vec<egui::ClippedPrimitive>,
        textures_delta: egui::TexturesDelta,
    ) {
        // Update textures
        for (id, image_delta) in &textures_delta.set {
            self.renderer
                .update_texture(device, queue, *id, image_delta);
        }

        // Update buffers
        self.renderer.update_buffers(
            device,
            queue,
            encoder,
            &paint_jobs,
            &screen_descriptor,
        );

        // Render
        {
            let mut render_pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("egui Render Pass"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Load, // Don't clear - render on top
                        store: wgpu::StoreOp::Store,
                    },
                })],
                depth_stencil_attachment: None,
                timestamp_writes: None,
                occlusion_query_set: None,
            });

            self.renderer.render(&mut render_pass, &paint_jobs, &screen_descriptor);
        }

        // Free textures
        for id in &textures_delta.free {
            self.renderer.free_texture(id);
        }
    }

    /// Get screen descriptor from window size
    pub fn screen_descriptor(&self, window: &Window) -> ScreenDescriptor {
        let size = window.inner_size();
        ScreenDescriptor {
            size_in_pixels: [size.width, size.height],
            pixels_per_point: window.scale_factor() as f32,
        }
    }

    // === State manipulation methods ===

    /// Toggle HUD visibility
    pub fn toggle_hud(&mut self) {
        self.state.toggle_hud();
    }

    /// Toggle help overlay
    pub fn toggle_help(&mut self) {
        self.state.toggle_help();
    }

    /// Toggle export menu
    pub fn toggle_export_menu(&mut self) {
        self.state.toggle_export_menu();
    }

    /// Toggle metabolites panel
    pub fn toggle_metabolites(&mut self) {
        self.state.toggle_metabolites();
    }

    /// Toggle disease panel
    pub fn toggle_disease(&mut self) {
        self.state.toggle_disease();
    }

    /// Request screenshot export
    pub fn request_screenshot(&mut self) {
        self.pending_export = Some(ExportAction::Screenshot);
    }

    /// Request JSON export
    pub fn request_json_export(&mut self) {
        self.pending_export = Some(ExportAction::JsonState);
    }

    /// Request CSV export
    pub fn request_csv_export(&mut self) {
        self.pending_export = Some(ExportAction::CsvTimeSeries);
    }

    /// Check if HUD wants to capture keyboard input
    pub fn wants_keyboard_input(&self) -> bool {
        self.ctx.wants_keyboard_input()
    }

    /// Check if HUD wants to capture mouse input
    pub fn wants_pointer_input(&self) -> bool {
        self.ctx.wants_pointer_input()
    }
}
