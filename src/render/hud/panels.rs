//! HUD panel definitions and rendering.

use egui::{Align2, Context, Window};

use super::state::HudState;
use super::theme::{HudColors, HudTypography};
use super::widgets::{key_value, section_header, DiseaseBadge, LabeledValue, O2Gauge, StatusBar};
use crate::state::{MetaboliteStatus, SimulationMetrics, SimulationMode};

/// Render all HUD panels
pub fn render_panels(ctx: &Context, state: &HudState, metrics: &SimulationMetrics) {
    if !state.hud_enabled {
        return;
    }

    if state.show_simulation_panel {
        render_simulation_panel(ctx, metrics);
    }

    if state.show_metabolites_panel {
        render_metabolites_panel(ctx, metrics);
    }

    if state.show_disease_panel {
        if let Some(ref disease) = metrics.disease {
            render_disease_panel(ctx, disease);
        }
    }

    if state.show_help {
        render_help_overlay(ctx);
    }

    if state.show_export_menu {
        render_export_menu(ctx);
    }
}

/// Simulation info panel (top-left)
fn render_simulation_panel(ctx: &Context, metrics: &SimulationMetrics) {
    Window::new("SIMULATION")
        .anchor(Align2::LEFT_TOP, [12.0, 12.0])
        .resizable(false)
        .collapsible(false)
        .title_bar(true)
        .show(ctx, |ui| {
            ui.set_min_width(140.0);

            // Time
            let time_str = if metrics.simulation_time_sec < 1.0 {
                format!("{:.2} ms", metrics.simulation_time_sec * 1000.0)
            } else if metrics.simulation_time_sec < 60.0 {
                format!("{:.2} s", metrics.simulation_time_sec)
            } else {
                format!("{:.1} min", metrics.simulation_time_sec / 60.0)
            };
            key_value(ui, "Time", &time_str, None);

            // FPS
            let fps_color = if metrics.fps >= 55.0 {
                HudColors::SUCCESS
            } else if metrics.fps >= 30.0 {
                HudColors::WARNING
            } else {
                HudColors::CRITICAL
            };
            key_value(ui, "FPS", &format!("{:.0}", metrics.fps), Some(fps_color));

            // Steps/second
            let steps_str = if metrics.steps_per_second >= 1000.0 {
                format!("{:.0}K/s", metrics.steps_per_second / 1000.0)
            } else {
                format!("{:.0}/s", metrics.steps_per_second)
            };
            key_value(ui, "Steps", &steps_str, None);

            ui.add_space(4.0);

            // Mode indicator
            let (mode_text, mode_color) = match metrics.mode {
                SimulationMode::Normal => ("Normal", HudColors::SUCCESS),
                SimulationMode::Stressed => ("Stressed", HudColors::WARNING),
                SimulationMode::Disease => ("Disease", HudColors::CRITICAL),
            };
            key_value(ui, "Mode", mode_text, Some(mode_color));

            // Physics status
            let physics_status = if metrics.physics_running { "ON" } else { "OFF" };
            let physics_color = if metrics.physics_running {
                HudColors::SUCCESS
            } else {
                HudColors::TEXT_SECONDARY
            };
            key_value(ui, "Physics", physics_status, Some(physics_color));

            // Substeps
            key_value(ui, "Substeps", &format!("{}", metrics.physics_substeps), None);
        });
}

/// Metabolites panel (top-right)
fn render_metabolites_panel(ctx: &Context, metrics: &SimulationMetrics) {
    Window::new("METABOLITES")
        .anchor(Align2::RIGHT_TOP, [-12.0, 12.0])
        .resizable(false)
        .collapsible(false)
        .title_bar(true)
        .show(ctx, |ui| {
            ui.set_min_width(180.0);

            // === Energy ===
            section_header(ui, "ENERGY");

            // ATP with status bar
            ui.add(
                LabeledValue::new("ATP", format!("{:.2}", metrics.atp_mM))
                    .unit("mM")
                    .status(metrics.atp_status),
            );
            ui.add(StatusBar::new(metrics.atp_mM, 0.5, 2.5, metrics.atp_status).target_range(1.5, 2.0));

            ui.add_space(2.0);

            // 2,3-DPG with status bar
            ui.add(
                LabeledValue::new("2,3-DPG", format!("{:.2}", metrics.dpg_2_3_mM))
                    .unit("mM")
                    .status(metrics.dpg_status),
            );
            ui.add(StatusBar::new(metrics.dpg_2_3_mM, 2.0, 8.0, metrics.dpg_status).target_range(4.5, 5.5));

            // === Oxygen ===
            section_header(ui, "OXYGEN");

            ui.horizontal(|ui| {
                // O2 saturation gauge
                ui.add(O2Gauge::new(metrics.o2_saturation).size(48.0));

                ui.vertical(|ui| {
                    // pO2
                    ui.add(LabeledValue::new("pO\u{2082}", format!("{:.0}", metrics.po2_mmHg)).unit("mmHg"));

                    // P50
                    ui.add(LabeledValue::new("P50", format!("{:.1}", metrics.p50_mmHg)).unit("mmHg"));
                });
            });

            // === pH ===
            section_header(ui, "ACID-BASE");

            ui.add(
                LabeledValue::new("pH", format!("{:.2}", metrics.ph)).status(metrics.ph_status),
            );

            // === Redox ===
            section_header(ui, "REDOX");

            // NADPH/NADP+
            ui.add(
                LabeledValue::new("NADPH/NADP\u{207A}", format!("{:.1}", metrics.nadph_nadp_ratio))
                    .status(metrics.redox_status),
            );

            // GSH/GSSG
            let gsh_display = if metrics.gsh_gssg_ratio > 999.0 {
                ">999".to_string()
            } else {
                format!("{:.0}", metrics.gsh_gssg_ratio)
            };
            ui.add(LabeledValue::new("GSH/GSSG", gsh_display).status(metrics.gsh_status));

            // H2O2
            ui.add(
                LabeledValue::new("H\u{2082}O\u{2082}", format!("{:.1}", metrics.h2o2_uM))
                    .unit("\u{03BC}M")
                    .status(metrics.h2o2_status),
            );

            // === Ions ===
            section_header(ui, "IONS");

            ui.add(
                LabeledValue::new("Na\u{207A}", format!("{:.1}", metrics.na_plus_mM))
                    .unit("mM")
                    .status(metrics.na_status),
            );

            ui.add(
                LabeledValue::new("K\u{207A}", format!("{:.0}", metrics.k_plus_mM))
                    .unit("mM")
                    .status(metrics.k_status),
            );

            ui.add(
                LabeledValue::new("Ca\u{00B2}\u{207A}", format!("{:.0}", metrics.ca_nM)).unit("nM"),
            );
        });
}

/// Disease panel (bottom-right)
fn render_disease_panel(ctx: &Context, disease: &crate::state::DiseaseIndicator) {
    Window::new("DISEASE")
        .anchor(Align2::RIGHT_BOTTOM, [-12.0, -12.0])
        .resizable(false)
        .collapsible(false)
        .title_bar(true)
        .show(ctx, |ui| {
            ui.set_min_width(200.0);

            // Disease badge
            ui.add(DiseaseBadge::new(&disease.name, disease.severity_percent));

            // Parameter
            ui.label(
                egui::RichText::new(&disease.parameter)
                    .size(HudTypography::LABEL_SIZE)
                    .color(HudColors::TEXT_SECONDARY),
            );

            ui.add_space(4.0);

            // Affected parameters
            section_header(ui, "AFFECTED");
            for (param, modifier) in &disease.affected_params {
                let modifier_str = format!("{:.2}x", modifier);
                let color = if *modifier < 0.8 {
                    HudColors::CRITICAL
                } else if *modifier < 1.0 {
                    HudColors::WARNING
                } else if *modifier > 1.2 {
                    HudColors::WARNING
                } else {
                    HudColors::TEXT_PRIMARY
                };
                key_value(ui, param, &modifier_str, Some(color));
            }
        });
}

/// Help overlay (center)
fn render_help_overlay(ctx: &Context) {
    Window::new("KEYBOARD SHORTCUTS")
        .anchor(Align2::CENTER_CENTER, [0.0, 0.0])
        .resizable(false)
        .collapsible(false)
        .title_bar(true)
        .show(ctx, |ui| {
            ui.set_min_width(280.0);

            section_header(ui, "CAMERA");
            key_value(ui, "Mouse Drag", "Orbit camera", None);
            key_value(ui, "R", "Reset camera", None);

            section_header(ui, "SIMULATION");
            key_value(ui, "P", "Toggle physics", None);
            key_value(ui, "F", "Apply/release force", None);
            key_value(ui, "+/-", "Adjust substeps", None);
            key_value(ui, "S", "Toggle spectrin", None);

            section_header(ui, "HUD");
            key_value(ui, "H", "Toggle help", None);
            key_value(ui, "E", "Export menu", None);
            key_value(ui, "M", "Toggle metabolites", None);
            key_value(ui, "D", "Toggle disease panel", None);
            key_value(ui, "Tab", "Toggle HUD", None);

            section_header(ui, "EXPORT");
            key_value(ui, "F12", "Screenshot", None);

            ui.add_space(8.0);
            ui.label(
                egui::RichText::new("Press H to close")
                    .size(HudTypography::SMALL_SIZE)
                    .color(HudColors::TEXT_SECONDARY),
            );
        });
}

/// Export menu (center)
fn render_export_menu(ctx: &Context) {
    Window::new("EXPORT")
        .anchor(Align2::CENTER_CENTER, [0.0, 0.0])
        .resizable(false)
        .collapsible(false)
        .title_bar(true)
        .show(ctx, |ui| {
            ui.set_min_width(200.0);

            ui.label(
                egui::RichText::new("Export Options")
                    .size(HudTypography::TITLE_SIZE)
                    .color(HudColors::TEXT_PRIMARY),
            );

            ui.add_space(8.0);

            if ui.button("Screenshot (PNG)").clicked() {
                // Will be handled by main
                log::info!("Screenshot requested");
            }

            if ui.button("State (JSON)").clicked() {
                log::info!("JSON export requested");
            }

            if ui.button("Time Series (CSV)").clicked() {
                log::info!("CSV export requested");
            }

            ui.add_space(8.0);
            ui.label(
                egui::RichText::new("Press E to close")
                    .size(HudTypography::SMALL_SIZE)
                    .color(HudColors::TEXT_SECONDARY),
            );
        });
}
