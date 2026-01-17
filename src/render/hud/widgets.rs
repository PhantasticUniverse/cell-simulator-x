//! Custom HUD widgets for scientific data display.

use egui::{Color32, Pos2, Rect, Response, Sense, Stroke, Ui, Vec2, Widget};

use super::theme::{HudColors, HudTheme, HudTypography};
use crate::state::MetaboliteStatus;

/// A horizontal status bar showing a value within a range
pub struct StatusBar {
    /// Current value
    value: f64,
    /// Minimum value (left edge)
    min: f64,
    /// Maximum value (right edge)
    max: f64,
    /// Target range (shown as normal zone)
    target_min: f64,
    target_max: f64,
    /// Status color
    status: MetaboliteStatus,
    /// Width of the bar
    width: f32,
    /// Height of the bar
    height: f32,
}

impl StatusBar {
    /// Create a new status bar
    pub fn new(value: f64, min: f64, max: f64, status: MetaboliteStatus) -> Self {
        Self {
            value,
            min,
            max,
            target_min: min + (max - min) * 0.3,
            target_max: min + (max - min) * 0.7,
            status,
            width: 100.0,
            height: 8.0,
        }
    }

    /// Set target range (normal zone)
    pub fn target_range(mut self, target_min: f64, target_max: f64) -> Self {
        self.target_min = target_min;
        self.target_max = target_max;
        self
    }

    /// Set bar dimensions
    pub fn size(mut self, width: f32, height: f32) -> Self {
        self.width = width;
        self.height = height;
        self
    }
}

impl Widget for StatusBar {
    fn ui(self, ui: &mut Ui) -> Response {
        let (rect, response) = ui.allocate_exact_size(Vec2::new(self.width, self.height), Sense::hover());

        if ui.is_rect_visible(rect) {
            let painter = ui.painter();

            // Background
            painter.rect_filled(rect, 3.0, HudColors::BAR_BG);

            // Target zone indicator (subtle)
            let target_left = ((self.target_min - self.min) / (self.max - self.min)) as f32;
            let target_right = ((self.target_max - self.min) / (self.max - self.min)) as f32;
            let target_rect = Rect::from_min_max(
                Pos2::new(rect.left() + rect.width() * target_left, rect.top()),
                Pos2::new(rect.left() + rect.width() * target_right, rect.bottom()),
            );
            painter.rect_filled(target_rect, 3.0, HudColors::SUCCESS.gamma_multiply(0.15));

            // Fill based on value
            let fill_pct = ((self.value - self.min) / (self.max - self.min)).clamp(0.0, 1.0) as f32;
            let fill_rect = Rect::from_min_max(
                rect.left_top(),
                Pos2::new(rect.left() + rect.width() * fill_pct, rect.bottom()),
            );

            let fill_color = HudTheme::status_color(self.status);
            painter.rect_filled(fill_rect, 3.0, fill_color.gamma_multiply(0.8));

            // Border
            painter.rect_stroke(rect, 3.0, Stroke::new(1.0, HudColors::BORDER));
        }

        response
    }
}

/// Labeled value display with status indicator
pub struct LabeledValue {
    label: String,
    value: String,
    unit: Option<String>,
    status: Option<MetaboliteStatus>,
}

impl LabeledValue {
    /// Create a new labeled value
    pub fn new(label: impl Into<String>, value: impl Into<String>) -> Self {
        Self {
            label: label.into(),
            value: value.into(),
            unit: None,
            status: None,
        }
    }

    /// Add unit suffix
    pub fn unit(mut self, unit: impl Into<String>) -> Self {
        self.unit = Some(unit.into());
        self
    }

    /// Add status color
    pub fn status(mut self, status: MetaboliteStatus) -> Self {
        self.status = Some(status);
        self
    }
}

impl Widget for LabeledValue {
    fn ui(self, ui: &mut Ui) -> Response {
        let response = ui.horizontal(|ui| {
            // Label (dim)
            ui.label(
                egui::RichText::new(&self.label)
                    .size(HudTypography::LABEL_SIZE)
                    .color(HudColors::TEXT_LABEL),
            );

            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                // Status indicator (small colored circle)
                if let Some(status) = self.status {
                    let color = HudTheme::status_color(status);
                    let (indicator_rect, _) = ui.allocate_exact_size(egui::vec2(8.0, 8.0), egui::Sense::hover());
                    ui.painter().circle_filled(indicator_rect.center(), 3.0, color);
                }

                // Unit (small, muted)
                if let Some(ref unit) = self.unit {
                    ui.label(
                        egui::RichText::new(unit)
                            .size(HudTypography::UNIT_SIZE)
                            .color(HudColors::TEXT_SECONDARY),
                    );
                }

                // Value (monospace, bright)
                let value_color = self.status.map(HudTheme::status_color).unwrap_or(HudColors::TEXT_PRIMARY);
                ui.label(
                    egui::RichText::new(&self.value)
                        .size(HudTypography::VALUE_SIZE)
                        .family(egui::FontFamily::Monospace)
                        .color(value_color),
                );
            });
        });

        response.response
    }
}

/// Circular oxygen saturation gauge
pub struct O2Gauge {
    saturation: f64,
    size: f32,
}

impl O2Gauge {
    pub fn new(saturation: f64) -> Self {
        Self {
            saturation,
            size: 48.0,
        }
    }

    pub fn size(mut self, size: f32) -> Self {
        self.size = size;
        self
    }
}

impl Widget for O2Gauge {
    fn ui(self, ui: &mut Ui) -> Response {
        let (rect, response) = ui.allocate_exact_size(Vec2::splat(self.size), Sense::hover());

        if ui.is_rect_visible(rect) {
            let painter = ui.painter();
            let center = rect.center();
            let radius = self.size / 2.0 - 2.0;

            // Background arc
            painter.circle_stroke(center, radius, Stroke::new(4.0, HudColors::BAR_BG));

            // Saturation arc
            let sat_color = HudTheme::o2_saturation_color(self.saturation);
            let sweep = (self.saturation * 360.0) as f32;

            // Draw arc using line segments
            let start_angle = -90.0_f32.to_radians();
            let segments = 32;
            for i in 0..segments {
                let t = i as f32 / segments as f32;
                let angle_deg = t * sweep;
                if angle_deg > sweep {
                    break;
                }
                let angle1 = start_angle + (t * sweep).to_radians();
                let angle2 = start_angle + ((t + 1.0 / segments as f32) * sweep).min(sweep).to_radians();

                let p1 = center + Vec2::new(angle1.cos(), angle1.sin()) * radius;
                let p2 = center + Vec2::new(angle2.cos(), angle2.sin()) * radius;

                painter.line_segment([p1, p2], Stroke::new(4.0, sat_color));
            }

            // Center text
            let text = format!("{:.0}%", self.saturation * 100.0);
            painter.text(
                center,
                egui::Align2::CENTER_CENTER,
                text,
                egui::FontId::new(HudTypography::VALUE_SIZE, egui::FontFamily::Monospace),
                HudColors::TEXT_PRIMARY,
            );
        }

        response
    }
}

/// Section header with subtle line
pub fn section_header(ui: &mut Ui, text: &str) {
    ui.add_space(4.0);
    ui.horizontal(|ui| {
        ui.label(
            egui::RichText::new(text)
                .size(HudTypography::LABEL_SIZE)
                .color(HudColors::TEXT_LABEL)
                .strong(),
        );
        ui.add_space(4.0);
        let rect = ui.available_rect_before_wrap();
        ui.painter().line_segment(
            [
                Pos2::new(rect.left(), rect.center().y),
                Pos2::new(rect.right(), rect.center().y),
            ],
            Stroke::new(1.0, HudColors::BORDER),
        );
    });
    ui.add_space(2.0);
}

/// Key-value pair in compact format
pub fn key_value(ui: &mut Ui, key: &str, value: &str, color: Option<Color32>) {
    ui.horizontal(|ui| {
        ui.label(
            egui::RichText::new(key)
                .size(HudTypography::LABEL_SIZE)
                .color(HudColors::TEXT_LABEL),
        );
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            ui.label(
                egui::RichText::new(value)
                    .size(HudTypography::VALUE_SIZE)
                    .family(egui::FontFamily::Monospace)
                    .color(color.unwrap_or(HudColors::TEXT_PRIMARY)),
            );
        });
    });
}

/// Disease status badge
pub struct DiseaseBadge {
    name: String,
    severity_percent: f64,
}

impl DiseaseBadge {
    pub fn new(name: impl Into<String>, severity_percent: f64) -> Self {
        Self {
            name: name.into(),
            severity_percent,
        }
    }
}

impl Widget for DiseaseBadge {
    fn ui(self, ui: &mut Ui) -> Response {
        let response = ui.horizontal(|ui| {
            // Red indicator dot
            let (dot_rect, _) = ui.allocate_exact_size(Vec2::splat(10.0), Sense::hover());
            ui.painter().circle_filled(dot_rect.center(), 4.0, HudColors::CRITICAL);

            // Disease name
            ui.label(
                egui::RichText::new(&self.name)
                    .size(HudTypography::TITLE_SIZE)
                    .color(HudColors::CRITICAL)
                    .strong(),
            );

            // Severity
            ui.label(
                egui::RichText::new(format!("{:.0}%", self.severity_percent))
                    .size(HudTypography::VALUE_SIZE)
                    .family(egui::FontFamily::Monospace)
                    .color(HudColors::WARNING),
            );
        });

        response.response
    }
}
