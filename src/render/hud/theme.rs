//! HUD theme system - "Laboratory Precision" dark theme.
//!
//! A professional scientific instrument interface with clinical precision
//! and biological color accents.

use egui::{Color32, FontFamily, FontId, Rounding, Stroke, Style, TextStyle, Visuals};

/// HUD color palette - scientific instrument aesthetic
pub struct HudColors;

impl HudColors {
    // === Background Colors ===
    /// Near-black with blue tint - main background
    pub const BACKGROUND: Color32 = Color32::from_rgb(10, 10, 15);
    /// Dark slate panel background (95% opacity)
    pub const PANEL_BG: Color32 = Color32::from_rgba_premultiplied(20, 22, 30, 242);
    /// Slightly lighter for hover states
    pub const PANEL_BG_HOVER: Color32 = Color32::from_rgb(30, 33, 42);

    // === Text Colors ===
    /// Cool white - primary text
    pub const TEXT_PRIMARY: Color32 = Color32::from_rgb(230, 235, 240);
    /// Muted gray - secondary text
    pub const TEXT_SECONDARY: Color32 = Color32::from_rgb(160, 170, 180);
    /// Dim gray - labels
    pub const TEXT_LABEL: Color32 = Color32::from_rgb(110, 120, 128);
    /// Very dim for disabled
    pub const TEXT_DISABLED: Color32 = Color32::from_rgb(80, 85, 90);

    // === Status Colors ===
    /// Scientific blue - accent color
    pub const ACCENT: Color32 = Color32::from_rgb(100, 180, 255);
    /// Emerald green - healthy/normal status
    pub const SUCCESS: Color32 = Color32::from_rgb(80, 200, 120);
    /// Amber - warning status
    pub const WARNING: Color32 = Color32::from_rgb(255, 200, 50);
    /// Alarm red - critical status
    pub const CRITICAL: Color32 = Color32::from_rgb(255, 80, 80);

    // === Biological Colors ===
    /// Oxygenated hemoglobin red
    pub const HEMOGLOBIN_OXY: Color32 = Color32::from_rgb(212, 32, 32);
    /// Deoxygenated hemoglobin (dark red)
    pub const HEMOGLOBIN_DEOXY: Color32 = Color32::from_rgb(107, 32, 32);
    /// ATP energy color (phosphorescent green)
    pub const ATP_GLOW: Color32 = Color32::from_rgb(144, 238, 144);
    /// Cellular membrane color
    pub const MEMBRANE: Color32 = Color32::from_rgb(255, 200, 150);

    // === UI Element Colors ===
    /// Progress bar background
    pub const BAR_BG: Color32 = Color32::from_rgb(30, 35, 45);
    /// Progress bar fill (normal)
    pub const BAR_FILL: Color32 = Color32::from_rgb(80, 150, 200);
    /// Border color
    pub const BORDER: Color32 = Color32::from_rgb(50, 55, 65);
    /// Border highlight
    pub const BORDER_HIGHLIGHT: Color32 = Color32::from_rgb(80, 100, 120);
}

/// Typography settings for scientific precision feel
pub struct HudTypography;

impl HudTypography {
    /// Large title size
    pub const TITLE_SIZE: f32 = 14.0;
    /// Main value display size
    pub const VALUE_SIZE: f32 = 16.0;
    /// Label size
    pub const LABEL_SIZE: f32 = 11.0;
    /// Small annotation size
    pub const SMALL_SIZE: f32 = 10.0;
    /// Unit suffix size
    pub const UNIT_SIZE: f32 = 10.0;
}

/// HUD theme configuration
pub struct HudTheme {
    /// Panel corner rounding
    pub panel_rounding: f32,
    /// Button rounding
    pub button_rounding: f32,
    /// Panel padding
    pub panel_padding: f32,
    /// Spacing between elements
    pub item_spacing: f32,
    /// Progress bar height
    pub bar_height: f32,
    /// Progress bar rounding
    pub bar_rounding: f32,
}

impl Default for HudTheme {
    fn default() -> Self {
        Self {
            panel_rounding: 6.0,
            button_rounding: 4.0,
            panel_padding: 12.0,
            item_spacing: 6.0,
            bar_height: 8.0,
            bar_rounding: 3.0,
        }
    }
}

impl HudTheme {
    /// Apply theme to egui context
    pub fn apply(&self, ctx: &egui::Context) {
        let mut style = Style::default();

        // Dark visuals
        let mut visuals = Visuals::dark();

        // Override colors
        visuals.panel_fill = HudColors::PANEL_BG;
        visuals.window_fill = HudColors::PANEL_BG;
        visuals.extreme_bg_color = HudColors::BACKGROUND;
        visuals.faint_bg_color = HudColors::PANEL_BG_HOVER;

        // Text colors
        visuals.override_text_color = Some(HudColors::TEXT_PRIMARY);

        // Widget colors
        visuals.widgets.noninteractive.bg_fill = HudColors::PANEL_BG;
        visuals.widgets.noninteractive.fg_stroke = Stroke::new(1.0, HudColors::TEXT_SECONDARY);
        visuals.widgets.noninteractive.rounding = Rounding::same(self.panel_rounding);

        visuals.widgets.inactive.bg_fill = HudColors::BAR_BG;
        visuals.widgets.inactive.fg_stroke = Stroke::new(1.0, HudColors::TEXT_SECONDARY);
        visuals.widgets.inactive.rounding = Rounding::same(self.button_rounding);

        visuals.widgets.hovered.bg_fill = HudColors::PANEL_BG_HOVER;
        visuals.widgets.hovered.fg_stroke = Stroke::new(1.0, HudColors::TEXT_PRIMARY);
        visuals.widgets.hovered.rounding = Rounding::same(self.button_rounding);

        visuals.widgets.active.bg_fill = HudColors::ACCENT;
        visuals.widgets.active.fg_stroke = Stroke::new(1.0, HudColors::TEXT_PRIMARY);

        // Selection
        visuals.selection.bg_fill = HudColors::ACCENT.gamma_multiply(0.3);
        visuals.selection.stroke = Stroke::new(1.0, HudColors::ACCENT);

        // Window stroke
        visuals.window_stroke = Stroke::new(1.0, HudColors::BORDER);
        visuals.window_rounding = Rounding::same(self.panel_rounding);

        style.visuals = visuals;

        // Spacing
        style.spacing.item_spacing = egui::vec2(self.item_spacing, self.item_spacing);
        style.spacing.window_margin = egui::Margin::same(self.panel_padding);
        style.spacing.button_padding = egui::vec2(8.0, 4.0);

        // Text styles
        style.text_styles.insert(
            TextStyle::Heading,
            FontId::new(HudTypography::TITLE_SIZE, FontFamily::Proportional),
        );
        style.text_styles.insert(
            TextStyle::Body,
            FontId::new(HudTypography::VALUE_SIZE, FontFamily::Proportional),
        );
        style.text_styles.insert(
            TextStyle::Small,
            FontId::new(HudTypography::SMALL_SIZE, FontFamily::Proportional),
        );
        style.text_styles.insert(
            TextStyle::Monospace,
            FontId::new(HudTypography::VALUE_SIZE, FontFamily::Monospace),
        );

        ctx.set_style(style);
    }

    /// Get color for metabolite status
    pub fn status_color(status: crate::state::MetaboliteStatus) -> Color32 {
        match status {
            crate::state::MetaboliteStatus::Normal => HudColors::SUCCESS,
            crate::state::MetaboliteStatus::Warning => HudColors::WARNING,
            crate::state::MetaboliteStatus::Critical => HudColors::CRITICAL,
        }
    }

    /// Get color for O2 saturation (gradient from deoxy to oxy)
    pub fn o2_saturation_color(saturation: f64) -> Color32 {
        let t = saturation.clamp(0.0, 1.0) as f32;

        // Interpolate from deoxygenated (dark red) to oxygenated (bright red)
        let r = HudColors::HEMOGLOBIN_DEOXY.r() as f32
            + t * (HudColors::HEMOGLOBIN_OXY.r() as f32 - HudColors::HEMOGLOBIN_DEOXY.r() as f32);
        let g = HudColors::HEMOGLOBIN_DEOXY.g() as f32
            + t * (HudColors::HEMOGLOBIN_OXY.g() as f32 - HudColors::HEMOGLOBIN_DEOXY.g() as f32);
        let b = HudColors::HEMOGLOBIN_DEOXY.b() as f32
            + t * (HudColors::HEMOGLOBIN_OXY.b() as f32 - HudColors::HEMOGLOBIN_DEOXY.b() as f32);

        Color32::from_rgb(r as u8, g as u8, b as u8)
    }

    /// Get bar fill color for a value within a range
    pub fn value_bar_color(value: f64, min: f64, max: f64) -> Color32 {
        let norm = ((value - min) / (max - min)).clamp(0.0, 1.0);

        // Gradient from critical (low) through warning (mid) to success (high)
        if norm < 0.3 {
            HudColors::CRITICAL
        } else if norm < 0.7 {
            HudColors::WARNING
        } else {
            HudColors::SUCCESS
        }
    }
}
