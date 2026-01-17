//! HUD panel visibility state management.

/// Panel visibility toggles
#[derive(Debug, Clone)]
pub struct HudState {
    /// Show simulation info panel (top-left)
    pub show_simulation_panel: bool,
    /// Show metabolites panel (top-right)
    pub show_metabolites_panel: bool,
    /// Show disease panel (bottom-right, only when disease active)
    pub show_disease_panel: bool,
    /// Show help overlay (center, toggle with H key)
    pub show_help: bool,
    /// Show export menu (toggle with E key)
    pub show_export_menu: bool,
    /// Show FPS counter
    pub show_fps: bool,
    /// HUD enabled at all
    pub hud_enabled: bool,
}

impl Default for HudState {
    fn default() -> Self {
        Self {
            show_simulation_panel: true,
            show_metabolites_panel: true,
            show_disease_panel: true,
            show_help: false,
            show_export_menu: false,
            show_fps: true,
            hud_enabled: true,
        }
    }
}

impl HudState {
    /// Create new HUD state with all panels visible
    pub fn new() -> Self {
        Self::default()
    }

    /// Toggle HUD visibility entirely
    pub fn toggle_hud(&mut self) {
        self.hud_enabled = !self.hud_enabled;
    }

    /// Toggle help overlay
    pub fn toggle_help(&mut self) {
        self.show_help = !self.show_help;
        // Close export menu when showing help
        if self.show_help {
            self.show_export_menu = false;
        }
    }

    /// Toggle export menu
    pub fn toggle_export_menu(&mut self) {
        self.show_export_menu = !self.show_export_menu;
        // Close help when showing export
        if self.show_export_menu {
            self.show_help = false;
        }
    }

    /// Toggle metabolites panel
    pub fn toggle_metabolites(&mut self) {
        self.show_metabolites_panel = !self.show_metabolites_panel;
    }

    /// Toggle disease panel
    pub fn toggle_disease(&mut self) {
        self.show_disease_panel = !self.show_disease_panel;
    }

    /// Toggle simulation panel
    pub fn toggle_simulation(&mut self) {
        self.show_simulation_panel = !self.show_simulation_panel;
    }
}
