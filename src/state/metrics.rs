//! Unified simulation metrics for HUD display and export.
//!
//! This module provides a single data structure that aggregates all displayable
//! metrics from the various solvers into a format suitable for rendering in the
//! HUD overlay and exporting to CSV/JSON.

use serde::{Deserialize, Serialize};

/// Status indicator for metabolite values relative to physiological range
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum MetaboliteStatus {
    /// Value within normal physiological range (green)
    #[default]
    Normal,
    /// Value approaching boundaries (yellow)
    Warning,
    /// Value outside safe range (red)
    Critical,
}

impl MetaboliteStatus {
    /// Determine status based on value and normal range
    pub fn from_value(value: f64, normal_min: f64, normal_max: f64) -> Self {
        // Add 20% margin for warning zone
        let margin = (normal_max - normal_min) * 0.2;

        if value < normal_min - margin || value > normal_max + margin {
            MetaboliteStatus::Critical
        } else if value < normal_min || value > normal_max {
            MetaboliteStatus::Warning
        } else {
            MetaboliteStatus::Normal
        }
    }
}

/// Simulation mode indicator
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum SimulationMode {
    #[default]
    Normal,
    Stressed,
    Disease,
}

/// Disease indicator for HUD display
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DiseaseIndicator {
    /// Disease name (e.g., "Storage Lesion")
    pub name: String,
    /// Disease parameter (e.g., "Day 21" or "12 mM glucose")
    pub parameter: String,
    /// Severity percentage (0-100)
    pub severity_percent: f64,
    /// Affected parameters for display
    pub affected_params: Vec<(String, f64)>,
}

/// Unified metrics struct for HUD display and export
///
/// This aggregates all relevant simulation data into a single structure
/// that can be easily rendered in the UI or serialized to JSON/CSV.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimulationMetrics {
    // === Timing ===
    /// Current simulation time in seconds
    pub simulation_time_sec: f64,
    /// Frames per second (rendering performance)
    pub fps: f32,
    /// Physics steps per second
    pub steps_per_second: f32,
    /// Total physics steps executed
    pub total_steps: u64,

    // === Mode ===
    /// Current simulation mode
    pub mode: SimulationMode,
    /// Whether physics simulation is running
    pub physics_running: bool,
    /// Number of physics substeps per frame
    pub physics_substeps: usize,

    // === Core Metabolites ===
    /// ATP concentration (mM) - target: 1.5-2.0 mM
    pub atp_mM: f64,
    /// ATP status indicator
    pub atp_status: MetaboliteStatus,

    /// ADP concentration (mM)
    pub adp_mM: f64,

    /// 2,3-DPG concentration (mM) - target: 4.5-5.5 mM
    pub dpg_2_3_mM: f64,
    /// 2,3-DPG status indicator
    pub dpg_status: MetaboliteStatus,

    /// Glucose concentration (mM) - target: 4.0-6.0 mM
    pub glucose_mM: f64,
    /// Glucose status indicator
    pub glucose_status: MetaboliteStatus,

    /// Lactate concentration (mM)
    pub lactate_mM: f64,

    // === Oxygen Transport ===
    /// Oxygen saturation (0.0-1.0 fraction)
    pub o2_saturation: f64,
    /// P50 value (mmHg)
    pub p50_mmHg: f64,
    /// Current pO2 (mmHg)
    pub po2_mmHg: f64,

    // === pH ===
    /// Intracellular pH - target: 7.2
    pub ph: f64,
    /// pH status
    pub ph_status: MetaboliteStatus,

    // === Redox State ===
    /// NADPH/NADP+ ratio - target: 10-20
    pub nadph_nadp_ratio: f64,
    /// Redox ratio status
    pub redox_status: MetaboliteStatus,

    /// GSH/GSSG ratio - target: 100-400
    pub gsh_gssg_ratio: f64,
    /// GSH status
    pub gsh_status: MetaboliteStatus,

    /// H2O2 concentration (uM) - target: <5 uM
    pub h2o2_uM: f64,
    /// H2O2 status
    pub h2o2_status: MetaboliteStatus,

    // === Ion Homeostasis ===
    /// Cytosolic Na+ (mM) - target: 8-12 mM
    pub na_plus_mM: f64,
    /// Na+ status
    pub na_status: MetaboliteStatus,

    /// Cytosolic K+ (mM) - target: 135-145 mM
    pub k_plus_mM: f64,
    /// K+ status
    pub k_status: MetaboliteStatus,

    /// Cytosolic Ca2+ (nM) - target: ~100 nM
    pub ca_nM: f64,

    // === Mechanics ===
    /// Membrane tension (pN/nm)
    pub membrane_tension_pN_per_nm: f64,
    /// Stiffness modifier from ATP depletion (1.0 = normal)
    pub stiffness_modifier: f64,
    /// Total energy (pJ)
    pub total_energy_pJ: f64,
    /// Max velocity (um/s)
    pub max_velocity_um_per_sec: f32,

    // === Disease State ===
    /// Active disease indicator (if any)
    pub disease: Option<DiseaseIndicator>,
}

impl Default for SimulationMetrics {
    fn default() -> Self {
        Self {
            // Timing
            simulation_time_sec: 0.0,
            fps: 0.0,
            steps_per_second: 0.0,
            total_steps: 0,

            // Mode
            mode: SimulationMode::Normal,
            physics_running: false,
            physics_substeps: 1,

            // Core metabolites (physiological defaults)
            atp_mM: 1.8,
            atp_status: MetaboliteStatus::Normal,
            adp_mM: 0.25,
            dpg_2_3_mM: 5.0,
            dpg_status: MetaboliteStatus::Normal,
            glucose_mM: 5.0,
            glucose_status: MetaboliteStatus::Normal,
            lactate_mM: 1.0,

            // Oxygen transport
            o2_saturation: 0.97,
            p50_mmHg: 26.8,
            po2_mmHg: 100.0,

            // pH
            ph: 7.2,
            ph_status: MetaboliteStatus::Normal,

            // Redox
            nadph_nadp_ratio: 15.0,
            redox_status: MetaboliteStatus::Normal,
            gsh_gssg_ratio: 200.0,
            gsh_status: MetaboliteStatus::Normal,
            h2o2_uM: 0.5,
            h2o2_status: MetaboliteStatus::Normal,

            // Ions
            na_plus_mM: 10.0,
            na_status: MetaboliteStatus::Normal,
            k_plus_mM: 140.0,
            k_status: MetaboliteStatus::Normal,
            ca_nM: 100.0,

            // Mechanics
            membrane_tension_pN_per_nm: 0.0,
            stiffness_modifier: 1.0,
            total_energy_pJ: 0.0,
            max_velocity_um_per_sec: 0.0,

            // Disease
            disease: None,
        }
    }
}

impl SimulationMetrics {
    /// Create new metrics with default physiological values
    pub fn new() -> Self {
        Self::default()
    }

    /// Update status indicators based on current values
    pub fn update_status(&mut self) {
        // ATP: target 1.5-2.0 mM
        self.atp_status = MetaboliteStatus::from_value(self.atp_mM, 1.5, 2.0);

        // 2,3-DPG: target 4.5-5.5 mM
        self.dpg_status = MetaboliteStatus::from_value(self.dpg_2_3_mM, 4.5, 5.5);

        // Glucose: target 4.0-6.0 mM
        self.glucose_status = MetaboliteStatus::from_value(self.glucose_mM, 4.0, 6.0);

        // pH: target 7.15-7.25
        self.ph_status = MetaboliteStatus::from_value(self.ph, 7.15, 7.25);

        // NADPH/NADP+: target 10-20
        self.redox_status = MetaboliteStatus::from_value(self.nadph_nadp_ratio, 10.0, 20.0);

        // GSH/GSSG: target 100-400
        self.gsh_status = MetaboliteStatus::from_value(self.gsh_gssg_ratio, 100.0, 400.0);

        // H2O2: should be <5 uM
        if self.h2o2_uM > 10.0 {
            self.h2o2_status = MetaboliteStatus::Critical;
        } else if self.h2o2_uM > 5.0 {
            self.h2o2_status = MetaboliteStatus::Warning;
        } else {
            self.h2o2_status = MetaboliteStatus::Normal;
        }

        // Na+: target 8-12 mM
        self.na_status = MetaboliteStatus::from_value(self.na_plus_mM, 8.0, 12.0);

        // K+: target 135-145 mM
        self.k_status = MetaboliteStatus::from_value(self.k_plus_mM, 135.0, 145.0);

        // Determine simulation mode
        self.mode = if self.disease.is_some() {
            SimulationMode::Disease
        } else if self.atp_status == MetaboliteStatus::Critical
            || self.h2o2_status == MetaboliteStatus::Critical
            || self.gsh_status == MetaboliteStatus::Critical
        {
            SimulationMode::Stressed
        } else {
            SimulationMode::Normal
        };
    }

    /// Update metrics from FullyIntegratedDiagnostics
    pub fn update_from_diagnostics(
        &mut self,
        diag: &crate::biochemistry::FullyIntegratedDiagnostics,
    ) {
        self.atp_mM = diag.atp_mM;
        self.adp_mM = diag.adp_mM;
        self.dpg_2_3_mM = diag.dpg_2_3_mM;
        self.glucose_mM = diag.glucose_mM;
        self.lactate_mM = diag.lactate_mM;
        self.ph = diag.ph;
        self.o2_saturation = diag.saturation;
        self.p50_mmHg = diag.p50_mmHg;
        self.po2_mmHg = diag.po2_mmHg;
        self.nadph_nadp_ratio = diag.nadph_nadp_ratio;
        self.gsh_gssg_ratio = diag.gsh_gssg_ratio;
        self.h2o2_uM = diag.h2o2_uM;
        self.na_plus_mM = diag.na_cytosolic_mM;
        self.k_plus_mM = diag.k_cytosolic_mM;
        self.ca_nM = diag.ca_cytosolic_uM * 1000.0; // Convert uM to nM
        self.membrane_tension_pN_per_nm = diag.membrane_tension_pN_per_nm;

        self.update_status();
    }

    /// Update from CoupledDiagnostics (includes mechanics)
    pub fn update_from_coupled_diagnostics(
        &mut self,
        diag: &crate::coupling::CoupledDiagnostics,
    ) {
        self.atp_mM = diag.atp_mM;
        self.adp_mM = diag.adp_mM;
        self.ca_nM = diag.ca_cytosolic_uM * 1000.0;
        self.membrane_tension_pN_per_nm = diag.membrane_tension_pN_per_nm;
        self.stiffness_modifier = diag.stiffness_modifier as f64;
        self.simulation_time_sec = diag.time_sec;
        self.h2o2_uM = diag.h2o2_uM;

        self.update_status();
    }
}
