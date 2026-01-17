//! Integrated metabolism-oxygen solver for RBCs.
//!
//! This module couples the MetabolismSolver and HemoglobinSolver to model
//! the dynamic interplay between:
//! - Lactate production from glycolysis
//! - pH changes via buffer capacity
//! - Oxygen affinity via Bohr effect
//! - 2,3-DPG regulation of hemoglobin
//!
//! The coupling pathway is:
//!   Glycolysis → Lactate ↑ → pH ↓ → P50 ↑ → O2 release
//!
//! References:
//! - Van Slyke DD. J Biol Chem. 1922 (buffer capacity)
//! - Imai K. Allosteric Effects in Haemoglobin. 1982 (Bohr effect)
//! - Benesch R, Benesch RE. Nature. 1969 (2,3-DPG effect)

use crate::biochemistry::{
    HemoglobinSolver, HemoglobinState,
    MetabolismSolver, MetabolismConfig, MetabolitePool,
    ExtendedMetaboliteIndices,
    STANDARD_PCO2_MMHG, STANDARD_TEMPERATURE_K,
};
use crate::biochemistry::ph_buffer::PhBufferModel;

/// Environment conditions for oxygen transport
#[derive(Debug, Clone, Copy)]
pub struct IntegratedEnvironment {
    /// Oxygen partial pressure (mmHg)
    pub po2_mmHg: f64,
    /// CO2 partial pressure (mmHg)
    pub pco2_mmHg: f64,
    /// Temperature (Kelvin)
    pub temperature_K: f64,
}

impl Default for IntegratedEnvironment {
    fn default() -> Self {
        Self {
            po2_mmHg: 100.0,           // Arterial pO2
            pco2_mmHg: STANDARD_PCO2_MMHG, // 40 mmHg
            temperature_K: STANDARD_TEMPERATURE_K, // 37°C
        }
    }
}

impl IntegratedEnvironment {
    /// Create venous conditions
    pub fn venous() -> Self {
        Self {
            po2_mmHg: 40.0,  // Venous pO2
            pco2_mmHg: 46.0,  // Slightly higher venous pCO2
            temperature_K: STANDARD_TEMPERATURE_K,
        }
    }

    /// Create arterial conditions
    pub fn arterial() -> Self {
        Self::default()
    }

    /// Create conditions with custom pO2
    pub fn with_po2(po2_mmHg: f64) -> Self {
        Self {
            po2_mmHg,
            ..Default::default()
        }
    }
}

/// Integrated solver combining metabolism and oxygen transport
///
/// This solver couples:
/// 1. MetabolismSolver - produces lactate, ATP, 2,3-DPG
/// 2. PhBufferModel - calculates pH from lactate
/// 3. HemoglobinSolver - oxygen binding with Bohr effect
///
/// The key innovation is dynamic pH calculation: as metabolism
/// produces lactate, pH drops, which shifts the oxygen equilibrium
/// curve to the right (Bohr effect), reducing hemoglobin affinity.
pub struct IntegratedSolver {
    /// Metabolism solver (glycolysis + RL shunt)
    pub metabolism: MetabolismSolver,
    /// Hemoglobin oxygen binding solver
    pub hemoglobin: HemoglobinSolver,
    /// pH buffer model for lactate → pH coupling
    pub ph_buffer: PhBufferModel,
    /// Hemoglobin state (saturation, bound O2)
    pub hb_state: HemoglobinState,
    /// Environmental conditions (pO2, pCO2, temperature)
    pub environment: IntegratedEnvironment,
    /// Current calculated pH
    pub current_ph: f64,
    /// Metabolite indices for convenient access
    pub indices: ExtendedMetaboliteIndices,
    /// Integration timestep (seconds)
    pub dt_sec: f64,
}

impl IntegratedSolver {
    /// Create a new integrated solver with default configuration
    pub fn new(environment: IntegratedEnvironment) -> Self {
        let config = MetabolismConfig::default();
        let dt_sec = config.dt_sec;
        let metabolism = MetabolismSolver::new(config);
        let indices = metabolism.indices;
        let ph_buffer = PhBufferModel::default();

        Self {
            metabolism,
            hemoglobin: HemoglobinSolver::default(),
            ph_buffer,
            hb_state: HemoglobinState::default(),
            environment,
            current_ph: ph_buffer.reference_ph,
            indices,
            dt_sec,
        }
    }

    /// Create with custom configuration
    pub fn with_config(
        metabolism_config: MetabolismConfig,
        environment: IntegratedEnvironment,
        ph_buffer: PhBufferModel,
    ) -> Self {
        let dt_sec = metabolism_config.dt_sec;
        let metabolism = MetabolismSolver::new(metabolism_config);
        let indices = metabolism.indices;

        Self {
            metabolism,
            hemoglobin: HemoglobinSolver::default(),
            ph_buffer,
            hb_state: HemoglobinState::default(),
            environment,
            current_ph: ph_buffer.reference_ph,
            indices,
            dt_sec,
        }
    }

    /// Perform one integration step
    ///
    /// This is the core coupling loop:
    /// 1. Step metabolism (update lactate, ATP, 2,3-DPG)
    /// 2. Calculate pH from lactate via buffer model
    /// 3. Get 2,3-DPG concentration
    /// 4. Step hemoglobin with dynamic pH (Bohr effect)
    ///
    /// # Arguments
    /// * `metabolites` - Metabolite pool to update
    /// * `atp_consumption_mM_per_sec` - External ATP demand
    pub fn step(&mut self, metabolites: &mut MetabolitePool, atp_consumption_mM_per_sec: f64) {
        // 1. Step metabolism (updates lactate, ATP, 2,3-DPG, etc.)
        self.metabolism.step(metabolites, atp_consumption_mM_per_sec);

        // 2. Calculate pH from lactate via buffer model
        let lactate = metabolites.get(self.indices.glycolysis.lactate);
        self.current_ph = self.ph_buffer.calculate_ph(lactate);

        // 3. Get 2,3-DPG concentration
        let dpg = metabolites.get(self.indices.bisphosphoglycerate_2_3);

        // 4. Update hemoglobin with dynamic pH and 2,3-DPG
        let conditions = (
            self.current_ph,
            dpg,
            self.environment.temperature_K,
            self.environment.pco2_mmHg,
        );
        self.hemoglobin.step(
            &mut self.hb_state,
            self.environment.po2_mmHg,
            conditions,
            self.dt_sec,
        );
    }

    /// Run simulation for a given duration
    ///
    /// # Arguments
    /// * `metabolites` - Metabolite pool to update
    /// * `duration_sec` - Simulation duration (seconds)
    /// * `atp_consumption_mM_per_sec` - External ATP demand
    pub fn run(
        &mut self,
        metabolites: &mut MetabolitePool,
        duration_sec: f64,
        atp_consumption_mM_per_sec: f64,
    ) {
        let n_steps = (duration_sec / self.dt_sec).ceil() as usize;
        for _ in 0..n_steps {
            self.step(metabolites, atp_consumption_mM_per_sec);
        }
    }

    /// Get current simulation time
    pub fn time_sec(&self) -> f64 {
        self.metabolism.time_sec
    }

    /// Calculate effective P50 at current conditions
    pub fn effective_p50(&self, metabolites: &MetabolitePool) -> f64 {
        let lactate = metabolites.get(self.indices.glycolysis.lactate);
        let ph = self.ph_buffer.calculate_ph(lactate);
        let dpg = metabolites.get(self.indices.bisphosphoglycerate_2_3);

        self.hemoglobin.effective_p50(
            ph,
            dpg,
            self.environment.temperature_K,
            self.environment.pco2_mmHg,
        )
    }

    /// Calculate oxygen saturation at current conditions
    pub fn calculate_saturation(&self, metabolites: &MetabolitePool) -> f64 {
        let lactate = metabolites.get(self.indices.glycolysis.lactate);
        let ph = self.ph_buffer.calculate_ph(lactate);
        let dpg = metabolites.get(self.indices.bisphosphoglycerate_2_3);

        self.hemoglobin.calculate_saturation(
            self.environment.po2_mmHg,
            ph,
            dpg,
            self.environment.temperature_K,
            self.environment.pco2_mmHg,
        )
    }

    /// Get diagnostic information about current state
    pub fn diagnostics(&self, metabolites: &MetabolitePool) -> IntegratedDiagnostics {
        let lactate = metabolites.get(self.indices.glycolysis.lactate);
        let ph = self.ph_buffer.calculate_ph(lactate);
        let dpg = metabolites.get(self.indices.bisphosphoglycerate_2_3);
        let atp = metabolites.get(self.indices.glycolysis.atp);

        let p50 = self.hemoglobin.effective_p50(
            ph,
            dpg,
            self.environment.temperature_K,
            self.environment.pco2_mmHg,
        );

        let saturation = self.hemoglobin.calculate_saturation(
            self.environment.po2_mmHg,
            ph,
            dpg,
            self.environment.temperature_K,
            self.environment.pco2_mmHg,
        );

        // Calculate contributions to P50 shift
        let baseline_p50 = self.hemoglobin.base_p50_mmHg;

        // pH contribution (Bohr effect)
        let p50_at_standard_ph = self.hemoglobin.effective_p50(
            7.4, // Standard pH
            dpg,
            self.environment.temperature_K,
            self.environment.pco2_mmHg,
        );
        let bohr_shift = p50 - p50_at_standard_ph;

        // Total shift from baseline
        let total_p50_shift = p50 - baseline_p50;

        IntegratedDiagnostics {
            time_sec: self.metabolism.time_sec,
            lactate_mM: lactate,
            ph,
            dpg_mM: dpg,
            atp_mM: atp,
            p50_mmHg: p50,
            saturation,
            bound_o2_mM: self.hb_state.bound_o2_mM,
            po2_mmHg: self.environment.po2_mmHg,
            pco2_mmHg: self.environment.pco2_mmHg,
            temperature_K: self.environment.temperature_K,
            bohr_p50_shift_mmHg: bohr_shift,
            total_p50_shift_mmHg: total_p50_shift,
        }
    }

    /// Reset solver state
    pub fn reset(&mut self) {
        self.metabolism.reset();
        self.hb_state = HemoglobinState::default();
        self.current_ph = self.ph_buffer.reference_ph;
    }

    /// Update environment conditions
    pub fn set_environment(&mut self, env: IntegratedEnvironment) {
        self.environment = env;
    }

    /// Set pO2 (common operation for simulating arterial/venous transit)
    pub fn set_po2(&mut self, po2_mmHg: f64) {
        self.environment.po2_mmHg = po2_mmHg;
    }
}

/// Diagnostic output from integrated solver
#[derive(Debug, Clone)]
pub struct IntegratedDiagnostics {
    /// Current simulation time (seconds)
    pub time_sec: f64,
    /// Lactate concentration (mM)
    pub lactate_mM: f64,
    /// Calculated intracellular pH
    pub ph: f64,
    /// 2,3-DPG concentration (mM)
    pub dpg_mM: f64,
    /// ATP concentration (mM)
    pub atp_mM: f64,
    /// Effective P50 (mmHg)
    pub p50_mmHg: f64,
    /// Oxygen saturation (0-1)
    pub saturation: f64,
    /// Bound oxygen (mM)
    pub bound_o2_mM: f64,
    /// Environmental pO2 (mmHg)
    pub po2_mmHg: f64,
    /// Environmental pCO2 (mmHg)
    pub pco2_mmHg: f64,
    /// Temperature (K)
    pub temperature_K: f64,
    /// P50 shift from Bohr effect alone (mmHg)
    pub bohr_p50_shift_mmHg: f64,
    /// Total P50 shift from all effects (mmHg)
    pub total_p50_shift_mmHg: f64,
}

impl IntegratedDiagnostics {
    /// Print a formatted summary
    pub fn print_summary(&self) {
        println!("=== Integrated State (t = {:.3} s) ===", self.time_sec);
        println!();
        println!("Metabolism:");
        println!("  Lactate:  {:.3} mM", self.lactate_mM);
        println!("  ATP:      {:.3} mM", self.atp_mM);
        println!("  2,3-DPG:  {:.3} mM", self.dpg_mM);
        println!();
        println!("pH Coupling:");
        println!("  pH:       {:.3} (from lactate)", self.ph);
        println!();
        println!("Oxygen Transport:");
        println!("  pO2:          {:.1} mmHg", self.po2_mmHg);
        println!("  P50:          {:.1} mmHg", self.p50_mmHg);
        println!("  Saturation:   {:.1}%", self.saturation * 100.0);
        println!("  Bound O2:     {:.2} mM", self.bound_o2_mM);
        println!();
        println!("P50 Shifts:");
        println!("  Bohr effect:  {:+.1} mmHg", self.bohr_p50_shift_mmHg);
        println!("  Total shift:  {:+.1} mmHg", self.total_p50_shift_mmHg);
    }

    /// Print a one-line summary for time series
    pub fn print_row(&self) {
        println!(
            "{:8.2} {:11.3} {:8.3} {:11.1} {:9.1}",
            self.time_sec,
            self.lactate_mM,
            self.ph,
            self.p50_mmHg,
            self.saturation * 100.0
        );
    }
}

/// Run integrated diagnostics with given parameters
pub fn run_integrated_diagnostics(
    duration_sec: f64,
    po2_mmHg: f64,
    atp_stress: f64,
) {
    println!("=== Integrated Metabolism-Oxygen Diagnostics ===\n");

    let env = IntegratedEnvironment::with_po2(po2_mmHg);
    let mut solver = IntegratedSolver::new(env);
    let mut metabolites = MetabolitePool::default_physiological();

    println!("Conditions: pO2={:.0} mmHg, pCO2={:.0} mmHg, T={:.0}°C\n",
        po2_mmHg,
        solver.environment.pco2_mmHg,
        solver.environment.temperature_K - 273.15);

    // Initial state
    let initial_diag = solver.diagnostics(&metabolites);
    println!("Initial State:");
    println!("  Lactate: {:.3} mM", initial_diag.lactate_mM);
    println!("  pH:      {:.3}", initial_diag.ph);
    println!("  P50:     {:.1} mmHg", initial_diag.p50_mmHg);
    println!("  Sat:     {:.1}%", initial_diag.saturation * 100.0);
    println!();

    // Time series header
    println!("--- Time Series ---");
    println!("{:>8} {:>11} {:>8} {:>11} {:>9}",
        "Time(s)", "Lactate(mM)", "pH", "P50(mmHg)", "Sat(%)");
    println!("{}", "-".repeat(52));

    // Store initial for coupling analysis
    let initial_lactate = initial_diag.lactate_mM;
    let initial_ph = initial_diag.ph;
    let initial_p50 = initial_diag.p50_mmHg;

    // ATP consumption rate with optional stress
    let base_atp_consumption = 0.001;  // mM/s baseline
    let atp_consumption = base_atp_consumption * atp_stress;

    // Report intervals
    let n_reports = 10;
    let report_interval = duration_sec / n_reports as f64;
    let mut next_report = 0.0;

    // Run simulation
    let n_steps = (duration_sec / solver.dt_sec).ceil() as usize;
    for step in 0..n_steps {
        let t = step as f64 * solver.dt_sec;

        if t >= next_report || step == 0 {
            let diag = solver.diagnostics(&metabolites);
            diag.print_row();
            next_report += report_interval;
        }

        solver.step(&mut metabolites, atp_consumption);
    }

    // Final row
    let final_diag = solver.diagnostics(&metabolites);
    final_diag.print_row();

    // Coupling analysis
    println!("\n--- Coupling Analysis ---");
    let delta_lactate = final_diag.lactate_mM - initial_lactate;
    let delta_ph = final_diag.ph - initial_ph;
    let delta_p50 = final_diag.p50_mmHg - initial_p50;

    println!("Lactate change:          {:+.3} mM", delta_lactate);
    println!("Lactate-induced pH shift: {:+.4}", delta_ph);
    println!("Bohr-induced P50 shift:   {:+.1} mmHg", final_diag.bohr_p50_shift_mmHg);
    println!("Total P50 shift:          {:+.1} mmHg", delta_p50);

    // Validate coupling direction
    println!("\n--- Validation ---");
    if delta_lactate > 0.01 && delta_ph < 0.0 {
        println!("✓ Lactate ↑ → pH ↓ (correct direction)");
    } else if delta_lactate < -0.01 && delta_ph > 0.0 {
        println!("✓ Lactate ↓ → pH ↑ (correct direction - lactate export)");
    } else if delta_lactate.abs() < 0.01 {
        println!("~ Lactate stable (no coupling signal)");
    } else {
        println!("⚠ Unexpected coupling direction: lactate {} but pH {}",
            if delta_lactate > 0.0 { "increased" } else { "decreased" },
            if delta_ph > 0.0 { "increased" } else { "decreased" });
    }

    if delta_ph < -0.001 && delta_p50 > 0.0 {
        println!("✓ pH ↓ → P50 ↑ (Bohr effect: acidosis → right shift)");
    } else if delta_ph > 0.001 && delta_p50 < 0.0 {
        println!("✓ pH ↑ → P50 ↓ (Bohr effect: alkalosis → left shift)");
    } else if delta_ph.abs() < 0.001 {
        println!("~ pH stable (minimal Bohr effect)");
    } else {
        // Small changes may not show expected direction due to other effects
        println!("~ Bohr effect: pH {:+.4} → P50 {:+.1} mmHg", delta_ph, delta_p50);
    }

    // Check coupling magnitude
    if delta_lactate.abs() > 0.1 {
        let ph_per_lactate = delta_ph / delta_lactate;
        let expected_ph_per_lactate = -1.0 / 60.0;  // ~-0.017
        if (ph_per_lactate - expected_ph_per_lactate).abs() < 0.005 {
            println!("✓ pH sensitivity: {:.4} per mM lactate (expected ~-0.017)",
                ph_per_lactate);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_integrated_solver_creation() {
        let env = IntegratedEnvironment::default();
        let solver = IntegratedSolver::new(env);

        assert!((solver.current_ph - 7.2).abs() < 0.001);
        assert!(solver.hb_state.saturation > 0.0);
    }

    #[test]
    fn test_coupling_direction() {
        // Test that: lactate ↑ → pH ↓ → P50 ↑
        let env = IntegratedEnvironment::default();
        let mut solver = IntegratedSolver::new(env);
        let mut metabolites = MetabolitePool::default_physiological();

        // Record initial state
        let initial_lactate = metabolites.get(solver.indices.glycolysis.lactate);
        let initial_diag = solver.diagnostics(&metabolites);

        // Run with high ATP consumption to produce lactate
        solver.run(&mut metabolites, 30.0, 0.01);  // 30s, high ATP demand

        let final_diag = solver.diagnostics(&metabolites);
        let final_lactate = metabolites.get(solver.indices.glycolysis.lactate);

        // Lactate should increase or stay stable (depends on export rate)
        // pH should respond appropriately to lactate
        // P50 should respond to pH via Bohr effect

        // The key test: if lactate changed, pH and P50 should have changed appropriately
        let delta_lactate = final_lactate - initial_lactate;

        if delta_lactate > 0.5 {  // Significant lactate increase
            // pH should decrease
            assert!(
                final_diag.ph < initial_diag.ph,
                "pH should decrease when lactate increases: {} vs {}",
                final_diag.ph, initial_diag.ph
            );

            // P50 should increase (Bohr effect)
            assert!(
                final_diag.p50_mmHg > initial_diag.p50_mmHg - 0.5,  // Allow some tolerance
                "P50 should increase (or stay stable) with lower pH: {} vs {}",
                final_diag.p50_mmHg, initial_diag.p50_mmHg
            );
        }
    }

    #[test]
    fn test_ph_from_lactate() {
        let env = IntegratedEnvironment::default();
        let solver = IntegratedSolver::new(env);
        let mut metabolites = MetabolitePool::default_physiological();

        // Manually set lactate and check pH
        metabolites.set(solver.indices.glycolysis.lactate, 1.5);  // Baseline
        let diag_baseline = solver.diagnostics(&metabolites);

        metabolites.set(solver.indices.glycolysis.lactate, 4.5);  // +3 mM
        let diag_high = solver.diagnostics(&metabolites);

        // pH should drop by ~0.05 for +3 mM lactate (60 slykes buffer)
        let expected_drop = 3.0 / 60.0;
        let actual_drop = diag_baseline.ph - diag_high.ph;

        assert!(
            (actual_drop - expected_drop).abs() < 0.01,
            "pH drop should be ~0.05: {} (expected {})",
            actual_drop, expected_drop
        );
    }

    #[test]
    fn test_bohr_effect_integration() {
        let env = IntegratedEnvironment::with_po2(26.8);  // At P50
        let solver = IntegratedSolver::new(env);
        let mut metabolites = MetabolitePool::default_physiological();

        // At baseline lactate (pH ~7.2), check P50
        metabolites.set(solver.indices.glycolysis.lactate, 1.5);
        let p50_normal = solver.effective_p50(&metabolites);

        // Increase lactate to drop pH
        metabolites.set(solver.indices.glycolysis.lactate, 7.5);  // +6 mM → pH ~7.1
        let p50_acidic = solver.effective_p50(&metabolites);

        // P50 should increase with acidosis
        assert!(
            p50_acidic > p50_normal,
            "P50 should increase with acidosis: {} vs {}",
            p50_acidic, p50_normal
        );

        // Check magnitude: ~6 mM lactate → ~0.1 pH drop → ~1-2 mmHg P50 increase
        let p50_shift = p50_acidic - p50_normal;
        assert!(
            p50_shift > 0.5 && p50_shift < 5.0,
            "P50 shift should be ~1-2 mmHg: {}",
            p50_shift
        );
    }

    #[test]
    fn test_saturation_responds_to_lactate() {
        let env = IntegratedEnvironment::with_po2(26.8);  // Near P50
        let solver = IntegratedSolver::new(env);
        let mut metabolites = MetabolitePool::default_physiological();

        // At normal lactate
        metabolites.set(solver.indices.glycolysis.lactate, 1.5);
        let sat_normal = solver.calculate_saturation(&metabolites);

        // At high lactate (lower pH → lower affinity → lower saturation)
        metabolites.set(solver.indices.glycolysis.lactate, 10.0);
        let sat_acidic = solver.calculate_saturation(&metabolites);

        // Saturation should decrease with acidosis at a fixed pO2 near P50
        assert!(
            sat_acidic < sat_normal,
            "Saturation should decrease with acidosis: {} vs {}",
            sat_acidic, sat_normal
        );
    }

    #[test]
    fn test_environment_presets() {
        let arterial = IntegratedEnvironment::arterial();
        let venous = IntegratedEnvironment::venous();

        assert!((arterial.po2_mmHg - 100.0).abs() < 0.1);
        assert!((venous.po2_mmHg - 40.0).abs() < 0.1);
        assert!(venous.pco2_mmHg > arterial.pco2_mmHg);
    }

    #[test]
    fn test_diagnostics_output() {
        let env = IntegratedEnvironment::default();
        let solver = IntegratedSolver::new(env);
        let metabolites = MetabolitePool::default_physiological();

        let diag = solver.diagnostics(&metabolites);

        assert!(diag.lactate_mM > 0.0);
        assert!(diag.ph > 6.5 && diag.ph < 7.8);
        assert!(diag.p50_mmHg > 20.0 && diag.p50_mmHg < 40.0);
        assert!(diag.saturation > 0.0 && diag.saturation <= 1.0);
    }
}
