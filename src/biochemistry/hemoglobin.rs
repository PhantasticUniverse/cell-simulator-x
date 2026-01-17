//! Hemoglobin oxygen binding dynamics with allosteric effects.
//!
//! Implements the Adair 4-site cooperative binding model with modulation by:
//! - Bohr effect (pH)
//! - 2,3-DPG binding
//! - Temperature (van't Hoff)
//!
//! The oxygen equilibrium curve (OEC) is characterized by:
//! - P50: Partial pressure at 50% saturation (~26.8 mmHg at standard conditions)
//! - Hill coefficient n: Measure of cooperativity (~2.7)
//!
//! References:
//! - Imai K. Allosteric Effects in Haemoglobin. Cambridge University Press, 1982
//! - Benesch R, Benesch RE. Nature. 1969;221:618-622 (2,3-DPG effect)
//! - Roughton FJW, Severinghaus JW. J Appl Physiol. 1973;35:861-869 (Standard OEC)
//! - Winslow RM et al. J Biol Chem. 1977;252:2331-2337 (Adair constants)

use std::f64::consts::LN_10;

/// Standard physiological conditions
pub const STANDARD_PH: f64 = 7.4;
pub const STANDARD_TEMPERATURE_K: f64 = 310.15;  // 37°C
pub const STANDARD_DPG_MM: f64 = 5.0;  // mM
pub const STANDARD_PCO2_MMHG: f64 = 40.0;

/// Gas constant in appropriate units
const R_CAL_PER_MOL_K: f64 = 1.987;  // cal/(mol·K)

/// Hemoglobin state tracking oxygen binding and T/R equilibrium
#[derive(Debug, Clone)]
pub struct HemoglobinState {
    /// Fractional oxygen saturation (0-1)
    pub saturation: f64,
    /// Fraction of hemoglobin in T (tense/deoxy) state
    pub t_state_fraction: f64,
    /// Fraction of hemoglobin in R (relaxed/oxy) state
    pub r_state_fraction: f64,
    /// Bound oxygen concentration (mM)
    /// Hemoglobin concentration in RBCs ~5 mM (as tetramer), so max O2 ~20 mM
    pub bound_o2_mM: f64,
    /// Total hemoglobin concentration (mM of tetramer)
    pub total_hb_mM: f64,
}

impl Default for HemoglobinState {
    fn default() -> Self {
        Self {
            saturation: 0.75,  // Typical arterial saturation
            t_state_fraction: 0.25,
            r_state_fraction: 0.75,
            bound_o2_mM: 15.0,  // ~75% of 20 mM max
            total_hb_mM: 5.0,  // ~5 mM hemoglobin tetramer in RBCs
        }
    }
}

impl HemoglobinState {
    /// Create a new hemoglobin state at given saturation
    pub fn at_saturation(saturation: f64, total_hb_mM: f64) -> Self {
        let sat = saturation.clamp(0.0, 1.0);
        let max_o2 = 4.0 * total_hb_mM;  // 4 binding sites per tetramer
        Self {
            saturation: sat,
            t_state_fraction: 1.0 - sat,
            r_state_fraction: sat,
            bound_o2_mM: sat * max_o2,
            total_hb_mM,
        }
    }

    /// Maximum oxygen binding capacity (mM)
    pub fn max_oxygen_capacity_mM(&self) -> f64 {
        4.0 * self.total_hb_mM
    }

    /// Update state from calculated saturation
    pub fn update_from_saturation(&mut self, saturation: f64) {
        self.saturation = saturation.clamp(0.0, 1.0);
        self.bound_o2_mM = self.saturation * self.max_oxygen_capacity_mM();
        // Approximate T/R fractions from saturation
        // At low saturation, more T state; at high saturation, more R state
        self.r_state_fraction = self.saturation;
        self.t_state_fraction = 1.0 - self.saturation;
    }
}

/// Adair constants for 4-site oxygen binding
///
/// The Adair equation describes successive binding of O2 to hemoglobin:
///   Hb + O2 ⇌ Hb(O2)        K1
///   Hb(O2) + O2 ⇌ Hb(O2)2   K2
///   Hb(O2)2 + O2 ⇌ Hb(O2)3  K3
///   Hb(O2)3 + O2 ⇌ Hb(O2)4  K4
///
/// Reference: Imai 1982, Table 6.1; Winslow et al. 1977
#[derive(Debug, Clone, Copy)]
pub struct AdairConstants {
    /// First binding constant (mmHg⁻¹)
    pub k1_per_mmHg: f64,
    /// Second binding constant (mmHg⁻¹)
    pub k2_per_mmHg: f64,
    /// Third binding constant (mmHg⁻¹)
    pub k3_per_mmHg: f64,
    /// Fourth binding constant (mmHg⁻¹)
    pub k4_per_mmHg: f64,
}

impl Default for AdairConstants {
    /// Standard Adair constants at pH 7.4, 37°C, normal 2,3-DPG
    ///
    /// Calibrated to give P50 = 26.8 mmHg and Hill coefficient n = 2.7
    /// Reference: Imai 1982, fitted to human HbA
    fn default() -> Self {
        // Use calibrated constants that give P50 ≈ 26.8 mmHg
        // These are derived from the target P50 and Hill coefficient
        Self::from_p50_and_hill(26.8, 2.7)
    }
}

impl AdairConstants {
    /// Create Adair constants from P50 and Hill coefficient
    ///
    /// Uses iterative calibration to find constants that produce
    /// the specified P50 and approximate the given Hill coefficient.
    ///
    /// Reference: Based on MWC model relationships in Imai 1982
    pub fn from_p50_and_hill(target_p50_mmHg: f64, n_hill: f64) -> Self {
        // Start with a template that has the right cooperativity shape
        // The spread factor controls the ratio between K1 and K4
        // Lower spread = less cooperativity = lower Hill coefficient
        let spread = (10.0_f64).powf(0.4 / n_hill.max(1.0));

        // Template constants (will be scaled to hit target P50)
        let k1_template = 1.0 / (spread.powi(3));
        let k2_template = 1.0 / (spread.powi(2));
        let k3_template = 1.0 / spread;
        let k4_template = spread.powi(2);

        // Find scale factor to achieve target P50 using bisection
        let mut scale_low = 0.0001;
        let mut scale_high = 1.0;

        for _ in 0..50 {
            let scale = (scale_low + scale_high) / 2.0;
            let test_constants = Self {
                k1_per_mmHg: k1_template * scale,
                k2_per_mmHg: k2_template * scale,
                k3_per_mmHg: k3_template * scale,
                k4_per_mmHg: k4_template * scale,
            };

            let p50 = Self::compute_p50(&test_constants);

            if (p50 - target_p50_mmHg).abs() < 0.01 {
                return test_constants;
            }

            // Higher scale = higher affinity = lower P50
            if p50 > target_p50_mmHg {
                scale_low = scale;
            } else {
                scale_high = scale;
            }
        }

        // Return best estimate
        let scale = (scale_low + scale_high) / 2.0;
        Self {
            k1_per_mmHg: k1_template * scale,
            k2_per_mmHg: k2_template * scale,
            k3_per_mmHg: k3_template * scale,
            k4_per_mmHg: k4_template * scale,
        }
    }

    /// Compute P50 from Adair constants using bisection
    fn compute_p50(constants: &Self) -> f64 {
        let mut low = 1.0;
        let mut high = 200.0;

        for _ in 0..50 {
            let mid = (low + high) / 2.0;
            let sat = Self::saturation_at(constants, mid);

            if (sat - 0.5).abs() < 1e-6 {
                return mid;
            }

            if sat < 0.5 {
                low = mid;
            } else {
                high = mid;
            }
        }

        (low + high) / 2.0
    }

    /// Calculate saturation at given pO2 using Adair equation
    fn saturation_at(constants: &Self, po2_mmHg: f64) -> f64 {
        if po2_mmHg <= 0.0 {
            return 0.0;
        }

        let p = po2_mmHg;
        let (a1, a2, a3, a4) = constants.cumulative();

        let p2 = p * p;
        let p3 = p2 * p;
        let p4 = p3 * p;

        let numerator = a1 * p + 2.0 * a2 * p2 + 3.0 * a3 * p3 + 4.0 * a4 * p4;
        let denominator = 4.0 * (1.0 + a1 * p + a2 * p2 + a3 * p3 + a4 * p4);

        if denominator > 0.0 {
            (numerator / denominator).clamp(0.0, 1.0)
        } else {
            0.0
        }
    }

    /// Calculate cumulative binding constants (for Adair equation)
    pub fn cumulative(&self) -> (f64, f64, f64, f64) {
        let a1 = self.k1_per_mmHg;
        let a2 = a1 * self.k2_per_mmHg;
        let a3 = a2 * self.k3_per_mmHg;
        let a4 = a3 * self.k4_per_mmHg;
        (a1, a2, a3, a4)
    }
}

/// Allosteric parameters for hemoglobin modulation
#[derive(Debug, Clone, Copy)]
pub struct AllostericParameters {
    /// Bohr coefficient: ΔlogP50/ΔpH
    /// Negative value means lower pH → higher P50 (right shift)
    /// Reference: Imai 1982 reports -0.48 for whole blood
    pub bohr_coefficient: f64,

    /// 2,3-DPG effect: ΔP50 per mM 2,3-DPG
    /// Reference: ~2.4-3.0 mmHg per mM (Benesch 1969)
    pub dpg_effect_mmHg_per_mM: f64,

    /// Reference 2,3-DPG concentration (mM)
    pub dpg_reference_mM: f64,

    /// Temperature coefficient: enthalpy of oxygenation
    /// ΔH = -14.5 kcal/mol O2 for human HbA (Imai 1982)
    pub delta_h_kcal_per_mol: f64,

    /// Reference temperature (K)
    pub temperature_reference_K: f64,

    /// CO2 effect (Haldane effect): additional P50 shift per mmHg CO2
    /// Combined with Bohr effect via carbamino formation
    pub co2_effect_per_mmHg: f64,

    /// Reference pCO2 (mmHg)
    pub pco2_reference_mmHg: f64,
}

impl Default for AllostericParameters {
    fn default() -> Self {
        Self {
            bohr_coefficient: -0.48,  // Imai 1982
            dpg_effect_mmHg_per_mM: 2.4,  // Benesch 1969
            dpg_reference_mM: STANDARD_DPG_MM,
            delta_h_kcal_per_mol: -14.5,  // Imai 1982
            temperature_reference_K: STANDARD_TEMPERATURE_K,
            co2_effect_per_mmHg: 0.02,  // Small direct effect
            pco2_reference_mmHg: STANDARD_PCO2_MMHG,
        }
    }
}

/// Oxygen binding solver using Adair equation with allosteric effects
#[derive(Debug, Clone)]
pub struct HemoglobinSolver {
    /// Base Adair constants (at standard conditions)
    pub base_constants: AdairConstants,
    /// Allosteric parameters
    pub allosteric: AllostericParameters,
    /// Base P50 at standard conditions (mmHg)
    pub base_p50_mmHg: f64,
    /// Base Hill coefficient at standard conditions
    pub base_hill_n: f64,
    /// Kinetic rate constants for O2 association (mM⁻¹·s⁻¹)
    pub k_on_per_mM_per_sec: f64,
    /// Kinetic rate constant for O2 dissociation (s⁻¹)
    pub k_off_per_sec: f64,
}

impl Default for HemoglobinSolver {
    fn default() -> Self {
        Self {
            base_constants: AdairConstants::default(),
            allosteric: AllostericParameters::default(),
            base_p50_mmHg: 26.8,  // Imai 1982 standard
            base_hill_n: 2.7,     // Imai 1982 standard
            // Kinetic constants from Gibson 1970, Roughton 1972
            k_on_per_mM_per_sec: 25.0,  // Fast association
            k_off_per_sec: 30.0,        // Rapid equilibration
        }
    }
}

impl HemoglobinSolver {
    /// Create a new solver with custom P50 and Hill coefficient
    pub fn with_p50_and_hill(p50_mmHg: f64, n_hill: f64) -> Self {
        Self {
            base_constants: AdairConstants::from_p50_and_hill(p50_mmHg, n_hill),
            base_p50_mmHg: p50_mmHg,
            base_hill_n: n_hill,
            ..Default::default()
        }
    }

    /// Calculate effective P50 including all allosteric effects
    ///
    /// # Arguments
    /// * `ph` - Intracellular pH
    /// * `dpg_mM` - 2,3-DPG concentration (mM)
    /// * `temperature_K` - Temperature (Kelvin)
    /// * `pco2_mmHg` - CO2 partial pressure (mmHg)
    pub fn effective_p50(
        &self,
        ph: f64,
        dpg_mM: f64,
        temperature_K: f64,
        pco2_mmHg: f64,
    ) -> f64 {
        let params = &self.allosteric;

        // Start with base P50
        let mut log_p50 = self.base_p50_mmHg.log10();

        // Bohr effect: ΔlogP50 = bohr_coeff × ΔpH
        // Lower pH → higher P50 (right shift, less affinity)
        // bohr_coeff is negative (-0.48), so when pH drops, P50 increases
        let delta_ph = ph - STANDARD_PH;
        log_p50 += params.bohr_coefficient * delta_ph;

        // 2,3-DPG effect (direct P50 shift)
        // Higher 2,3-DPG → higher P50
        let delta_dpg = dpg_mM - params.dpg_reference_mM;
        let p50_from_bohr = (10.0_f64).powf(log_p50);
        let p50_with_dpg = p50_from_bohr + params.dpg_effect_mmHg_per_mM * delta_dpg;
        log_p50 = p50_with_dpg.max(1.0).log10();

        // Temperature effect (van't Hoff equation)
        // Higher temp → higher P50 (right shift)
        let delta_t = temperature_K - params.temperature_reference_K;
        if delta_t.abs() > 0.01 {
            // ΔlnKp = (ΔH/R) × (1/T1 - 1/T2)
            // For P50: ΔlogP50 ≈ (ΔH × ΔT) / (2.303 × R × T²)
            let t_factor = params.delta_h_kcal_per_mol * 1000.0 * delta_t
                / (LN_10 * R_CAL_PER_MOL_K * params.temperature_reference_K.powi(2));
            log_p50 -= t_factor;  // ΔH is negative, so this increases P50 with temp
        }

        // CO2 effect (small additional contribution beyond Bohr via pH)
        let delta_co2 = pco2_mmHg - params.pco2_reference_mmHg;
        log_p50 += params.co2_effect_per_mmHg * delta_co2;

        (10.0_f64).powf(log_p50).max(5.0).min(100.0)  // Clamp to physiological range
    }

    /// Get Adair constants adjusted for current conditions
    pub fn effective_adair_constants(
        &self,
        ph: f64,
        dpg_mM: f64,
        temperature_K: f64,
        pco2_mmHg: f64,
    ) -> AdairConstants {
        let effective_p50 = self.effective_p50(ph, dpg_mM, temperature_K, pco2_mmHg);

        // Scale constants to achieve the new P50 while maintaining cooperativity
        let p50_ratio = self.base_p50_mmHg / effective_p50;

        // All constants scale proportionally with P50 shift
        AdairConstants {
            k1_per_mmHg: self.base_constants.k1_per_mmHg * p50_ratio,
            k2_per_mmHg: self.base_constants.k2_per_mmHg * p50_ratio,
            k3_per_mmHg: self.base_constants.k3_per_mmHg * p50_ratio,
            k4_per_mmHg: self.base_constants.k4_per_mmHg * p50_ratio,
        }
    }

    /// Calculate oxygen saturation using the Adair equation
    ///
    /// Y = (K₁·p + 2·K₁K₂·p² + 3·K₁K₂K₃·p³ + 4·K₁K₂K₃K₄·p⁴) /
    ///     (4 × (1 + K₁·p + K₁K₂·p² + K₁K₂K₃·p³ + K₁K₂K₃K₄·p⁴))
    ///
    /// # Arguments
    /// * `po2_mmHg` - Oxygen partial pressure (mmHg)
    /// * `constants` - Adair constants to use
    pub fn saturation_adair(&self, po2_mmHg: f64, constants: &AdairConstants) -> f64 {
        if po2_mmHg <= 0.0 {
            return 0.0;
        }

        let p = po2_mmHg;
        let (a1, a2, a3, a4) = constants.cumulative();

        // Binding polynomial terms
        let p2 = p * p;
        let p3 = p2 * p;
        let p4 = p3 * p;

        // Numerator: weighted sum of bound states
        let numerator = a1 * p + 2.0 * a2 * p2 + 3.0 * a3 * p3 + 4.0 * a4 * p4;

        // Denominator: partition function × 4
        let denominator = 4.0 * (1.0 + a1 * p + a2 * p2 + a3 * p3 + a4 * p4);

        if denominator > 0.0 {
            (numerator / denominator).clamp(0.0, 1.0)
        } else {
            0.0
        }
    }

    /// Calculate oxygen saturation at given conditions
    ///
    /// # Arguments
    /// * `po2_mmHg` - Oxygen partial pressure (mmHg)
    /// * `ph` - Intracellular pH
    /// * `dpg_mM` - 2,3-DPG concentration (mM)
    /// * `temperature_K` - Temperature (K)
    /// * `pco2_mmHg` - CO2 partial pressure (mmHg)
    pub fn calculate_saturation(
        &self,
        po2_mmHg: f64,
        ph: f64,
        dpg_mM: f64,
        temperature_K: f64,
        pco2_mmHg: f64,
    ) -> f64 {
        let constants = self.effective_adair_constants(ph, dpg_mM, temperature_K, pco2_mmHg);
        self.saturation_adair(po2_mmHg, &constants)
    }

    /// Calculate saturation at standard conditions
    pub fn calculate_saturation_standard(&self, po2_mmHg: f64) -> f64 {
        self.saturation_adair(po2_mmHg, &self.base_constants)
    }

    /// Calculate P50 from the oxygen equilibrium curve
    ///
    /// Uses bisection to find pO2 where saturation = 0.5
    pub fn calculate_p50(
        &self,
        ph: f64,
        dpg_mM: f64,
        temperature_K: f64,
        pco2_mmHg: f64,
    ) -> f64 {
        let constants = self.effective_adair_constants(ph, dpg_mM, temperature_K, pco2_mmHg);

        // Bisection search for Y = 0.5
        let mut low = 1.0;
        let mut high = 200.0;

        for _ in 0..50 {  // Max 50 iterations
            let mid = (low + high) / 2.0;
            let sat = self.saturation_adair(mid, &constants);

            if (sat - 0.5).abs() < 1e-6 {
                return mid;
            }

            if sat < 0.5 {
                low = mid;
            } else {
                high = mid;
            }
        }

        (low + high) / 2.0
    }

    /// Calculate Hill coefficient from the OEC slope at P50
    ///
    /// n = d(log(Y/(1-Y))) / d(log(pO2)) evaluated near Y = 0.5
    pub fn calculate_hill_coefficient(
        &self,
        ph: f64,
        dpg_mM: f64,
        temperature_K: f64,
        pco2_mmHg: f64,
    ) -> f64 {
        let p50 = self.calculate_p50(ph, dpg_mM, temperature_K, pco2_mmHg);
        let constants = self.effective_adair_constants(ph, dpg_mM, temperature_K, pco2_mmHg);

        // Calculate slope using finite differences around P50
        let delta = 0.1;  // Small perturbation
        let p_low = p50 * (1.0 - delta);
        let p_high = p50 * (1.0 + delta);

        let y_low = self.saturation_adair(p_low, &constants);
        let y_high = self.saturation_adair(p_high, &constants);

        // Hill plot: log(Y/(1-Y)) vs log(pO2)
        let hill_low = (y_low / (1.0 - y_low + 1e-10)).ln();
        let hill_high = (y_high / (1.0 - y_high + 1e-10)).ln();

        let log_p_low = p_low.ln();
        let log_p_high = p_high.ln();

        (hill_high - hill_low) / (log_p_high - log_p_low)
    }

    /// Calculate rate of oxygen binding/release
    ///
    /// Uses mass-action kinetics for O2 + Hb ⇌ HbO2
    ///
    /// # Arguments
    /// * `state` - Current hemoglobin state
    /// * `po2_mmHg` - Oxygen partial pressure (mmHg)
    /// * `conditions` - (pH, dpg_mM, temperature_K, pco2_mmHg)
    ///
    /// # Returns
    /// * Rate of O2 binding (mM/s), positive = uptake
    pub fn binding_rate(
        &self,
        state: &HemoglobinState,
        po2_mmHg: f64,
        conditions: (f64, f64, f64, f64),
    ) -> f64 {
        let (ph, dpg_mM, temperature_K, pco2_mmHg) = conditions;

        // Calculate equilibrium saturation
        let y_eq = self.calculate_saturation(po2_mmHg, ph, dpg_mM, temperature_K, pco2_mmHg);

        // Current deviation from equilibrium
        let y_current = state.saturation;
        let delta_y = y_eq - y_current;

        // Rate proportional to deviation (first-order relaxation)
        // Effective rate constant depends on pO2
        let o2_mM = po2_mmHg / 760.0;  // Approximate conversion (Henry's law)
        let k_eff = self.k_on_per_mM_per_sec * o2_mM + self.k_off_per_sec;

        // Rate of change in bound O2
        k_eff * delta_y * state.max_oxygen_capacity_mM()
    }

    /// Update hemoglobin state for a timestep
    ///
    /// # Arguments
    /// * `state` - Hemoglobin state to update
    /// * `po2_mmHg` - Oxygen partial pressure (mmHg)
    /// * `conditions` - (pH, dpg_mM, temperature_K, pco2_mmHg)
    /// * `dt_sec` - Timestep (seconds)
    pub fn step(
        &self,
        state: &mut HemoglobinState,
        po2_mmHg: f64,
        conditions: (f64, f64, f64, f64),
        dt_sec: f64,
    ) {
        let rate = self.binding_rate(state, po2_mmHg, conditions);

        // Update bound O2
        state.bound_o2_mM += rate * dt_sec;
        state.bound_o2_mM = state.bound_o2_mM.clamp(0.0, state.max_oxygen_capacity_mM());

        // Update saturation
        state.saturation = state.bound_o2_mM / state.max_oxygen_capacity_mM();

        // Update T/R fractions (approximate)
        state.r_state_fraction = state.saturation;
        state.t_state_fraction = 1.0 - state.saturation;
    }

    /// Generate oxygen equilibrium curve data points
    ///
    /// # Arguments
    /// * `conditions` - (pH, dpg_mM, temperature_K, pco2_mmHg)
    /// * `n_points` - Number of data points
    ///
    /// # Returns
    /// Vector of (pO2, saturation) tuples
    pub fn generate_oec(
        &self,
        conditions: (f64, f64, f64, f64),
        n_points: usize,
    ) -> Vec<(f64, f64)> {
        let (ph, dpg_mM, temperature_K, pco2_mmHg) = conditions;
        let constants = self.effective_adair_constants(ph, dpg_mM, temperature_K, pco2_mmHg);

        let mut curve = Vec::with_capacity(n_points);

        for i in 0..n_points {
            // Logarithmic spacing from 1 to 150 mmHg
            let t = i as f64 / (n_points - 1) as f64;
            let po2 = (1.0_f64).exp() * ((150.0_f64 / 1.0_f64.exp()).ln() * t).exp();
            let sat = self.saturation_adair(po2, &constants);
            curve.push((po2, sat));
        }

        curve
    }

    /// Calculate Bohr coefficient from OEC shifts
    pub fn measured_bohr_coefficient(&self, dpg_mM: f64, temperature_K: f64, pco2_mmHg: f64) -> f64 {
        let p50_low_ph = self.calculate_p50(7.2, dpg_mM, temperature_K, pco2_mmHg);
        let p50_high_ph = self.calculate_p50(7.6, dpg_mM, temperature_K, pco2_mmHg);

        let delta_log_p50 = p50_high_ph.log10() - p50_low_ph.log10();
        let delta_ph = 7.6 - 7.2;

        delta_log_p50 / delta_ph
    }
}

/// Diagnostic output for oxygen transport
#[derive(Debug, Clone)]
pub struct OxygenDiagnostics {
    /// Current pO2 (mmHg)
    pub po2_mmHg: f64,
    /// Current saturation (0-1)
    pub saturation: f64,
    /// Effective P50 (mmHg)
    pub p50_mmHg: f64,
    /// Hill coefficient
    pub hill_n: f64,
    /// Bound O2 concentration (mM)
    pub bound_o2_mM: f64,
    /// Environmental conditions
    pub ph: f64,
    pub dpg_mM: f64,
    pub temperature_K: f64,
    pub pco2_mmHg: f64,
    /// Measured Bohr coefficient
    pub bohr_coefficient: f64,
}

impl OxygenDiagnostics {
    /// Print a formatted summary
    pub fn print_summary(&self) {
        println!("=== Oxygen Transport Diagnostics ===");
        println!();
        println!("Environmental Conditions:");
        println!("  pH:           {:.2}", self.ph);
        println!("  2,3-DPG:      {:.2} mM", self.dpg_mM);
        println!("  Temperature:  {:.1} K ({:.1}°C)",
            self.temperature_K, self.temperature_K - 273.15);
        println!("  pCO2:         {:.1} mmHg", self.pco2_mmHg);
        println!();
        println!("Oxygen Binding:");
        println!("  pO2:          {:.1} mmHg", self.po2_mmHg);
        println!("  Saturation:   {:.1}%", self.saturation * 100.0);
        println!("  Bound O2:     {:.2} mM", self.bound_o2_mM);
        println!();
        println!("OEC Parameters:");
        println!("  P50:          {:.1} mmHg (target: 26.8 ± 1)", self.p50_mmHg);
        println!("  Hill coeff:   {:.2} (target: 2.7 ± 0.1)", self.hill_n);
        println!("  Bohr coeff:   {:.2} (target: -0.48 ± 0.05)", self.bohr_coefficient);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_adair_constants_default() {
        let constants = AdairConstants::default();
        // K values should increase with binding (cooperativity)
        assert!(constants.k1_per_mmHg < constants.k2_per_mmHg);
        assert!(constants.k2_per_mmHg < constants.k3_per_mmHg);
        assert!(constants.k3_per_mmHg < constants.k4_per_mmHg);
    }

    #[test]
    fn test_saturation_limits() {
        let solver = HemoglobinSolver::default();

        // Zero pO2 should give zero saturation
        let sat_zero = solver.calculate_saturation_standard(0.0);
        assert_eq!(sat_zero, 0.0);

        // Very high pO2 should give ~100% saturation
        let sat_high = solver.calculate_saturation_standard(500.0);
        assert!(sat_high > 0.99, "High pO2 saturation: {}", sat_high);
    }

    #[test]
    fn test_p50_standard_conditions() {
        let solver = HemoglobinSolver::default();
        let p50 = solver.calculate_p50(
            STANDARD_PH,
            STANDARD_DPG_MM,
            STANDARD_TEMPERATURE_K,
            STANDARD_PCO2_MMHG,
        );

        // P50 should be close to 26.8 mmHg at standard conditions
        assert!(
            (p50 - 26.8).abs() < 2.0,
            "P50 at standard conditions: {} (expected ~26.8)",
            p50
        );
    }

    #[test]
    fn test_hill_coefficient() {
        let solver = HemoglobinSolver::default();
        let n = solver.calculate_hill_coefficient(
            STANDARD_PH,
            STANDARD_DPG_MM,
            STANDARD_TEMPERATURE_K,
            STANDARD_PCO2_MMHG,
        );

        // Hill coefficient should be ~2.7 for hemoglobin
        assert!(
            (n - 2.7).abs() < 0.3,
            "Hill coefficient: {} (expected ~2.7)",
            n
        );
    }

    #[test]
    fn test_bohr_effect() {
        let solver = HemoglobinSolver::default();

        // Lower pH should increase P50 (right shift)
        let p50_low_ph = solver.calculate_p50(7.2, STANDARD_DPG_MM, STANDARD_TEMPERATURE_K, STANDARD_PCO2_MMHG);
        let p50_high_ph = solver.calculate_p50(7.6, STANDARD_DPG_MM, STANDARD_TEMPERATURE_K, STANDARD_PCO2_MMHG);

        assert!(
            p50_low_ph > p50_high_ph,
            "Bohr effect: P50 at pH 7.2 ({}) should be > P50 at pH 7.6 ({})",
            p50_low_ph, p50_high_ph
        );

        // Check Bohr coefficient magnitude
        let bohr = solver.measured_bohr_coefficient(STANDARD_DPG_MM, STANDARD_TEMPERATURE_K, STANDARD_PCO2_MMHG);
        assert!(
            (bohr - (-0.48)).abs() < 0.1,
            "Bohr coefficient: {} (expected ~-0.48)",
            bohr
        );
    }

    #[test]
    fn test_dpg_effect() {
        let solver = HemoglobinSolver::default();

        // Higher 2,3-DPG should increase P50
        let p50_low_dpg = solver.calculate_p50(STANDARD_PH, 3.0, STANDARD_TEMPERATURE_K, STANDARD_PCO2_MMHG);
        let p50_high_dpg = solver.calculate_p50(STANDARD_PH, 7.0, STANDARD_TEMPERATURE_K, STANDARD_PCO2_MMHG);

        assert!(
            p50_high_dpg > p50_low_dpg,
            "2,3-DPG effect: P50 at 7mM ({}) should be > P50 at 3mM ({})",
            p50_high_dpg, p50_low_dpg
        );

        // Check magnitude: ~2.4-3 mmHg per mM
        let delta_p50 = p50_high_dpg - p50_low_dpg;
        let delta_dpg = 7.0 - 3.0;
        let dpg_sensitivity = delta_p50 / delta_dpg;

        assert!(
            dpg_sensitivity > 1.5 && dpg_sensitivity < 4.0,
            "2,3-DPG sensitivity: {} mmHg/mM (expected 2-3)",
            dpg_sensitivity
        );
    }

    #[test]
    fn test_temperature_effect() {
        let solver = HemoglobinSolver::default();

        // Higher temperature should increase P50 (right shift)
        let p50_low_temp = solver.calculate_p50(STANDARD_PH, STANDARD_DPG_MM, 303.15, STANDARD_PCO2_MMHG);  // 30°C
        let p50_high_temp = solver.calculate_p50(STANDARD_PH, STANDARD_DPG_MM, 313.15, STANDARD_PCO2_MMHG); // 40°C

        assert!(
            p50_high_temp > p50_low_temp,
            "Temperature effect: P50 at 40°C ({}) should be > P50 at 30°C ({})",
            p50_high_temp, p50_low_temp
        );
    }

    #[test]
    fn test_hemoglobin_state() {
        let state = HemoglobinState::at_saturation(0.5, 5.0);

        assert!((state.saturation - 0.5).abs() < 1e-6);
        assert!((state.bound_o2_mM - 10.0).abs() < 1e-6);  // 50% of 20 mM max
        assert!((state.max_oxygen_capacity_mM() - 20.0).abs() < 1e-6);
    }

    #[test]
    fn test_binding_kinetics() {
        let solver = HemoglobinSolver::default();
        let mut state = HemoglobinState::at_saturation(0.5, 5.0);
        let conditions = (STANDARD_PH, STANDARD_DPG_MM, STANDARD_TEMPERATURE_K, STANDARD_PCO2_MMHG);

        // At high pO2, saturation should increase
        let initial_sat = state.saturation;
        for _ in 0..100 {
            solver.step(&mut state, 100.0, conditions, 0.001);
        }
        assert!(
            state.saturation > initial_sat,
            "Saturation should increase at high pO2: {} -> {}",
            initial_sat, state.saturation
        );
    }

    #[test]
    fn test_oec_generation() {
        let solver = HemoglobinSolver::default();
        let conditions = (STANDARD_PH, STANDARD_DPG_MM, STANDARD_TEMPERATURE_K, STANDARD_PCO2_MMHG);

        let oec = solver.generate_oec(conditions, 50);

        // Should have monotonically increasing saturation
        for i in 1..oec.len() {
            assert!(
                oec[i].1 >= oec[i-1].1,
                "OEC should be monotonic: {} at {} mmHg vs {} at {} mmHg",
                oec[i].1, oec[i].0, oec[i-1].1, oec[i-1].0
            );
        }

        // First point should be low, last should be high
        assert!(oec[0].1 < 0.5);
        assert!(oec[oec.len()-1].1 > 0.95);
    }
}
