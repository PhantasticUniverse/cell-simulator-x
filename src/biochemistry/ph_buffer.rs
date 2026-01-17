//! pH buffer model for red blood cell cytoplasm.
//!
//! RBCs have substantial intracellular buffering capacity due to:
//! - Hemoglobin (major buffer, ~50% of total)
//! - 2,3-DPG and organic phosphates
//! - Bicarbonate system
//!
//! The total buffer capacity is ~60 slykes (mM/pH unit), which means
//! a 1 mM increase in acid (like lactate) causes ~0.017 pH unit decrease.
//!
//! References:
//! - Van Slyke DD. J Biol Chem. 1922;52:525-570 (buffer capacity definition)
//! - Jacobs MH, Stewart DR. J Cell Comp Physiol. 1947;30:79-103 (RBC pH)
//! - Siggaard-Andersen O. The Acid-Base Status of the Blood. 1974 (buffer systems)

/// pH buffer model for RBC cytoplasm
///
/// Uses the Van Slyke buffer capacity to calculate pH changes
/// from lactate accumulation. The buffer capacity represents
/// the amount of strong acid (in mM) required to change pH by 1 unit.
///
/// # Example
/// ```
/// use cell_simulator_x::biochemistry::ph_buffer::PhBufferModel;
///
/// let buffer = PhBufferModel::default();
///
/// // At baseline lactate (1.5 mM), pH should be 7.2
/// let ph = buffer.calculate_ph(1.5);
/// assert!((ph - 7.2).abs() < 0.001);
///
/// // If lactate rises to 4.5 mM (+3 mM), pH drops by ~0.05
/// let ph_high_lactate = buffer.calculate_ph(4.5);
/// assert!(ph_high_lactate < 7.2);
/// assert!(ph_high_lactate > 7.1);
/// ```
#[derive(Debug, Clone, Copy)]
pub struct PhBufferModel {
    /// Total buffer capacity in slykes (mM/pH unit)
    /// Reference: ~60 slykes for human RBCs (Van Slyke 1922)
    pub buffer_capacity_slykes: f64,

    /// Reference pH at baseline lactate
    /// Reference: 7.2 for RBC cytoplasm (Jacobs 1947)
    pub reference_ph: f64,

    /// Baseline lactate concentration (mM)
    /// Reference: ~1.5 mM in normal conditions
    pub baseline_lactate_mM: f64,

    /// Minimum physiological pH (clamping)
    pub min_ph: f64,

    /// Maximum physiological pH (clamping)
    pub max_ph: f64,
}

impl Default for PhBufferModel {
    fn default() -> Self {
        Self {
            buffer_capacity_slykes: 60.0,  // Van Slyke 1922
            reference_ph: 7.2,              // Jacobs 1947
            baseline_lactate_mM: 1.5,       // Normal resting level
            min_ph: 6.8,                    // Extreme acidosis limit
            max_ph: 7.6,                    // Extreme alkalosis limit
        }
    }
}

impl PhBufferModel {
    /// Create a new pH buffer model with custom parameters
    pub fn new(
        buffer_capacity_slykes: f64,
        reference_ph: f64,
        baseline_lactate_mM: f64,
    ) -> Self {
        Self {
            buffer_capacity_slykes,
            reference_ph,
            baseline_lactate_mM,
            ..Default::default()
        }
    }

    /// Calculate intracellular pH based on lactate concentration
    ///
    /// Uses the relationship: ΔpH = -ΔLactate / β_total
    /// where β_total is the buffer capacity in slykes.
    ///
    /// Lactate is a strong acid (pKa ~3.9), so it fully dissociates
    /// at physiological pH and contributes H+ ions proportionally.
    ///
    /// # Arguments
    /// * `lactate_mM` - Current lactate concentration (mM)
    ///
    /// # Returns
    /// Calculated pH, clamped to physiological range
    pub fn calculate_ph(&self, lactate_mM: f64) -> f64 {
        // Calculate change in lactate from baseline
        let delta_lactate = lactate_mM - self.baseline_lactate_mM;

        // pH change: negative because lactate is an acid
        // Higher lactate → lower pH
        let delta_ph = -delta_lactate / self.buffer_capacity_slykes;

        // Apply change to reference pH and clamp
        let ph = self.reference_ph + delta_ph;
        ph.clamp(self.min_ph, self.max_ph)
    }

    /// Calculate lactate concentration from pH
    ///
    /// Inverse of calculate_ph - useful for setting initial conditions
    /// or validation.
    ///
    /// # Arguments
    /// * `ph` - Target pH value
    ///
    /// # Returns
    /// Calculated lactate concentration (mM)
    pub fn calculate_lactate_from_ph(&self, ph: f64) -> f64 {
        let delta_ph = ph - self.reference_ph;
        let delta_lactate = -delta_ph * self.buffer_capacity_slykes;
        (self.baseline_lactate_mM + delta_lactate).max(0.0)
    }

    /// Calculate pH change per mM lactate change
    ///
    /// This is the sensitivity coefficient for lactate → pH coupling.
    /// At standard buffer capacity (~60 slykes), this is ~-0.017 pH/mM.
    pub fn ph_sensitivity_per_mM_lactate(&self) -> f64 {
        -1.0 / self.buffer_capacity_slykes
    }

    /// Validate that pH is within physiological bounds
    pub fn is_physiological(&self, ph: f64) -> bool {
        ph >= self.min_ph && ph <= self.max_ph
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_values() {
        let buffer = PhBufferModel::default();
        assert!((buffer.buffer_capacity_slykes - 60.0).abs() < 0.001);
        assert!((buffer.reference_ph - 7.2).abs() < 0.001);
        assert!((buffer.baseline_lactate_mM - 1.5).abs() < 0.001);
    }

    #[test]
    fn test_ph_at_baseline_lactate() {
        let buffer = PhBufferModel::default();
        let ph = buffer.calculate_ph(buffer.baseline_lactate_mM);
        assert!(
            (ph - buffer.reference_ph).abs() < 1e-6,
            "pH at baseline lactate should equal reference pH: {} vs {}",
            ph, buffer.reference_ph
        );
    }

    #[test]
    fn test_lactate_increase_decreases_ph() {
        let buffer = PhBufferModel::default();

        // Increase lactate by 3 mM (1.5 → 4.5)
        let ph_baseline = buffer.calculate_ph(1.5);
        let ph_high_lactate = buffer.calculate_ph(4.5);

        assert!(
            ph_high_lactate < ph_baseline,
            "Higher lactate should decrease pH: {} vs {}",
            ph_high_lactate, ph_baseline
        );

        // Expected drop: 3 mM / 60 slykes = 0.05 pH units
        let expected_drop = 3.0 / 60.0;
        let actual_drop = ph_baseline - ph_high_lactate;
        assert!(
            (actual_drop - expected_drop).abs() < 0.001,
            "pH drop should be ~0.05: {} vs expected {}",
            actual_drop, expected_drop
        );
    }

    #[test]
    fn test_ph_sensitivity() {
        let buffer = PhBufferModel::default();
        let sensitivity = buffer.ph_sensitivity_per_mM_lactate();

        // At 60 slykes: -1/60 ≈ -0.0167
        assert!(
            (sensitivity - (-1.0 / 60.0)).abs() < 0.0001,
            "Sensitivity should be ~-0.017: {}",
            sensitivity
        );
    }

    #[test]
    fn test_ph_clamping() {
        let buffer = PhBufferModel::default();

        // Very high lactate should clamp to min pH
        // Need lactate high enough to push pH below 6.8
        // delta_lactate needed: (7.2 - 6.8) * 60 = 24 mM above baseline
        let ph_extreme = buffer.calculate_ph(30.0);  // +28.5 mM from baseline
        assert!(
            (ph_extreme - buffer.min_ph).abs() < 0.001,
            "Extreme lactate should clamp to min pH: {} vs {}",
            ph_extreme, buffer.min_ph
        );

        // Very negative lactate (hypothetical) should clamp to max pH
        // Need lactate low enough to push pH above 7.6
        // delta_lactate needed: (7.2 - 7.6) * 60 = -24 mM below baseline
        // So lactate = 1.5 - 24 = -22.5 (or lower)
        let ph_low_lactate = buffer.calculate_ph(-25.0);  // Would push pH up significantly
        assert!(
            (ph_low_lactate - buffer.max_ph).abs() < 0.001,
            "Very low lactate should clamp to max pH: {} vs {}",
            ph_low_lactate, buffer.max_ph
        );
    }

    #[test]
    fn test_round_trip_lactate_ph() {
        let buffer = PhBufferModel::default();

        // Test round trip: lactate → pH → lactate
        let original_lactate = 3.0;
        let ph = buffer.calculate_ph(original_lactate);
        let recovered_lactate = buffer.calculate_lactate_from_ph(ph);

        assert!(
            (recovered_lactate - original_lactate).abs() < 0.001,
            "Round trip should preserve lactate: {} vs {}",
            recovered_lactate, original_lactate
        );
    }

    #[test]
    fn test_physiological_validation() {
        let buffer = PhBufferModel::default();

        assert!(buffer.is_physiological(7.2));
        assert!(buffer.is_physiological(7.0));
        assert!(buffer.is_physiological(7.4));
        assert!(!buffer.is_physiological(6.5));
        assert!(!buffer.is_physiological(8.0));
    }

    #[test]
    fn test_coupling_magnitude() {
        // Validate the coupling produces reasonable physiological changes
        // From the plan: 1 mM lactate ≈ 0.017 pH drop
        let buffer = PhBufferModel::default();

        let ph_1 = buffer.calculate_ph(1.5);
        let ph_2 = buffer.calculate_ph(2.5);
        let ph_drop = ph_1 - ph_2;

        // 1 mM increase should cause ~0.017 pH drop
        assert!(
            (ph_drop - 0.0167).abs() < 0.002,
            "1 mM lactate increase should cause ~0.017 pH drop: {}",
            ph_drop
        );
    }
}
