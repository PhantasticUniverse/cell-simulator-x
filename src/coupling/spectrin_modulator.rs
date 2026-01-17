//! ATP-dependent spectrin stiffness modulation.
//!
//! Implements the reverse coupling from biochemistry to physics:
//! ATP depletion → reduced spectrin phosphorylation → increased stiffness.
//!
//! ## Biological Basis
//! Spectrin is phosphorylated by kinases (e.g., casein kinase I) which requires ATP.
//! Phosphorylation reduces spectrin-actin binding affinity, making the network
//! more compliant. Under ATP depletion:
//! - Dephosphorylation dominates
//! - Spectrin-actin junctions strengthen
//! - Membrane becomes stiffer
//!
//! This contributes to the characteristic shape changes in stored RBCs
//! (storage lesion) where ATP depletes over time.
//!
//! ## Implementation
//! The stiffness modifier is a function of ATP concentration:
//! - At normal ATP (2.0 mM): modifier = 1.0 (baseline stiffness)
//! - At zero ATP: modifier = 1.5 (50% increase in stiffness)
//! - Linear interpolation between these extremes
//!
//! ## References
//! - Manno S et al., PNAS 2002 - Spectrin phosphorylation and membrane flexibility
//! - Discher DE et al., Science 1994 - Spectrin mechanics
//! - Boal DH, Soft Matter 2012 - Red blood cell cytoskeleton

/// Modulates spectrin network stiffness based on ATP levels.
///
/// This creates reverse coupling from metabolism to mechanics:
/// low ATP → stiffer membrane.
#[derive(Debug, Clone)]
pub struct SpectrinModulator {
    /// Reference ATP concentration for normal stiffness (mM)
    /// Reference: ~2.0 mM in healthy RBCs
    pub atp_reference_mM: f64,

    /// Maximum stiffening factor at zero ATP.
    /// A value of 0.5 means 50% increase in stiffness when ATP = 0.
    /// Reference: Manno et al. 2002 showed ~40-60% increase in rigidity
    pub max_stiffening_factor: f64,

    /// Minimum ATP threshold below which maximum stiffening applies (mM)
    /// Prevents division by zero and caps extreme behavior
    pub min_atp_threshold_mM: f64,
}

impl SpectrinModulator {
    /// Create a new spectrin modulator with custom parameters.
    pub fn new(atp_reference_mM: f64, max_stiffening_factor: f64) -> Self {
        Self {
            atp_reference_mM,
            max_stiffening_factor,
            min_atp_threshold_mM: 0.1, // 0.1 mM minimum
        }
    }

    /// Calculate the stiffness modifier based on current ATP level.
    ///
    /// Returns a multiplier for spectrin stiffness:
    /// - 1.0 at normal ATP (2.0 mM)
    /// - Up to 1.5 at zero ATP (with max_stiffening_factor = 0.5)
    ///
    /// # Arguments
    /// * `atp_mM` - Current ATP concentration in mM
    ///
    /// # Formula
    /// ```text
    /// normalized = ATP / ATP_ref (clamped to [0, 1])
    /// modifier = 1.0 + max_stiffening * (1.0 - normalized)
    /// ```
    pub fn stiffness_modifier(&self, atp_mM: f64) -> f32 {
        // Normalize ATP relative to reference
        let normalized = (atp_mM / self.atp_reference_mM).clamp(0.0, 1.0);

        // Linear interpolation from 1.0 (normal) to (1 + max_stiffening) (depleted)
        let modifier = 1.0 + self.max_stiffening_factor * (1.0 - normalized);

        modifier as f32
    }

    /// Calculate modified persistence length for WLC model.
    ///
    /// Lower persistence length = stiffer spectrin.
    /// L_p' = L_p / modifier
    ///
    /// # Arguments
    /// * `base_persistence_length_um` - Normal persistence length (typically 0.020 μm = 20 nm)
    /// * `atp_mM` - Current ATP concentration
    pub fn modified_persistence_length(&self, base_persistence_length_um: f32, atp_mM: f64) -> f32 {
        let modifier = self.stiffness_modifier(atp_mM);
        // Higher modifier = lower persistence length = stiffer
        base_persistence_length_um / modifier
    }

    /// Calculate modified shear modulus for Skalak model.
    ///
    /// Higher shear modulus = stiffer membrane.
    /// G_s' = G_s × modifier
    ///
    /// # Arguments
    /// * `base_shear_modulus` - Normal shear modulus (typically 5.5 μN/m)
    /// * `atp_mM` - Current ATP concentration
    pub fn modified_shear_modulus(&self, base_shear_modulus: f32, atp_mM: f64) -> f32 {
        let modifier = self.stiffness_modifier(atp_mM);
        base_shear_modulus * modifier
    }

    /// Get descriptive status for diagnostics.
    pub fn status_description(&self, atp_mM: f64) -> &'static str {
        let modifier = self.stiffness_modifier(atp_mM);
        if modifier < 1.1 {
            "normal"
        } else if modifier < 1.25 {
            "mildly stiffened"
        } else if modifier < 1.4 {
            "moderately stiffened"
        } else {
            "severely stiffened"
        }
    }
}

impl Default for SpectrinModulator {
    fn default() -> Self {
        Self {
            atp_reference_mM: 2.0,       // Normal RBC ATP
            max_stiffening_factor: 0.5,   // 50% increase at zero ATP
            min_atp_threshold_mM: 0.1,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_stiffness_modifier_normal_atp() {
        let modulator = SpectrinModulator::default();

        // At reference ATP, modifier should be 1.0
        let modifier = modulator.stiffness_modifier(2.0);
        assert!(
            (modifier - 1.0).abs() < 0.01,
            "Modifier at normal ATP should be 1.0, got {}",
            modifier
        );
    }

    #[test]
    fn test_stiffness_modifier_half_atp() {
        let modulator = SpectrinModulator::default();

        // At half reference ATP (1.0 mM), modifier should be ~1.25
        let modifier = modulator.stiffness_modifier(1.0);
        assert!(
            (modifier - 1.25).abs() < 0.01,
            "Modifier at half ATP should be 1.25, got {}",
            modifier
        );
    }

    #[test]
    fn test_stiffness_modifier_zero_atp() {
        let modulator = SpectrinModulator::default();

        // At zero ATP, modifier should be 1.5 (1.0 + 0.5 max stiffening)
        let modifier = modulator.stiffness_modifier(0.0);
        assert!(
            (modifier - 1.5).abs() < 0.01,
            "Modifier at zero ATP should be 1.5, got {}",
            modifier
        );
    }

    #[test]
    fn test_stiffness_modifier_high_atp() {
        let modulator = SpectrinModulator::default();

        // At higher-than-normal ATP, modifier should still be 1.0 (clamped)
        let modifier = modulator.stiffness_modifier(3.0);
        assert!(
            (modifier - 1.0).abs() < 0.01,
            "Modifier at high ATP should be clamped to 1.0, got {}",
            modifier
        );
    }

    #[test]
    fn test_modified_persistence_length() {
        let modulator = SpectrinModulator::default();
        let base_lp = 0.020; // 20 nm

        // At normal ATP, persistence length unchanged
        let lp_normal = modulator.modified_persistence_length(base_lp, 2.0);
        assert!((lp_normal - base_lp).abs() < 0.001);

        // At zero ATP, persistence length reduced by 1.5x
        let lp_depleted = modulator.modified_persistence_length(base_lp, 0.0);
        let expected = base_lp / 1.5;
        assert!(
            (lp_depleted - expected).abs() < 0.001,
            "Expected {}, got {}",
            expected,
            lp_depleted
        );
    }

    #[test]
    fn test_modified_shear_modulus() {
        let modulator = SpectrinModulator::default();
        let base_gs = 5.5; // μN/m

        // At normal ATP, shear modulus unchanged
        let gs_normal = modulator.modified_shear_modulus(base_gs, 2.0);
        assert!((gs_normal - base_gs).abs() < 0.01);

        // At zero ATP, shear modulus increased by 1.5x
        let gs_depleted = modulator.modified_shear_modulus(base_gs, 0.0);
        let expected = base_gs * 1.5;
        assert!(
            (gs_depleted - expected).abs() < 0.01,
            "Expected {}, got {}",
            expected,
            gs_depleted
        );
    }

    #[test]
    fn test_monotonic_increase() {
        let modulator = SpectrinModulator::default();

        // As ATP decreases, modifier should increase
        let m1 = modulator.stiffness_modifier(2.0);
        let m2 = modulator.stiffness_modifier(1.5);
        let m3 = modulator.stiffness_modifier(1.0);
        let m4 = modulator.stiffness_modifier(0.5);
        let m5 = modulator.stiffness_modifier(0.0);

        assert!(m1 <= m2);
        assert!(m2 <= m3);
        assert!(m3 <= m4);
        assert!(m4 <= m5);
    }
}
