//! Rapoport-Luebering shunt for 2,3-BPG regulation.
//!
//! The Rapoport-Luebering (R-L) shunt is a bypass of glycolysis unique to
//! red blood cells that regulates 2,3-bisphosphoglycerate (2,3-BPG) levels.
//! 2,3-BPG is a critical regulator of hemoglobin oxygen affinity.
//!
//! Pathway:
//!   1,3-BPG ─(BPGM)→ 2,3-BPG ─(BPGP)→ 3-PG
//!
//! At steady state, ~20% of 1,3-BPG flows through the shunt.
//! This bypass sacrifices 1 ATP per glucose but provides oxygen affinity control.
//!
//! References:
//! - Rose ZB, Liebowitz J. J Biol Chem. 1970;245:3232-3241 (Discovery)
//! - Benesch R, Benesch RE. Nature. 1969;221:618-622 (2,3-BPG function)
//! - Mulquiney PJ, Kuchel PW. Biochem J. 1999;342:567-580 (Kinetic model)

use super::enzyme::{Enzyme, ReactionStoichiometry, michaelis_menten};
use super::MetabolitePool;
use super::glycolysis::MetaboliteIndices;

/// Extended indices including 2,3-BPG
#[derive(Debug, Clone, Copy)]
pub struct ShuntIndices {
    /// 1,3-bisphosphoglycerate (from glycolysis)
    pub bisphosphoglycerate_1_3: usize,
    /// 2,3-bisphosphoglycerate (the shunt intermediate)
    pub bisphosphoglycerate_2_3: usize,
    /// 3-phosphoglycerate (feeds back to glycolysis)
    pub phosphoglycerate_3: usize,
    /// Inorganic phosphate (released by BPGP)
    pub pi: usize,
}

impl ShuntIndices {
    /// Create shunt indices compatible with glycolysis indices
    pub fn from_glycolysis(glycolysis_indices: &MetaboliteIndices, dpg_idx: usize) -> Self {
        Self {
            bisphosphoglycerate_1_3: glycolysis_indices.bisphosphoglycerate_1_3,
            bisphosphoglycerate_2_3: dpg_idx,
            phosphoglycerate_3: glycolysis_indices.phosphoglycerate_3,
            pi: glycolysis_indices.pi,
        }
    }
}

/// Rapoport-Luebering shunt solver
///
/// Contains the two enzymes of the shunt:
/// - BPGM (bisphosphoglycerate mutase): 1,3-BPG → 2,3-BPG
/// - BPGP (bisphosphoglycerate phosphatase): 2,3-BPG → 3-PG + Pi
pub struct RapoportLueberingSolver {
    pub indices: ShuntIndices,
    pub bpgm: Bpgm,
    pub bpgp: Bpgp,
}

impl RapoportLueberingSolver {
    /// Create a new shunt solver
    ///
    /// The 2,3-BPG index must be provided as it extends the glycolysis indices
    pub fn new(glycolysis_indices: &MetaboliteIndices, dpg_index: usize) -> Self {
        let indices = ShuntIndices::from_glycolysis(glycolysis_indices, dpg_index);
        Self {
            bpgm: Bpgm::new(&indices),
            bpgp: Bpgp::new(&indices),
            indices,
        }
    }

    /// Compute derivatives for shunt reactions
    pub fn compute_derivatives(&self, metabolites: &MetabolitePool, dydt: &mut [f64]) {
        self.bpgm.apply_to_derivatives(metabolites, dydt);
        self.bpgp.apply_to_derivatives(metabolites, dydt);
    }

    /// Get reaction rates for diagnostics
    pub fn get_rates(&self, metabolites: &MetabolitePool) -> Vec<(&'static str, f64)> {
        vec![
            (self.bpgm.name(), self.bpgm.rate(metabolites)),
            (self.bpgp.name(), self.bpgp.rate(metabolites)),
        ]
    }

    /// Calculate the shunt flux ratio (fraction of 1,3-BPG going through shunt)
    ///
    /// At steady state, this should be ~0.2 (20%)
    pub fn shunt_ratio(&self, metabolites: &MetabolitePool, pgk_rate: f64) -> f64 {
        let bpgm_rate = self.bpgm.rate(metabolites);
        let total_13bpg_flux = bpgm_rate + pgk_rate;
        if total_13bpg_flux > 1e-12 {
            bpgm_rate / total_13bpg_flux
        } else {
            0.0
        }
    }
}

// ============================================================================
// Shunt Enzymes
// ============================================================================

/// Bisphosphoglycerate Mutase (BPGM) - EC 5.4.2.4
///
/// 1,3-BPG → 2,3-BPG
///
/// This enzyme converts 1,3-BPG to 2,3-BPG, diverting carbon from
/// the ATP-generating PGK reaction. The 2,3-BPG produced binds
/// deoxyhemoglobin and decreases oxygen affinity.
///
/// Reference: Rose & Liebowitz 1970, Mulquiney-Kuchel 1999
pub struct Bpgm {
    /// Maximum velocity (mM/s)
    /// Reference: ~0.045 mM/s at 37°C in RBCs
    pub vmax_mM_per_sec: f64,
    /// Km for 1,3-BPG (mM)
    /// Reference: ~0.004 mM (very high affinity)
    pub km_13bpg_mM: f64,
    /// Product inhibition constant for 2,3-BPG (mM)
    /// 2,3-BPG inhibits its own synthesis
    pub ki_23bpg_mM: f64,
    stoichiometry: ReactionStoichiometry,
    idx_13bpg: usize,
    idx_23bpg: usize,
}

impl Bpgm {
    pub fn new(indices: &ShuntIndices) -> Self {
        Self {
            vmax_mM_per_sec: 0.045,
            km_13bpg_mM: 0.004,
            ki_23bpg_mM: 2.0,  // Inhibited by high 2,3-BPG
            stoichiometry: ReactionStoichiometry::new(
                vec![(indices.bisphosphoglycerate_1_3, 1.0)],
                vec![(indices.bisphosphoglycerate_2_3, 1.0)],
            ),
            idx_13bpg: indices.bisphosphoglycerate_1_3,
            idx_23bpg: indices.bisphosphoglycerate_2_3,
        }
    }
}

impl Enzyme for Bpgm {
    fn name(&self) -> &'static str { "BPGM" }
    fn ec_number(&self) -> &'static str { "5.4.2.4" }

    fn rate(&self, metabolites: &MetabolitePool) -> f64 {
        let bpg_13 = metabolites.get(self.idx_13bpg);
        let bpg_23 = metabolites.get(self.idx_23bpg);

        // Michaelis-Menten with product inhibition
        let base_rate = michaelis_menten(self.vmax_mM_per_sec, self.km_13bpg_mM, bpg_13);

        // Product inhibition by 2,3-BPG (feedback regulation)
        base_rate / (1.0 + bpg_23 / self.ki_23bpg_mM)
    }

    fn stoichiometry(&self) -> &ReactionStoichiometry {
        &self.stoichiometry
    }
}

/// Bisphosphoglycerate Phosphatase (BPGP) - EC 3.1.3.13
///
/// 2,3-BPG + H2O → 3-PG + Pi
///
/// This enzyme hydrolyzes 2,3-BPG back to 3-phosphoglycerate,
/// which re-enters glycolysis. The rate of this enzyme determines
/// the steady-state 2,3-BPG concentration.
///
/// The enzyme is activated by Cl⁻ and inhibited by Pi and 3-PG.
///
/// Reference: Rose & Liebowitz 1970, Mulquiney-Kuchel 1999
pub struct Bpgp {
    /// Maximum velocity (mM/s)
    /// Reference: ~0.0015 mM/s (slow, maintaining high 2,3-BPG)
    pub vmax_mM_per_sec: f64,
    /// Km for 2,3-BPG (mM)
    /// Reference: ~2.0 mM (relatively high, given [2,3-BPG] ~5 mM)
    pub km_23bpg_mM: f64,
    /// Ki for Pi (competitive inhibition)
    pub ki_pi_mM: f64,
    /// Ki for 3-PG (product inhibition)
    pub ki_3pg_mM: f64,
    stoichiometry: ReactionStoichiometry,
    idx_23bpg: usize,
    idx_3pg: usize,
    idx_pi: usize,
}

impl Bpgp {
    pub fn new(indices: &ShuntIndices) -> Self {
        Self {
            vmax_mM_per_sec: 0.0015,
            km_23bpg_mM: 2.0,
            ki_pi_mM: 0.5,
            ki_3pg_mM: 1.0,
            stoichiometry: ReactionStoichiometry::new(
                vec![(indices.bisphosphoglycerate_2_3, 1.0)],
                vec![
                    (indices.phosphoglycerate_3, 1.0),
                    (indices.pi, 1.0),
                ],
            ),
            idx_23bpg: indices.bisphosphoglycerate_2_3,
            idx_3pg: indices.phosphoglycerate_3,
            idx_pi: indices.pi,
        }
    }
}

impl Enzyme for Bpgp {
    fn name(&self) -> &'static str { "BPGP" }
    fn ec_number(&self) -> &'static str { "3.1.3.13" }

    fn rate(&self, metabolites: &MetabolitePool) -> f64 {
        let bpg_23 = metabolites.get(self.idx_23bpg);
        let pi = metabolites.get(self.idx_pi);
        let pg_3 = metabolites.get(self.idx_3pg);

        // Effective Km increased by competitive inhibition from Pi
        let km_apparent = self.km_23bpg_mM * (1.0 + pi / self.ki_pi_mM);

        // Base rate with apparent Km
        let base_rate = michaelis_menten(self.vmax_mM_per_sec, km_apparent, bpg_23);

        // Product inhibition by 3-PG
        base_rate / (1.0 + pg_3 / self.ki_3pg_mM)
    }

    fn stoichiometry(&self) -> &ReactionStoichiometry {
        &self.stoichiometry
    }
}

// ============================================================================
// Regulatory Functions
// ============================================================================

/// Calculate the expected P50 shift from 2,3-BPG concentration
///
/// 2,3-BPG binds to deoxyhemoglobin and stabilizes the T-state,
/// shifting the oxygen dissociation curve to the right (increasing P50).
///
/// # Arguments
/// * `dpg_mM` - 2,3-BPG concentration in mM
///
/// # Returns
/// * Approximate P50 in mmHg
///
/// Reference: Benesch & Benesch 1969, Bunn & Forget 1986
pub fn calculate_p50_from_dpg(dpg_mM: f64) -> f64 {
    // Empirical relationship: P50 increases ~3-4 mmHg per mM 2,3-BPG
    // Normal: ~5 mM 2,3-BPG → P50 ~27 mmHg
    // Base P50 without 2,3-BPG would be ~15-17 mmHg
    const BASE_P50_MMHG: f64 = 15.0;
    const DPG_EFFECT: f64 = 2.4;  // mmHg per mM 2,3-BPG

    BASE_P50_MMHG + DPG_EFFECT * dpg_mM
}

/// Estimate 2,3-BPG steady-state concentration from pH
///
/// Acidosis (low pH) inhibits BPGM and activates BPGP, lowering 2,3-BPG.
/// Alkalosis has the opposite effect.
///
/// This approximates the Bohr effect's influence on 2,3-BPG regulation.
///
/// # Arguments
/// * `ph` - Intracellular pH
///
/// # Returns
/// * Expected steady-state 2,3-BPG in mM
pub fn estimate_dpg_from_ph(ph: f64) -> f64 {
    // Normal pH 7.2 → 5.0 mM 2,3-BPG
    // Each 0.1 pH unit change alters 2,3-BPG by ~0.5 mM
    const NORMAL_PH: f64 = 7.2;
    const NORMAL_DPG_MM: f64 = 5.0;
    const DPG_PER_PH: f64 = 5.0;  // mM per pH unit

    let dpg = NORMAL_DPG_MM + DPG_PER_PH * (ph - NORMAL_PH);
    dpg.clamp(2.0, 8.0)  // Physiological range
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::biochemistry::glycolysis::MetaboliteIndices;

    fn create_test_metabolites() -> (MetabolitePool, ShuntIndices) {
        let glyco_idx = MetaboliteIndices::default();
        let dpg_idx = 17;  // Add 2,3-BPG after glycolysis metabolites
        let shunt_idx = ShuntIndices::from_glycolysis(&glyco_idx, dpg_idx);

        let mut pool = MetabolitePool::new(18);

        // Set concentrations
        pool.set(shunt_idx.bisphosphoglycerate_1_3, 0.001);  // 1 μM
        pool.set(shunt_idx.bisphosphoglycerate_2_3, 5.0);    // 5 mM (normal)
        pool.set(shunt_idx.phosphoglycerate_3, 0.1);
        pool.set(shunt_idx.pi, 1.0);

        (pool, shunt_idx)
    }

    #[test]
    fn test_bpgm_rate() {
        let (metabolites, indices) = create_test_metabolites();
        let bpgm = Bpgm::new(&indices);

        let rate = bpgm.rate(&metabolites);
        assert!(rate > 0.0, "BPGM should have positive rate with 1,3-BPG present");
        assert!(rate < bpgm.vmax_mM_per_sec, "Rate should be below Vmax");
    }

    #[test]
    fn test_bpgm_product_inhibition() {
        let (mut metabolites, indices) = create_test_metabolites();
        let bpgm = Bpgm::new(&indices);

        // Rate at low 2,3-BPG
        metabolites.set(indices.bisphosphoglycerate_2_3, 0.5);
        let rate_low_dpg = bpgm.rate(&metabolites);

        // Rate at high 2,3-BPG
        metabolites.set(indices.bisphosphoglycerate_2_3, 10.0);
        let rate_high_dpg = bpgm.rate(&metabolites);

        assert!(
            rate_low_dpg > rate_high_dpg,
            "BPGM should be inhibited by high 2,3-BPG: {} vs {}",
            rate_low_dpg, rate_high_dpg
        );
    }

    #[test]
    fn test_bpgp_rate() {
        let (metabolites, indices) = create_test_metabolites();
        let bpgp = Bpgp::new(&indices);

        let rate = bpgp.rate(&metabolites);
        assert!(rate > 0.0, "BPGP should have positive rate with 2,3-BPG present");
    }

    #[test]
    fn test_bpgp_pi_inhibition() {
        let (mut metabolites, indices) = create_test_metabolites();
        let bpgp = Bpgp::new(&indices);

        // Rate at low Pi
        metabolites.set(indices.pi, 0.1);
        let rate_low_pi = bpgp.rate(&metabolites);

        // Rate at high Pi
        metabolites.set(indices.pi, 5.0);
        let rate_high_pi = bpgp.rate(&metabolites);

        assert!(
            rate_low_pi > rate_high_pi,
            "BPGP should be inhibited by high Pi: {} vs {}",
            rate_low_pi, rate_high_pi
        );
    }

    #[test]
    fn test_shunt_solver() {
        let glyco_idx = MetaboliteIndices::default();
        let solver = RapoportLueberingSolver::new(&glyco_idx, 17);
        let (metabolites, _) = create_test_metabolites();

        let mut dydt = vec![0.0; 18];
        solver.compute_derivatives(&metabolites, &mut dydt);

        // 2,3-BPG should have net flux (production - consumption)
        assert!(dydt[17].is_finite());
    }

    #[test]
    fn test_p50_calculation() {
        // Normal 2,3-BPG (~5 mM) should give P50 ~27 mmHg
        let p50_normal = calculate_p50_from_dpg(5.0);
        assert!(
            (p50_normal - 27.0).abs() < 3.0,
            "P50 at normal 2,3-BPG should be ~27 mmHg, got {}",
            p50_normal
        );

        // Lower 2,3-BPG should give lower P50
        let p50_low = calculate_p50_from_dpg(2.0);
        assert!(
            p50_low < p50_normal,
            "Lower 2,3-BPG should decrease P50: {} vs {}",
            p50_low, p50_normal
        );

        // Higher 2,3-BPG should give higher P50
        let p50_high = calculate_p50_from_dpg(8.0);
        assert!(
            p50_high > p50_normal,
            "Higher 2,3-BPG should increase P50: {} vs {}",
            p50_high, p50_normal
        );
    }

    #[test]
    fn test_ph_effect_on_dpg() {
        // Normal pH should give normal 2,3-BPG
        let dpg_normal = estimate_dpg_from_ph(7.2);
        assert!(
            (dpg_normal - 5.0).abs() < 0.1,
            "Normal pH should give ~5 mM 2,3-BPG, got {}",
            dpg_normal
        );

        // Acidosis should decrease 2,3-BPG
        let dpg_acidosis = estimate_dpg_from_ph(7.0);
        assert!(
            dpg_acidosis < dpg_normal,
            "Acidosis should decrease 2,3-BPG: {} vs {}",
            dpg_acidosis, dpg_normal
        );

        // Alkalosis should increase 2,3-BPG
        let dpg_alkalosis = estimate_dpg_from_ph(7.4);
        assert!(
            dpg_alkalosis > dpg_normal,
            "Alkalosis should increase 2,3-BPG: {} vs {}",
            dpg_alkalosis, dpg_normal
        );
    }
}
