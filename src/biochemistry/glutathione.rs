//! Glutathione metabolism for RBC antioxidant defense.
//!
//! The glutathione system is the primary defense against oxidative stress in RBCs.
//! GSH (reduced glutathione) detoxifies H2O2 via glutathione peroxidase, and is
//! regenerated from GSSG by glutathione reductase using NADPH from the PPP.
//!
//! Key reactions:
//! - GPx: 2 GSH + H2O2 -> GSSG + 2 H2O
//! - GR: GSSG + NADPH + H+ -> 2 GSH + NADP+
//! - gamma-GCS: Glu + Cys + ATP -> gamma-Glu-Cys + ADP + Pi (GSH synthesis)
//! - GS: gamma-Glu-Cys + Gly + ATP -> GSH + ADP + Pi
//!
//! Targets:
//! - Total GSH: 2-3 mM (Beutler 1969)
//! - GSH/GSSG ratio: 100-400 (Meister 1983)
//! - Cytosolic H2O2: <5 uM (Chance 1979)
//!
//! References:
//! - Meister A, Anderson ME. Annu Rev Biochem. 1983;52:711-760
//! - Beutler E. Blood. 1969;33:637-644
//! - Chance B et al. Methods Enzymol. 1979;105:121-126

use super::enzyme::{Enzyme, ReactionStoichiometry};
use super::MetabolitePool;
use super::pentose_phosphate::RedoxIndices;

/// Basal H2O2 production rate from hemoglobin autoxidation
///
/// Met-Hb formation produces superoxide which dismutates to H2O2.
/// Reference: Winterbourn 1990, ~1-5 uM/s production
/// Increased to 5 µM/s to balance with NADPH consumption requirements.
pub const BASAL_H2O2_PRODUCTION_MM_PER_SEC: f64 = 0.005;  // 5 uM/s

// ============================================================================
// Glutathione Cycle Enzymes
// ============================================================================

/// Glutathione Peroxidase (GPx) - EC 1.11.1.9
///
/// 2 GSH + H2O2 -> GSSG + 2 H2O
///
/// Primary H2O2 detoxification enzyme. Uses ping-pong mechanism.
/// Reference: Flohé L, Günzler WA. Methods Enzymol. 1984;105:114-121
/// Vmax = 50 mM/s, Km_GSH = 1.0 mM, Km_H2O2 = 0.001 mM
pub struct GlutathionePeroxidase {
    pub vmax_mM_per_sec: f64,
    pub km_gsh_mM: f64,
    pub km_h2o2_mM: f64,
    stoichiometry: ReactionStoichiometry,
    idx_gsh: usize,
    idx_gssg: usize,
    idx_h2o2: usize,
}

impl GlutathionePeroxidase {
    pub fn new(indices: &RedoxIndices) -> Self {
        Self {
            // Reference: Flohé 1984, Board 2013
            // Tuned to achieve H2O2 steady state of ~2-5 µM with 5 µM/s production
            // At H2O2=2µM: rate ≈ 0.02 * 0.25 = 0.005 mM/s ≈ H2O2 production
            vmax_mM_per_sec: 0.02,
            km_gsh_mM: 1.0,
            km_h2o2_mM: 0.002,  // 2 µM - at this H2O2 level, rate matches production
            stoichiometry: ReactionStoichiometry::new(
                vec![
                    (indices.gsh, 2.0),
                    (indices.h2o2, 1.0),
                ],
                vec![
                    (indices.gssg, 1.0),
                    // H2O not tracked
                ],
            ),
            idx_gsh: indices.gsh,
            idx_gssg: indices.gssg,
            idx_h2o2: indices.h2o2,
        }
    }
}

impl Enzyme for GlutathionePeroxidase {
    fn name(&self) -> &'static str { "Glutathione Peroxidase" }
    fn ec_number(&self) -> &'static str { "1.11.1.9" }

    fn rate(&self, metabolites: &MetabolitePool) -> f64 {
        let gsh = metabolites.get(self.idx_gsh);
        let h2o2 = metabolites.get(self.idx_h2o2);

        if gsh <= 0.0 || h2o2 <= 0.0 {
            return 0.0;
        }

        // Ping-pong bi-bi mechanism:
        // v = Vmax * [GSH] * [H2O2] / (Km_GSH * [H2O2] + Km_H2O2 * [GSH] + [GSH] * [H2O2])
        let numerator = self.vmax_mM_per_sec * gsh * h2o2;
        let denominator = self.km_gsh_mM * h2o2 + self.km_h2o2_mM * gsh + gsh * h2o2;

        if denominator <= 0.0 {
            return 0.0;
        }

        numerator / denominator
    }

    fn stoichiometry(&self) -> &ReactionStoichiometry {
        &self.stoichiometry
    }
}

/// Glutathione Reductase (GR) - EC 1.8.1.7
///
/// GSSG + NADPH + H+ -> 2 GSH + NADP+
///
/// Regenerates GSH from GSSG using NADPH. Critical link between
/// PPP and antioxidant defense.
/// Reference: Worthington DJ, Rosemeyer MA. Eur J Biochem. 1976;67:231-238
pub struct GlutathioneReductase {
    pub vmax_mM_per_sec: f64,
    pub km_gssg_mM: f64,
    pub km_nadph_mM: f64,
    stoichiometry: ReactionStoichiometry,
    idx_gssg: usize,
    idx_gsh: usize,
    idx_nadph: usize,
    idx_nadp: usize,
}

impl GlutathioneReductase {
    pub fn new(indices: &RedoxIndices) -> Self {
        Self {
            // Reference: Worthington 1976, Beutler 1984
            // RBC GR activity: ~15 U/g Hb = ~0.085 mM/s at 340 g/L Hb
            // Km_GSSG adjusted to allow GSSG accumulation to physiological levels
            vmax_mM_per_sec: 0.15,
            km_gssg_mM: 0.015,  // Allows GSSG to accumulate for GSH/GSSG ~200
            km_nadph_mM: 0.015,
            stoichiometry: ReactionStoichiometry::new(
                vec![
                    (indices.gssg, 1.0),
                    (indices.nadph, 1.0),
                ],
                vec![
                    (indices.gsh, 2.0),
                    (indices.nadp_plus, 1.0),
                ],
            ),
            idx_gssg: indices.gssg,
            idx_gsh: indices.gsh,
            idx_nadph: indices.nadph,
            idx_nadp: indices.nadp_plus,
        }
    }
}

impl Enzyme for GlutathioneReductase {
    fn name(&self) -> &'static str { "Glutathione Reductase" }
    fn ec_number(&self) -> &'static str { "1.8.1.7" }

    fn rate(&self, metabolites: &MetabolitePool) -> f64 {
        let gssg = metabolites.get(self.idx_gssg);
        let nadph = metabolites.get(self.idx_nadph);

        if gssg <= 0.0 || nadph <= 0.0 {
            return 0.0;
        }

        // Ordered bi-bi mechanism
        let gssg_term = gssg / (self.km_gssg_mM + gssg);
        let nadph_term = nadph / (self.km_nadph_mM + nadph);

        self.vmax_mM_per_sec * gssg_term * nadph_term
    }

    fn stoichiometry(&self) -> &ReactionStoichiometry {
        &self.stoichiometry
    }
}

/// gamma-Glutamylcysteine Synthetase (gamma-GCS) - EC 6.3.2.2
///
/// Glu + Cys + ATP -> gamma-Glu-Cys + ADP + Pi
///
/// Rate-limiting step for GSH synthesis. Strong feedback inhibition by GSH.
/// Reference: Richman PG, Meister A. J Biol Chem. 1975;250:1422-1426
pub struct GammaGlutamylcysteineSynthetase {
    pub vmax_mM_per_sec: f64,
    pub km_glu_mM: f64,
    pub km_cys_mM: f64,
    pub km_atp_mM: f64,
    pub ki_gsh_mM: f64,  // Feedback inhibition by GSH
    stoichiometry: ReactionStoichiometry,
    idx_glu: usize,
    idx_cys: usize,
    idx_atp: usize,
    idx_adp: usize,
    idx_gamma_glu_cys: usize,
    idx_gsh: usize,
}

impl GammaGlutamylcysteineSynthetase {
    pub fn new(indices: &RedoxIndices) -> Self {
        Self {
            // Reference: Richman 1975, Meister 1983
            vmax_mM_per_sec: 0.01,  // Rate-limiting, slow
            km_glu_mM: 1.8,
            km_cys_mM: 0.3,
            km_atp_mM: 0.4,
            ki_gsh_mM: 2.3,  // Feedback inhibition
            stoichiometry: ReactionStoichiometry::new(
                vec![
                    (indices.glutamate, 1.0),
                    (indices.cysteine, 1.0),
                    (indices.atp, 1.0),
                ],
                vec![
                    (indices.gamma_glu_cys, 1.0),
                    (indices.adp, 1.0),
                    // Pi not tracked
                ],
            ),
            idx_glu: indices.glutamate,
            idx_cys: indices.cysteine,
            idx_atp: indices.atp,
            idx_adp: indices.adp,
            idx_gamma_glu_cys: indices.gamma_glu_cys,
            idx_gsh: indices.gsh,
        }
    }
}

impl Enzyme for GammaGlutamylcysteineSynthetase {
    fn name(&self) -> &'static str { "gamma-GCS" }
    fn ec_number(&self) -> &'static str { "6.3.2.2" }

    fn rate(&self, metabolites: &MetabolitePool) -> f64 {
        let glu = metabolites.get(self.idx_glu);
        let cys = metabolites.get(self.idx_cys);
        let atp = metabolites.get(self.idx_atp);
        let gsh = metabolites.get(self.idx_gsh);

        if glu <= 0.0 || cys <= 0.0 || atp <= 0.0 {
            return 0.0;
        }

        // Tri-substrate kinetics with GSH feedback inhibition
        let glu_term = glu / (self.km_glu_mM + glu);
        let cys_term = cys / (self.km_cys_mM + cys);
        let atp_term = atp / (self.km_atp_mM + atp);

        // Feedback inhibition by GSH (non-competitive)
        let inhibition = 1.0 / (1.0 + gsh / self.ki_gsh_mM);

        self.vmax_mM_per_sec * glu_term * cys_term * atp_term * inhibition
    }

    fn stoichiometry(&self) -> &ReactionStoichiometry {
        &self.stoichiometry
    }
}

/// Glutathione Synthetase (GS) - EC 6.3.2.3
///
/// gamma-Glu-Cys + Gly + ATP -> GSH + ADP + Pi
///
/// Final step in GSH synthesis.
/// Reference: Meister A. Methods Enzymol. 1985;113:393-399
pub struct GlutathioneSynthetase {
    pub vmax_mM_per_sec: f64,
    pub km_gamma_glu_cys_mM: f64,
    pub km_gly_mM: f64,
    pub km_atp_mM: f64,
    stoichiometry: ReactionStoichiometry,
    idx_gamma_glu_cys: usize,
    idx_gly: usize,
    idx_atp: usize,
    idx_adp: usize,
    idx_gsh: usize,
}

impl GlutathioneSynthetase {
    pub fn new(indices: &RedoxIndices) -> Self {
        Self {
            // Reference: Meister 1985
            vmax_mM_per_sec: 0.02,
            km_gamma_glu_cys_mM: 0.15,
            km_gly_mM: 0.4,
            km_atp_mM: 0.3,
            stoichiometry: ReactionStoichiometry::new(
                vec![
                    (indices.gamma_glu_cys, 1.0),
                    (indices.glycine, 1.0),
                    (indices.atp, 1.0),
                ],
                vec![
                    (indices.gsh, 1.0),
                    (indices.adp, 1.0),
                    // Pi not tracked
                ],
            ),
            idx_gamma_glu_cys: indices.gamma_glu_cys,
            idx_gly: indices.glycine,
            idx_atp: indices.atp,
            idx_adp: indices.adp,
            idx_gsh: indices.gsh,
        }
    }
}

impl Enzyme for GlutathioneSynthetase {
    fn name(&self) -> &'static str { "Glutathione Synthetase" }
    fn ec_number(&self) -> &'static str { "6.3.2.3" }

    fn rate(&self, metabolites: &MetabolitePool) -> f64 {
        let gamma_glu_cys = metabolites.get(self.idx_gamma_glu_cys);
        let gly = metabolites.get(self.idx_gly);
        let atp = metabolites.get(self.idx_atp);

        if gamma_glu_cys <= 0.0 || gly <= 0.0 || atp <= 0.0 {
            return 0.0;
        }

        let gamma_glu_cys_term = gamma_glu_cys / (self.km_gamma_glu_cys_mM + gamma_glu_cys);
        let gly_term = gly / (self.km_gly_mM + gly);
        let atp_term = atp / (self.km_atp_mM + atp);

        self.vmax_mM_per_sec * gamma_glu_cys_term * gly_term * atp_term
    }

    fn stoichiometry(&self) -> &ReactionStoichiometry {
        &self.stoichiometry
    }
}

// ============================================================================
// Glutathione Cycle Solver
// ============================================================================

/// Glutathione cycle solver combining detoxification and regeneration
pub struct GlutathioneCycle {
    pub indices: RedoxIndices,
    pub gpx: GlutathionePeroxidase,
    pub gr: GlutathioneReductase,
    pub gamma_gcs: GammaGlutamylcysteineSynthetase,
    pub gs: GlutathioneSynthetase,
    /// H2O2 production rate from autoxidation
    pub h2o2_production_rate_mM_per_sec: f64,
}

impl GlutathioneCycle {
    /// Create a new glutathione cycle solver
    pub fn new(indices: &RedoxIndices) -> Self {
        Self {
            gpx: GlutathionePeroxidase::new(indices),
            gr: GlutathioneReductase::new(indices),
            gamma_gcs: GammaGlutamylcysteineSynthetase::new(indices),
            gs: GlutathioneSynthetase::new(indices),
            indices: *indices,
            h2o2_production_rate_mM_per_sec: BASAL_H2O2_PRODUCTION_MM_PER_SEC,
        }
    }

    /// Compute all glutathione cycle derivatives
    pub fn compute_derivatives(&self, metabolites: &MetabolitePool, dydt: &mut [f64]) {
        // H2O2 detoxification
        self.gpx.apply_to_derivatives(metabolites, dydt);

        // GSH regeneration
        self.gr.apply_to_derivatives(metabolites, dydt);

        // GSH de novo synthesis
        self.gamma_gcs.apply_to_derivatives(metabolites, dydt);
        self.gs.apply_to_derivatives(metabolites, dydt);

        // Basal H2O2 production from Hb autoxidation
        dydt[self.indices.h2o2] += self.h2o2_production_rate_mM_per_sec;
    }

    /// Set oxidative stress level (multiplier for H2O2 production)
    pub fn set_oxidative_stress(&mut self, multiplier: f64) {
        self.h2o2_production_rate_mM_per_sec = BASAL_H2O2_PRODUCTION_MM_PER_SEC * multiplier;
    }

    /// Get all reaction rates for diagnostics
    pub fn get_rates(&self, metabolites: &MetabolitePool) -> Vec<(&'static str, f64)> {
        vec![
            (self.gpx.name(), self.gpx.rate(metabolites)),
            (self.gr.name(), self.gr.rate(metabolites)),
            (self.gamma_gcs.name(), self.gamma_gcs.rate(metabolites)),
            (self.gs.name(), self.gs.rate(metabolites)),
            ("H2O2 Production", self.h2o2_production_rate_mM_per_sec),
        ]
    }

    /// Calculate GSH/GSSG ratio
    pub fn gsh_gssg_ratio(&self, metabolites: &MetabolitePool) -> f64 {
        let gsh = metabolites.get(self.indices.gsh);
        let gssg = metabolites.get(self.indices.gssg);

        if gssg > 1e-9 {
            gsh / gssg
        } else {
            f64::INFINITY
        }
    }

    /// Get total glutathione (GSH + 2*GSSG)
    pub fn total_glutathione_mM(&self, metabolites: &MetabolitePool) -> f64 {
        let gsh = metabolites.get(self.indices.gsh);
        let gssg = metabolites.get(self.indices.gssg);
        gsh + 2.0 * gssg
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::glycolysis::MetaboliteIndices;

    fn create_test_indices() -> RedoxIndices {
        let glycolysis = MetaboliteIndices::default();
        RedoxIndices::new(&glycolysis, 18)
    }

    fn create_test_metabolites(indices: &RedoxIndices) -> MetabolitePool {
        let total = 18 + RedoxIndices::new_metabolite_count();
        let mut pool = MetabolitePool::new(total);

        // Glutathione - physiological concentrations
        pool.set(indices.gsh, 2.5);   // 2.5 mM GSH (target 2-3 mM)
        pool.set(indices.gssg, 0.01); // 0.01 mM GSSG (GSH/GSSG ~250)
        pool.set(indices.h2o2, 0.001); // 1 uM H2O2

        // NADPH/NADP+ for GR
        pool.set(indices.nadph, 0.3);
        pool.set(indices.nadp_plus, 0.02);

        // GSH synthesis precursors
        pool.set(indices.glutamate, 0.5);
        pool.set(indices.cysteine, 0.05);
        pool.set(indices.glycine, 0.5);
        pool.set(indices.gamma_glu_cys, 0.001);

        // ATP for synthesis
        pool.set(indices.atp, 2.0);
        pool.set(indices.adp, 0.25);

        pool
    }

    #[test]
    fn test_gpx_rate() {
        let indices = create_test_indices();
        let gpx = GlutathionePeroxidase::new(&indices);
        let metabolites = create_test_metabolites(&indices);

        let rate = gpx.rate(&metabolites);
        assert!(rate > 0.0, "GPx should have positive rate");
    }

    #[test]
    fn test_gpx_h2o2_affinity() {
        let indices = create_test_indices();
        let gpx = GlutathionePeroxidase::new(&indices);
        let mut metabolites = create_test_metabolites(&indices);

        // GPx should still work efficiently at very low H2O2 (high affinity)
        metabolites.set(indices.h2o2, 0.0001);  // 0.1 uM
        let rate_low = gpx.rate(&metabolites);

        metabolites.set(indices.h2o2, 0.01);  // 10 uM
        let rate_high = gpx.rate(&metabolites);

        assert!(rate_low > 0.0, "GPx should work at low H2O2");
        assert!(rate_high > rate_low, "Higher H2O2 should give higher rate");
    }

    #[test]
    fn test_gr_rate() {
        let indices = create_test_indices();
        let gr = GlutathioneReductase::new(&indices);
        let metabolites = create_test_metabolites(&indices);

        let rate = gr.rate(&metabolites);
        assert!(rate > 0.0, "GR should have positive rate");
    }

    #[test]
    fn test_gamma_gcs_feedback_inhibition() {
        let indices = create_test_indices();
        let gamma_gcs = GammaGlutamylcysteineSynthetase::new(&indices);
        let mut metabolites = create_test_metabolites(&indices);

        // Rate at low GSH
        metabolites.set(indices.gsh, 0.5);
        let rate_low_gsh = gamma_gcs.rate(&metabolites);

        // Rate at high GSH (should be inhibited)
        metabolites.set(indices.gsh, 5.0);
        let rate_high_gsh = gamma_gcs.rate(&metabolites);

        assert!(
            rate_high_gsh < rate_low_gsh,
            "gamma-GCS should be inhibited by high GSH"
        );
    }

    #[test]
    fn test_glutathione_cycle() {
        let indices = create_test_indices();
        let cycle = GlutathioneCycle::new(&indices);
        let metabolites = create_test_metabolites(&indices);

        let mut dydt = vec![0.0; 18 + RedoxIndices::new_metabolite_count()];
        cycle.compute_derivatives(&metabolites, &mut dydt);

        // H2O2 should be consumed by GPx but produced by autoxidation
        // Net should depend on balance
        assert!(dydt[indices.h2o2].is_finite(), "H2O2 derivative should be finite");

        // GSH should be both consumed (GPx) and regenerated (GR, synthesis)
        assert!(dydt[indices.gsh].is_finite(), "GSH derivative should be finite");
    }

    #[test]
    fn test_gsh_gssg_ratio() {
        let indices = create_test_indices();
        let cycle = GlutathioneCycle::new(&indices);
        let metabolites = create_test_metabolites(&indices);

        let ratio = cycle.gsh_gssg_ratio(&metabolites);
        assert!(ratio > 100.0, "GSH/GSSG ratio should be > 100");
        assert!(ratio < 500.0, "GSH/GSSG ratio should be < 500");
    }

    #[test]
    fn test_total_glutathione() {
        let indices = create_test_indices();
        let cycle = GlutathioneCycle::new(&indices);
        let metabolites = create_test_metabolites(&indices);

        let total = cycle.total_glutathione_mM(&metabolites);
        assert!(total > 2.0, "Total glutathione should be > 2 mM");
        assert!(total < 4.0, "Total glutathione should be < 4 mM");
    }
}
