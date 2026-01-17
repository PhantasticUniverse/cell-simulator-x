//! Pentose Phosphate Pathway (PPP) for RBC metabolism.
//!
//! The PPP provides NADPH for antioxidant defense and ribose-5-phosphate
//! for nucleotide synthesis. In mature RBCs, the primary function is NADPH
//! production for glutathione reduction.
//!
//! Pathway structure:
//! - Oxidative branch (irreversible): G6P -> NADPH + CO2 + Ru5P
//! - Non-oxidative branch (reversible): Ru5P <-> F6P + G3P (recycling)
//!
//! References:
//! - Beutler E. Red Cell Metabolism: A Manual of Biochemical Methods. 1984
//! - Veech RL et al. Biochem J. 1969;115:609-619 (NADPH/NADP+ ratios)
//! - Wood T. The Pentose Phosphate Pathway. Academic Press, 1985

use super::enzyme::{
    Enzyme, ReactionStoichiometry,
    reversible_mm,
};
use super::MetabolitePool;
use super::glycolysis::MetaboliteIndices;

/// Indices for redox metabolites in the concentration vector
#[derive(Debug, Clone, Copy)]
pub struct RedoxIndices {
    // Shared with glycolysis
    pub glucose_6_phosphate: usize,
    pub fructose_6_phosphate: usize,
    pub glyceraldehyde_3_phosphate: usize,

    // PPP intermediates
    pub phosphogluconolactone_6: usize,   // 6-PGL
    pub phosphogluconate_6: usize,         // 6-PG
    pub ribulose_5_phosphate: usize,       // Ru5P
    pub ribose_5_phosphate: usize,         // R5P
    pub xylulose_5_phosphate: usize,       // Xu5P
    pub sedoheptulose_7_phosphate: usize,  // S7P
    pub erythrose_4_phosphate: usize,      // E4P

    // Redox pairs
    pub nadph: usize,
    pub nadp_plus: usize,

    // Glutathione
    pub gsh: usize,      // Reduced glutathione
    pub gssg: usize,     // Oxidized glutathione
    pub h2o2: usize,     // Hydrogen peroxide

    // Piezo1/Calcium
    pub ca2_plus_cytosolic: usize,  // Cytosolic Ca2+ (in uM for numerical stability)

    // Glutathione synthesis precursors
    pub glutamate: usize,
    pub cysteine: usize,
    pub glycine: usize,
    pub gamma_glu_cys: usize,  // gamma-glutamylcysteine intermediate

    // ATP/ADP shared with glycolysis
    pub atp: usize,
    pub adp: usize,
}

impl RedoxIndices {
    /// Create redox indices starting after glycolysis indices
    pub fn new(glycolysis: &MetaboliteIndices, start_idx: usize) -> Self {
        Self {
            // Shared indices from glycolysis
            glucose_6_phosphate: glycolysis.glucose_6_phosphate,
            fructose_6_phosphate: glycolysis.fructose_6_phosphate,
            glyceraldehyde_3_phosphate: glycolysis.glyceraldehyde_3_phosphate,
            atp: glycolysis.atp,
            adp: glycolysis.adp,

            // New PPP intermediates (starting at start_idx)
            phosphogluconolactone_6: start_idx,
            phosphogluconate_6: start_idx + 1,
            ribulose_5_phosphate: start_idx + 2,
            ribose_5_phosphate: start_idx + 3,
            xylulose_5_phosphate: start_idx + 4,
            sedoheptulose_7_phosphate: start_idx + 5,
            erythrose_4_phosphate: start_idx + 6,

            // Redox pairs
            nadph: start_idx + 7,
            nadp_plus: start_idx + 8,

            // Glutathione
            gsh: start_idx + 9,
            gssg: start_idx + 10,
            h2o2: start_idx + 11,

            // Calcium
            ca2_plus_cytosolic: start_idx + 12,

            // GSH synthesis precursors
            glutamate: start_idx + 13,
            cysteine: start_idx + 14,
            glycine: start_idx + 15,
            gamma_glu_cys: start_idx + 16,
        }
    }

    /// Number of new metabolites added by redox module
    pub fn new_metabolite_count() -> usize {
        17  // All indices from phosphogluconolactone_6 through gamma_glu_cys
    }
}

// ============================================================================
// Oxidative Branch Enzymes
// ============================================================================

/// Glucose-6-Phosphate Dehydrogenase (G6PDH) - EC 1.1.1.49
///
/// G6P + NADP+ -> 6-PGL + NADPH + H+
///
/// Rate-limiting step of PPP. Deficiency causes hemolytic anemia under
/// oxidative stress. Strong product inhibition by NADPH maintains low
/// flux under normal conditions.
///
/// Reference: Beutler 1984, Vmax = 0.84 mM/s, Km_G6P = 0.039 mM
pub struct Glucose6PhosphateDehydrogenase {
    pub vmax_mM_per_sec: f64,
    pub km_g6p_mM: f64,
    pub km_nadp_mM: f64,
    pub ki_nadph_mM: f64,  // Competitive inhibition constant
    stoichiometry: ReactionStoichiometry,
    idx_g6p: usize,
    idx_nadp: usize,
    idx_nadph: usize,
}

impl Glucose6PhosphateDehydrogenase {
    pub fn new(indices: &RedoxIndices) -> Self {
        Self {
            // Reference: Beutler 1984, Kirkman 1986
            // Vmax set to maintain NADPH/NADP+ ratio 10-20 at steady state
            vmax_mM_per_sec: 0.06,
            km_g6p_mM: 0.039,
            km_nadp_mM: 0.005,
            ki_nadph_mM: 0.005, // Strong NADPH inhibition
            stoichiometry: ReactionStoichiometry::new(
                vec![
                    (indices.glucose_6_phosphate, 1.0),
                    (indices.nadp_plus, 1.0),
                ],
                vec![
                    (indices.phosphogluconolactone_6, 1.0),
                    (indices.nadph, 1.0),
                ],
            ),
            idx_g6p: indices.glucose_6_phosphate,
            idx_nadp: indices.nadp_plus,
            idx_nadph: indices.nadph,
        }
    }
}

impl Enzyme for Glucose6PhosphateDehydrogenase {
    fn name(&self) -> &'static str { "G6P Dehydrogenase" }
    fn ec_number(&self) -> &'static str { "1.1.1.49" }

    fn rate(&self, metabolites: &MetabolitePool) -> f64 {
        let g6p = metabolites.get(self.idx_g6p);
        let nadp = metabolites.get(self.idx_nadp);
        let nadph = metabolites.get(self.idx_nadph);

        if g6p <= 0.0 || nadp <= 0.0 {
            return 0.0;
        }

        // Bi-substrate kinetics with competitive inhibition by NADPH
        // v = Vmax * [G6P] * [NADP+] / ((Km_G6P + [G6P]) * (Km_NADP * (1 + [NADPH]/Ki) + [NADP+]))
        let km_nadp_apparent = self.km_nadp_mM * (1.0 + nadph / self.ki_nadph_mM);
        let g6p_term = g6p / (self.km_g6p_mM + g6p);
        let nadp_term = nadp / (km_nadp_apparent + nadp);

        self.vmax_mM_per_sec * g6p_term * nadp_term
    }

    fn stoichiometry(&self) -> &ReactionStoichiometry {
        &self.stoichiometry
    }
}

/// 6-Phosphogluconolactonase (PGL) - EC 3.1.1.31
///
/// 6-PGL -> 6-PG + H2O
///
/// Very fast, essentially spontaneous reaction.
/// Reference: Rate constant k = 10/s (Wood 1985)
pub struct Phosphogluconolactonase {
    pub k_sec: f64,  // First-order rate constant
    stoichiometry: ReactionStoichiometry,
    idx_pgl: usize,
}

impl Phosphogluconolactonase {
    pub fn new(indices: &RedoxIndices) -> Self {
        Self {
            k_sec: 10.0,  // Fast hydrolysis
            stoichiometry: ReactionStoichiometry::new(
                vec![(indices.phosphogluconolactone_6, 1.0)],
                vec![(indices.phosphogluconate_6, 1.0)],
            ),
            idx_pgl: indices.phosphogluconolactone_6,
        }
    }
}

impl Enzyme for Phosphogluconolactonase {
    fn name(&self) -> &'static str { "Phosphogluconolactonase" }
    fn ec_number(&self) -> &'static str { "3.1.1.31" }

    fn rate(&self, metabolites: &MetabolitePool) -> f64 {
        let pgl = metabolites.get(self.idx_pgl);
        // First-order kinetics: v = k * [S]
        self.k_sec * pgl
    }

    fn stoichiometry(&self) -> &ReactionStoichiometry {
        &self.stoichiometry
    }
}

/// 6-Phosphogluconate Dehydrogenase (6PGDH) - EC 1.1.1.44
///
/// 6-PG + NADP+ -> Ru5P + CO2 + NADPH
///
/// Second NADPH-producing step in oxidative branch.
/// Reference: Beutler 1984, Vmax = 0.42 mM/s
pub struct Phosphogluconate6Dehydrogenase {
    pub vmax_mM_per_sec: f64,
    pub km_6pg_mM: f64,
    pub km_nadp_mM: f64,
    stoichiometry: ReactionStoichiometry,
    idx_6pg: usize,
    idx_nadp: usize,
}

impl Phosphogluconate6Dehydrogenase {
    pub fn new(indices: &RedoxIndices) -> Self {
        Self {
            // Vmax scaled with G6PDH to maintain pathway balance
            vmax_mM_per_sec: 0.04,
            km_6pg_mM: 0.035,
            km_nadp_mM: 0.008,
            stoichiometry: ReactionStoichiometry::new(
                vec![
                    (indices.phosphogluconate_6, 1.0),
                    (indices.nadp_plus, 1.0),
                ],
                vec![
                    (indices.ribulose_5_phosphate, 1.0),
                    (indices.nadph, 1.0),
                    // CO2 is released (not tracked)
                ],
            ),
            idx_6pg: indices.phosphogluconate_6,
            idx_nadp: indices.nadp_plus,
        }
    }
}

impl Enzyme for Phosphogluconate6Dehydrogenase {
    fn name(&self) -> &'static str { "6PG Dehydrogenase" }
    fn ec_number(&self) -> &'static str { "1.1.1.44" }

    fn rate(&self, metabolites: &MetabolitePool) -> f64 {
        let pg6 = metabolites.get(self.idx_6pg);
        let nadp = metabolites.get(self.idx_nadp);

        if pg6 <= 0.0 || nadp <= 0.0 {
            return 0.0;
        }

        // Bi-substrate Michaelis-Menten
        let pg6_term = pg6 / (self.km_6pg_mM + pg6);
        let nadp_term = nadp / (self.km_nadp_mM + nadp);

        self.vmax_mM_per_sec * pg6_term * nadp_term
    }

    fn stoichiometry(&self) -> &ReactionStoichiometry {
        &self.stoichiometry
    }
}

// ============================================================================
// Non-oxidative Branch Enzymes
// ============================================================================

/// Ribulose-5-Phosphate Epimerase (RPE) - EC 5.1.3.1
///
/// Ru5P <-> Xu5P
///
/// Reversible isomerization, Keq = 1.4 ([Xu5P]/[Ru5P])
/// Reference: Wood 1985
pub struct RibulosePhosphateEpimerase {
    pub vmax_f_mM_per_sec: f64,
    pub vmax_r_mM_per_sec: f64,
    pub km_ru5p_mM: f64,
    pub km_xu5p_mM: f64,
    stoichiometry: ReactionStoichiometry,
    idx_ru5p: usize,
    idx_xu5p: usize,
}

impl RibulosePhosphateEpimerase {
    pub fn new(indices: &RedoxIndices) -> Self {
        // Keq = 1.4, so Vmax_r = Vmax_f * Km_xu5p / (Keq * Km_ru5p)
        let keq = 1.4;
        let km_ru5p = 0.5;
        let km_xu5p = 0.5;
        let vmax_f = 2.0;  // Fast equilibrating
        let vmax_r = vmax_f * km_xu5p / (keq * km_ru5p);

        Self {
            vmax_f_mM_per_sec: vmax_f,
            vmax_r_mM_per_sec: vmax_r,
            km_ru5p_mM: km_ru5p,
            km_xu5p_mM: km_xu5p,
            stoichiometry: ReactionStoichiometry::new(
                vec![(indices.ribulose_5_phosphate, 1.0)],
                vec![(indices.xylulose_5_phosphate, 1.0)],
            ),
            idx_ru5p: indices.ribulose_5_phosphate,
            idx_xu5p: indices.xylulose_5_phosphate,
        }
    }
}

impl Enzyme for RibulosePhosphateEpimerase {
    fn name(&self) -> &'static str { "Ru5P Epimerase" }
    fn ec_number(&self) -> &'static str { "5.1.3.1" }

    fn rate(&self, metabolites: &MetabolitePool) -> f64 {
        let ru5p = metabolites.get(self.idx_ru5p);
        let xu5p = metabolites.get(self.idx_xu5p);

        reversible_mm(
            self.vmax_f_mM_per_sec,
            self.km_ru5p_mM,
            ru5p,
            self.vmax_r_mM_per_sec,
            self.km_xu5p_mM,
            xu5p,
        )
    }

    fn stoichiometry(&self) -> &ReactionStoichiometry {
        &self.stoichiometry
    }
}

/// Ribose-5-Phosphate Isomerase (RPI) - EC 5.3.1.6
///
/// Ru5P <-> R5P
///
/// Reversible isomerization, Keq = 3.0 ([R5P]/[Ru5P])
/// Reference: Wood 1985
pub struct RibosePhosphateIsomerase {
    pub vmax_f_mM_per_sec: f64,
    pub vmax_r_mM_per_sec: f64,
    pub km_ru5p_mM: f64,
    pub km_r5p_mM: f64,
    stoichiometry: ReactionStoichiometry,
    idx_ru5p: usize,
    idx_r5p: usize,
}

impl RibosePhosphateIsomerase {
    pub fn new(indices: &RedoxIndices) -> Self {
        let keq = 3.0;
        let km_ru5p = 0.8;
        let km_r5p = 1.5;
        let vmax_f = 2.0;
        let vmax_r = vmax_f * km_r5p / (keq * km_ru5p);

        Self {
            vmax_f_mM_per_sec: vmax_f,
            vmax_r_mM_per_sec: vmax_r,
            km_ru5p_mM: km_ru5p,
            km_r5p_mM: km_r5p,
            stoichiometry: ReactionStoichiometry::new(
                vec![(indices.ribulose_5_phosphate, 1.0)],
                vec![(indices.ribose_5_phosphate, 1.0)],
            ),
            idx_ru5p: indices.ribulose_5_phosphate,
            idx_r5p: indices.ribose_5_phosphate,
        }
    }
}

impl Enzyme for RibosePhosphateIsomerase {
    fn name(&self) -> &'static str { "R5P Isomerase" }
    fn ec_number(&self) -> &'static str { "5.3.1.6" }

    fn rate(&self, metabolites: &MetabolitePool) -> f64 {
        let ru5p = metabolites.get(self.idx_ru5p);
        let r5p = metabolites.get(self.idx_r5p);

        reversible_mm(
            self.vmax_f_mM_per_sec,
            self.km_ru5p_mM,
            ru5p,
            self.vmax_r_mM_per_sec,
            self.km_r5p_mM,
            r5p,
        )
    }

    fn stoichiometry(&self) -> &ReactionStoichiometry {
        &self.stoichiometry
    }
}

/// Transketolase (TK) - EC 2.2.1.1
///
/// Catalyzes two reactions in the non-oxidative branch:
/// 1. Xu5P + R5P <-> G3P + S7P
/// 2. Xu5P + E4P <-> G3P + F6P
///
/// We model the net flux as primarily reaction 1 for simplicity.
/// Reference: Wood 1985, Beutler 1984
pub struct Transketolase {
    pub vmax_f_mM_per_sec: f64,
    pub km_xu5p_mM: f64,
    pub km_r5p_mM: f64,
    pub keq: f64,
    stoichiometry: ReactionStoichiometry,
    idx_xu5p: usize,
    idx_r5p: usize,
    idx_g3p: usize,
    idx_s7p: usize,
}

impl Transketolase {
    pub fn new(indices: &RedoxIndices) -> Self {
        Self {
            vmax_f_mM_per_sec: 1.0,
            km_xu5p_mM: 0.15,
            km_r5p_mM: 0.5,
            keq: 1.2,  // Slightly favors products
            stoichiometry: ReactionStoichiometry::new(
                vec![
                    (indices.xylulose_5_phosphate, 1.0),
                    (indices.ribose_5_phosphate, 1.0),
                ],
                vec![
                    (indices.glyceraldehyde_3_phosphate, 1.0),
                    (indices.sedoheptulose_7_phosphate, 1.0),
                ],
            ),
            idx_xu5p: indices.xylulose_5_phosphate,
            idx_r5p: indices.ribose_5_phosphate,
            idx_g3p: indices.glyceraldehyde_3_phosphate,
            idx_s7p: indices.sedoheptulose_7_phosphate,
        }
    }
}

impl Enzyme for Transketolase {
    fn name(&self) -> &'static str { "Transketolase" }
    fn ec_number(&self) -> &'static str { "2.2.1.1" }

    fn rate(&self, metabolites: &MetabolitePool) -> f64 {
        let xu5p = metabolites.get(self.idx_xu5p);
        let r5p = metabolites.get(self.idx_r5p);
        let g3p = metabolites.get(self.idx_g3p);
        let s7p = metabolites.get(self.idx_s7p);

        if xu5p <= 0.0 || r5p <= 0.0 {
            // Check for reverse reaction
            if g3p <= 0.0 || s7p <= 0.0 {
                return 0.0;
            }
        }

        // Bi-bi kinetics with mass action ratio correction
        let forward = self.vmax_f_mM_per_sec * xu5p * r5p /
            ((self.km_xu5p_mM + xu5p) * (self.km_r5p_mM + r5p));

        // Mass action ratio for equilibrium
        let substrate_product = xu5p * r5p;
        let product_product = g3p * s7p;
        let mass_action = if substrate_product > 1e-12 {
            product_product / (substrate_product * self.keq)
        } else {
            0.0
        };

        forward * (1.0 - mass_action).max(-1.0).min(1.0)
    }

    fn stoichiometry(&self) -> &ReactionStoichiometry {
        &self.stoichiometry
    }
}

/// Transaldolase (TA) - EC 2.2.1.2
///
/// S7P + G3P <-> E4P + F6P
///
/// Links transketolase to glycolysis by regenerating F6P.
/// Reference: Wood 1985
pub struct Transaldolase {
    pub vmax_f_mM_per_sec: f64,
    pub km_s7p_mM: f64,
    pub km_g3p_mM: f64,
    pub keq: f64,
    stoichiometry: ReactionStoichiometry,
    idx_s7p: usize,
    idx_g3p: usize,
    idx_e4p: usize,
    idx_f6p: usize,
}

impl Transaldolase {
    pub fn new(indices: &RedoxIndices) -> Self {
        Self {
            vmax_f_mM_per_sec: 0.8,
            km_s7p_mM: 0.3,
            km_g3p_mM: 0.1,
            keq: 1.05,
            stoichiometry: ReactionStoichiometry::new(
                vec![
                    (indices.sedoheptulose_7_phosphate, 1.0),
                    (indices.glyceraldehyde_3_phosphate, 1.0),
                ],
                vec![
                    (indices.erythrose_4_phosphate, 1.0),
                    (indices.fructose_6_phosphate, 1.0),
                ],
            ),
            idx_s7p: indices.sedoheptulose_7_phosphate,
            idx_g3p: indices.glyceraldehyde_3_phosphate,
            idx_e4p: indices.erythrose_4_phosphate,
            idx_f6p: indices.fructose_6_phosphate,
        }
    }
}

impl Enzyme for Transaldolase {
    fn name(&self) -> &'static str { "Transaldolase" }
    fn ec_number(&self) -> &'static str { "2.2.1.2" }

    fn rate(&self, metabolites: &MetabolitePool) -> f64 {
        let s7p = metabolites.get(self.idx_s7p);
        let g3p = metabolites.get(self.idx_g3p);
        let e4p = metabolites.get(self.idx_e4p);
        let f6p = metabolites.get(self.idx_f6p);

        if s7p <= 0.0 || g3p <= 0.0 {
            if e4p <= 0.0 || f6p <= 0.0 {
                return 0.0;
            }
        }

        let forward = self.vmax_f_mM_per_sec * s7p * g3p /
            ((self.km_s7p_mM + s7p) * (self.km_g3p_mM + g3p));

        let substrate_product = s7p * g3p;
        let product_product = e4p * f6p;
        let mass_action = if substrate_product > 1e-12 {
            product_product / (substrate_product * self.keq)
        } else {
            0.0
        };

        forward * (1.0 - mass_action).max(-1.0).min(1.0)
    }

    fn stoichiometry(&self) -> &ReactionStoichiometry {
        &self.stoichiometry
    }
}

// ============================================================================
// Pentose Phosphate Pathway Solver
// ============================================================================

/// Pentose Phosphate Pathway solver combining all 7 enzymes
pub struct PentosePhosphatePathway {
    pub indices: RedoxIndices,
    pub g6pdh: Glucose6PhosphateDehydrogenase,
    pub pgl: Phosphogluconolactonase,
    pub pgdh_6: Phosphogluconate6Dehydrogenase,
    pub rpe: RibulosePhosphateEpimerase,
    pub rpi: RibosePhosphateIsomerase,
    pub tk: Transketolase,
    pub ta: Transaldolase,
}

impl PentosePhosphatePathway {
    /// Create a new PPP solver
    pub fn new(indices: &RedoxIndices) -> Self {
        Self {
            g6pdh: Glucose6PhosphateDehydrogenase::new(indices),
            pgl: Phosphogluconolactonase::new(indices),
            pgdh_6: Phosphogluconate6Dehydrogenase::new(indices),
            rpe: RibulosePhosphateEpimerase::new(indices),
            rpi: RibosePhosphateIsomerase::new(indices),
            tk: Transketolase::new(indices),
            ta: Transaldolase::new(indices),
            indices: *indices,
        }
    }

    /// Compute all PPP reaction rates and apply to derivatives
    pub fn compute_derivatives(&self, metabolites: &MetabolitePool, dydt: &mut [f64]) {
        // Oxidative branch
        self.g6pdh.apply_to_derivatives(metabolites, dydt);
        self.pgl.apply_to_derivatives(metabolites, dydt);
        self.pgdh_6.apply_to_derivatives(metabolites, dydt);

        // Non-oxidative branch
        self.rpe.apply_to_derivatives(metabolites, dydt);
        self.rpi.apply_to_derivatives(metabolites, dydt);
        self.tk.apply_to_derivatives(metabolites, dydt);
        self.ta.apply_to_derivatives(metabolites, dydt);
    }

    /// Get all reaction rates for diagnostics
    pub fn get_rates(&self, metabolites: &MetabolitePool) -> Vec<(&'static str, f64)> {
        vec![
            (self.g6pdh.name(), self.g6pdh.rate(metabolites)),
            (self.pgl.name(), self.pgl.rate(metabolites)),
            (self.pgdh_6.name(), self.pgdh_6.rate(metabolites)),
            (self.rpe.name(), self.rpe.rate(metabolites)),
            (self.rpi.name(), self.rpi.rate(metabolites)),
            (self.tk.name(), self.tk.rate(metabolites)),
            (self.ta.name(), self.ta.rate(metabolites)),
        ]
    }

    /// Get total NADPH production rate (from oxidative branch)
    pub fn nadph_production_rate(&self, metabolites: &MetabolitePool) -> f64 {
        self.g6pdh.rate(metabolites) + self.pgdh_6.rate(metabolites)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_test_indices() -> RedoxIndices {
        let glycolysis = MetaboliteIndices::default();
        RedoxIndices::new(&glycolysis, 18)  // Start after 2,3-BPG
    }

    fn create_test_metabolites(indices: &RedoxIndices) -> MetabolitePool {
        let total = 18 + RedoxIndices::new_metabolite_count();
        let mut pool = MetabolitePool::new(total);

        // Set glycolysis-shared metabolites
        pool.set(indices.glucose_6_phosphate, 0.05);
        pool.set(indices.fructose_6_phosphate, 0.02);
        pool.set(indices.glyceraldehyde_3_phosphate, 0.005);

        // Set PPP intermediates
        pool.set(indices.phosphogluconolactone_6, 0.001);
        pool.set(indices.phosphogluconate_6, 0.01);
        pool.set(indices.ribulose_5_phosphate, 0.01);
        pool.set(indices.ribose_5_phosphate, 0.01);
        pool.set(indices.xylulose_5_phosphate, 0.01);
        pool.set(indices.sedoheptulose_7_phosphate, 0.01);
        pool.set(indices.erythrose_4_phosphate, 0.01);

        // Set redox pairs - physiological ratio NADPH/NADP+ = 10-20
        pool.set(indices.nadph, 0.3);  // ~0.3 mM NADPH
        pool.set(indices.nadp_plus, 0.02);  // ~0.02 mM NADP+

        pool
    }

    #[test]
    fn test_g6pdh_rate() {
        let indices = create_test_indices();
        let g6pdh = Glucose6PhosphateDehydrogenase::new(&indices);
        let metabolites = create_test_metabolites(&indices);

        let rate = g6pdh.rate(&metabolites);
        assert!(rate > 0.0, "G6PDH should have positive rate");
        assert!(rate < g6pdh.vmax_mM_per_sec, "Rate should be below Vmax");
    }

    #[test]
    fn test_g6pdh_nadph_inhibition() {
        let indices = create_test_indices();
        let g6pdh = Glucose6PhosphateDehydrogenase::new(&indices);
        let mut metabolites = create_test_metabolites(&indices);

        // Rate at low NADPH
        metabolites.set(indices.nadph, 0.01);
        let rate_low = g6pdh.rate(&metabolites);

        // Rate at high NADPH (should be inhibited)
        metabolites.set(indices.nadph, 0.5);
        let rate_high = g6pdh.rate(&metabolites);

        assert!(rate_high < rate_low, "G6PDH should be inhibited by NADPH");
    }

    #[test]
    fn test_ppp_solver() {
        let indices = create_test_indices();
        let ppp = PentosePhosphatePathway::new(&indices);
        let metabolites = create_test_metabolites(&indices);

        let mut dydt = vec![0.0; 18 + RedoxIndices::new_metabolite_count()];
        ppp.compute_derivatives(&metabolites, &mut dydt);

        // NADPH should be produced
        assert!(dydt[indices.nadph] > 0.0, "NADPH should be produced");

        // NADP+ should be consumed
        assert!(dydt[indices.nadp_plus] < 0.0, "NADP+ should be consumed");
    }

    #[test]
    fn test_nadph_production_rate() {
        let indices = create_test_indices();
        let ppp = PentosePhosphatePathway::new(&indices);
        let metabolites = create_test_metabolites(&indices);

        let rate = ppp.nadph_production_rate(&metabolites);
        assert!(rate > 0.0, "NADPH production should be positive");
    }
}
