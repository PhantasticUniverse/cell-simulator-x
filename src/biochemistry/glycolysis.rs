//! Glycolysis pathway enzymes for RBC metabolism.
//!
//! Implements the 10 reactions of glycolysis plus lactate dehydrogenase.
//! RBCs are unique in that they lack mitochondria and rely entirely on
//! anaerobic glycolysis for ATP production.
//!
//! Pathway: Glucose → 2 Lactate + 2 ATP (net)
//!
//! References:
//! - Rapoport TA et al. Eur J Biochem. 1976;69:571-584 (RBC glycolysis model)
//! - Joshi A, Palsson BO. J Theor Biol. 1989;141:515-528 (Detailed kinetics)
//! - Mulquiney PJ, Kuchel PW. Biochem J. 1999;342:567-580 (RBC model)

use super::enzyme::{
    Enzyme, ReactionStoichiometry,
    michaelis_menten, hill_kinetics, reversible_mm, ordered_bi_bi,
    reversible_ordered_bi_bi,
};
use super::MetabolitePool;

/// Indices for metabolites in the concentration vector
#[derive(Debug, Clone, Copy)]
pub struct MetaboliteIndices {
    // External/transport
    pub glucose: usize,
    pub lactate: usize,

    // Glycolytic intermediates
    pub glucose_6_phosphate: usize,
    pub fructose_6_phosphate: usize,
    pub fructose_1_6_bisphosphate: usize,
    pub dihydroxyacetone_phosphate: usize,
    pub glyceraldehyde_3_phosphate: usize,
    pub bisphosphoglycerate_1_3: usize,
    pub phosphoglycerate_3: usize,
    pub phosphoglycerate_2: usize,
    pub phosphoenolpyruvate: usize,
    pub pyruvate: usize,

    // Cofactors
    pub atp: usize,
    pub adp: usize,
    pub nad: usize,
    pub nadh: usize,
    pub pi: usize,  // Inorganic phosphate
}

impl Default for MetaboliteIndices {
    fn default() -> Self {
        Self {
            glucose: 0,
            glucose_6_phosphate: 1,
            fructose_6_phosphate: 2,
            fructose_1_6_bisphosphate: 3,
            dihydroxyacetone_phosphate: 4,
            glyceraldehyde_3_phosphate: 5,
            bisphosphoglycerate_1_3: 6,
            phosphoglycerate_3: 7,
            phosphoglycerate_2: 8,
            phosphoenolpyruvate: 9,
            pyruvate: 10,
            lactate: 11,
            atp: 12,
            adp: 13,
            nad: 14,
            nadh: 15,
            pi: 16,
        }
    }
}

/// Glycolysis solver containing all 11 enzymes
pub struct GlycolysisSolver {
    pub indices: MetaboliteIndices,
    pub hexokinase: Hexokinase,
    pub glucose_6_p_isomerase: Glucose6PIsomerase,
    pub phosphofructokinase: Phosphofructokinase,
    pub aldolase: Aldolase,
    pub triose_p_isomerase: TriosePIsomerase,
    pub gapdh: GAPDH,
    pub phosphoglycerate_kinase: PhosphoglycerateKinase,
    pub phosphoglycerate_mutase: PhosphoglycerateMutase,
    pub enolase: Enolase,
    pub pyruvate_kinase: PyruvateKinase,
    pub lactate_dehydrogenase: LactateDehydrogenase,
}

impl GlycolysisSolver {
    /// Create a new glycolysis solver with default parameters
    pub fn new() -> Self {
        let indices = MetaboliteIndices::default();
        Self {
            hexokinase: Hexokinase::new(&indices),
            glucose_6_p_isomerase: Glucose6PIsomerase::new(&indices),
            phosphofructokinase: Phosphofructokinase::new(&indices),
            aldolase: Aldolase::new(&indices),
            triose_p_isomerase: TriosePIsomerase::new(&indices),
            gapdh: GAPDH::new(&indices),
            phosphoglycerate_kinase: PhosphoglycerateKinase::new(&indices),
            phosphoglycerate_mutase: PhosphoglycerateMutase::new(&indices),
            enolase: Enolase::new(&indices),
            pyruvate_kinase: PyruvateKinase::new(&indices),
            lactate_dehydrogenase: LactateDehydrogenase::new(&indices),
            indices,
        }
    }

    /// Compute all reaction rates and apply to derivatives
    pub fn compute_derivatives(&self, metabolites: &MetabolitePool, dydt: &mut [f64]) {
        self.hexokinase.apply_to_derivatives(metabolites, dydt);
        self.glucose_6_p_isomerase.apply_to_derivatives(metabolites, dydt);
        self.phosphofructokinase.apply_to_derivatives(metabolites, dydt);
        self.aldolase.apply_to_derivatives(metabolites, dydt);
        self.triose_p_isomerase.apply_to_derivatives(metabolites, dydt);
        self.gapdh.apply_to_derivatives(metabolites, dydt);
        self.phosphoglycerate_kinase.apply_to_derivatives(metabolites, dydt);
        self.phosphoglycerate_mutase.apply_to_derivatives(metabolites, dydt);
        self.enolase.apply_to_derivatives(metabolites, dydt);
        self.pyruvate_kinase.apply_to_derivatives(metabolites, dydt);
        self.lactate_dehydrogenase.apply_to_derivatives(metabolites, dydt);
    }

    /// Get all reaction rates for diagnostics
    pub fn get_rates(&self, metabolites: &MetabolitePool) -> Vec<(&'static str, f64)> {
        vec![
            (self.hexokinase.name(), self.hexokinase.rate(metabolites)),
            (self.glucose_6_p_isomerase.name(), self.glucose_6_p_isomerase.rate(metabolites)),
            (self.phosphofructokinase.name(), self.phosphofructokinase.rate(metabolites)),
            (self.aldolase.name(), self.aldolase.rate(metabolites)),
            (self.triose_p_isomerase.name(), self.triose_p_isomerase.rate(metabolites)),
            (self.gapdh.name(), self.gapdh.rate(metabolites)),
            (self.phosphoglycerate_kinase.name(), self.phosphoglycerate_kinase.rate(metabolites)),
            (self.phosphoglycerate_mutase.name(), self.phosphoglycerate_mutase.rate(metabolites)),
            (self.enolase.name(), self.enolase.rate(metabolites)),
            (self.pyruvate_kinase.name(), self.pyruvate_kinase.rate(metabolites)),
            (self.lactate_dehydrogenase.name(), self.lactate_dehydrogenase.rate(metabolites)),
        ]
    }
}

impl Default for GlycolysisSolver {
    fn default() -> Self {
        Self::new()
    }
}

// ============================================================================
// Individual Enzyme Implementations
// ============================================================================

/// Hexokinase (HK) - EC 2.7.1.1
///
/// Glucose + ATP → Glucose-6-P + ADP
///
/// Rate-limiting step, inhibited by G6P (product inhibition)
/// Reference: Rapoport 1976, Vmax = 0.033 mM/s, Km_glc = 0.1 mM
pub struct Hexokinase {
    pub vmax_mM_per_sec: f64,
    pub km_glucose_mM: f64,
    pub km_atp_mM: f64,
    pub ki_g6p_mM: f64,
    stoichiometry: ReactionStoichiometry,
    idx_glucose: usize,
    idx_atp: usize,
    idx_g6p: usize,
}

impl Hexokinase {
    pub fn new(indices: &MetaboliteIndices) -> Self {
        Self {
            // Reference: Rapoport 1976, 0.033 mM/s
            // Slightly elevated to support ATP homeostasis with pump load
            vmax_mM_per_sec: 0.035,
            km_glucose_mM: 0.1,
            km_atp_mM: 1.0,
            // Original G6P inhibition - maintains feedback control
            ki_g6p_mM: 0.5,
            stoichiometry: ReactionStoichiometry::new(
                vec![(indices.glucose, 1.0), (indices.atp, 1.0)],
                vec![(indices.glucose_6_phosphate, 1.0), (indices.adp, 1.0)],
            ),
            idx_glucose: indices.glucose,
            idx_atp: indices.atp,
            idx_g6p: indices.glucose_6_phosphate,
        }
    }
}

impl Enzyme for Hexokinase {
    fn name(&self) -> &'static str { "Hexokinase" }
    fn ec_number(&self) -> &'static str { "2.7.1.1" }

    fn rate(&self, metabolites: &MetabolitePool) -> f64 {
        let glc = metabolites.get(self.idx_glucose);
        let atp = metabolites.get(self.idx_atp);
        let g6p = metabolites.get(self.idx_g6p);

        // Ordered bi-bi with product inhibition by G6P
        let base_rate = ordered_bi_bi(
            self.vmax_mM_per_sec,
            self.km_glucose_mM,
            self.km_atp_mM,
            self.km_glucose_mM, // Ki_glc ≈ Km_glc
            glc,
            atp,
        );

        // Product inhibition by G6P
        base_rate / (1.0 + g6p / self.ki_g6p_mM)
    }

    fn stoichiometry(&self) -> &ReactionStoichiometry {
        &self.stoichiometry
    }
}

/// Glucose-6-phosphate Isomerase (GPI) - EC 5.3.1.9
///
/// Glucose-6-P ⇌ Fructose-6-P
///
/// Near-equilibrium reaction, Keq ≈ 0.4
/// Reference: Rapoport 1976
pub struct Glucose6PIsomerase {
    pub vmax_f_mM_per_sec: f64,
    pub vmax_r_mM_per_sec: f64,
    pub km_g6p_mM: f64,
    pub km_f6p_mM: f64,
    stoichiometry: ReactionStoichiometry,
    idx_g6p: usize,
    idx_f6p: usize,
}

impl Glucose6PIsomerase {
    pub fn new(indices: &MetaboliteIndices) -> Self {
        Self {
            vmax_f_mM_per_sec: 1.0,  // Fast equilibrating
            vmax_r_mM_per_sec: 2.5,  // To satisfy Keq = 0.4
            km_g6p_mM: 0.5,
            km_f6p_mM: 0.1,
            stoichiometry: ReactionStoichiometry::new(
                vec![(indices.glucose_6_phosphate, 1.0)],
                vec![(indices.fructose_6_phosphate, 1.0)],
            ),
            idx_g6p: indices.glucose_6_phosphate,
            idx_f6p: indices.fructose_6_phosphate,
        }
    }
}

impl Enzyme for Glucose6PIsomerase {
    fn name(&self) -> &'static str { "Glucose-6-P Isomerase" }
    fn ec_number(&self) -> &'static str { "5.3.1.9" }

    fn rate(&self, metabolites: &MetabolitePool) -> f64 {
        let g6p = metabolites.get(self.idx_g6p);
        let f6p = metabolites.get(self.idx_f6p);

        reversible_mm(
            self.vmax_f_mM_per_sec,
            self.km_g6p_mM,
            g6p,
            self.vmax_r_mM_per_sec,
            self.km_f6p_mM,
            f6p,
        )
    }

    fn stoichiometry(&self) -> &ReactionStoichiometry {
        &self.stoichiometry
    }
}

/// Phosphofructokinase (PFK) - EC 2.7.1.11
///
/// Fructose-6-P + ATP → Fructose-1,6-bisP + ADP
///
/// Key regulatory enzyme with cooperative kinetics
/// Activated by ADP/AMP, inhibited by ATP, citrate
/// Reference: Rapoport 1976, Joshi-Palsson 1989
pub struct Phosphofructokinase {
    pub vmax_mM_per_sec: f64,
    pub k_half_f6p_mM: f64,
    pub km_atp_mM: f64,
    pub hill_coefficient: f64,
    stoichiometry: ReactionStoichiometry,
    idx_f6p: usize,
    idx_atp: usize,
}

impl Phosphofructokinase {
    pub fn new(indices: &MetaboliteIndices) -> Self {
        Self {
            // Increased Vmax to support flux through lower glycolysis
            vmax_mM_per_sec: 0.4,
            k_half_f6p_mM: 0.15,
            km_atp_mM: 0.1,
            hill_coefficient: 2.5,
            stoichiometry: ReactionStoichiometry::new(
                vec![(indices.fructose_6_phosphate, 1.0), (indices.atp, 1.0)],
                vec![(indices.fructose_1_6_bisphosphate, 1.0), (indices.adp, 1.0)],
            ),
            idx_f6p: indices.fructose_6_phosphate,
            idx_atp: indices.atp,
        }
    }
}

impl Enzyme for Phosphofructokinase {
    fn name(&self) -> &'static str { "Phosphofructokinase" }
    fn ec_number(&self) -> &'static str { "2.7.1.11" }

    fn rate(&self, metabolites: &MetabolitePool) -> f64 {
        let f6p = metabolites.get(self.idx_f6p);
        let atp = metabolites.get(self.idx_atp);

        // Hill kinetics for F6P (cooperative), MM for ATP
        let f6p_term = hill_kinetics(1.0, self.k_half_f6p_mM, f6p, self.hill_coefficient);
        let atp_term = michaelis_menten(1.0, self.km_atp_mM, atp);

        self.vmax_mM_per_sec * f6p_term * atp_term
    }

    fn stoichiometry(&self) -> &ReactionStoichiometry {
        &self.stoichiometry
    }
}

/// Aldolase - EC 4.1.2.13
///
/// Fructose-1,6-bisP ⇌ DHAP + Glyceraldehyde-3-P
///
/// Reversible cleavage, Keq = 7.8e-5 M (thermodynamically unfavorable)
/// Reference: Rapoport 1976
pub struct Aldolase {
    pub vmax_f_mM_per_sec: f64,
    pub km_fbp_mM: f64,
    pub keq_mM: f64,
    stoichiometry: ReactionStoichiometry,
    idx_fbp: usize,
    idx_dhap: usize,
    idx_gap: usize,
}

impl Aldolase {
    pub fn new(indices: &MetaboliteIndices) -> Self {
        Self {
            vmax_f_mM_per_sec: 0.3,
            km_fbp_mM: 0.01,
            keq_mM: 7.8e-5,  // Equilibrium constant in mM
            stoichiometry: ReactionStoichiometry::new(
                vec![(indices.fructose_1_6_bisphosphate, 1.0)],
                vec![
                    (indices.dihydroxyacetone_phosphate, 1.0),
                    (indices.glyceraldehyde_3_phosphate, 1.0),
                ],
            ),
            idx_fbp: indices.fructose_1_6_bisphosphate,
            idx_dhap: indices.dihydroxyacetone_phosphate,
            idx_gap: indices.glyceraldehyde_3_phosphate,
        }
    }
}

impl Enzyme for Aldolase {
    fn name(&self) -> &'static str { "Aldolase" }
    fn ec_number(&self) -> &'static str { "4.1.2.13" }

    fn rate(&self, metabolites: &MetabolitePool) -> f64 {
        let fbp = metabolites.get(self.idx_fbp);
        let dhap = metabolites.get(self.idx_dhap);
        let gap = metabolites.get(self.idx_gap);

        // Haldane relationship: Keq = Vf*Km_products / Vr*Km_substrates
        // Simplified reversible kinetics
        let mass_action_ratio = if fbp > 1e-9 {
            (dhap * gap) / (fbp * self.keq_mM)
        } else {
            1.0  // At equilibrium if no substrate
        };

        // Net rate: forward - reverse approaching equilibrium
        let v_forward = self.vmax_f_mM_per_sec * fbp / (self.km_fbp_mM + fbp);
        v_forward * (1.0 - mass_action_ratio).max(-1.0).min(1.0)
    }

    fn stoichiometry(&self) -> &ReactionStoichiometry {
        &self.stoichiometry
    }
}

/// Triose Phosphate Isomerase (TPI) - EC 5.3.1.1
///
/// DHAP ⇌ Glyceraldehyde-3-P
///
/// Very fast equilibrating enzyme, near-diffusion limited
/// Keq = 0.045 ([GAP]/[DHAP])
/// Reference: Rapoport 1976
pub struct TriosePIsomerase {
    pub k_forward: f64,  // Rate constant for mass action
    pub k_reverse: f64,
    pub keq: f64,
    stoichiometry: ReactionStoichiometry,
    idx_dhap: usize,
    idx_gap: usize,
}

impl TriosePIsomerase {
    pub fn new(indices: &MetaboliteIndices) -> Self {
        let keq = 0.045;
        Self {
            k_forward: 10.0,  // Fast equilibrating
            k_reverse: 10.0 / keq,
            keq,
            stoichiometry: ReactionStoichiometry::new(
                vec![(indices.dihydroxyacetone_phosphate, 1.0)],
                vec![(indices.glyceraldehyde_3_phosphate, 1.0)],
            ),
            idx_dhap: indices.dihydroxyacetone_phosphate,
            idx_gap: indices.glyceraldehyde_3_phosphate,
        }
    }
}

impl Enzyme for TriosePIsomerase {
    fn name(&self) -> &'static str { "Triose-P Isomerase" }
    fn ec_number(&self) -> &'static str { "5.3.1.1" }

    fn rate(&self, metabolites: &MetabolitePool) -> f64 {
        let dhap = metabolites.get(self.idx_dhap);
        let gap = metabolites.get(self.idx_gap);

        // Mass action kinetics for fast equilibrating
        self.k_forward * dhap - self.k_reverse * gap
    }

    fn stoichiometry(&self) -> &ReactionStoichiometry {
        &self.stoichiometry
    }
}

/// Glyceraldehyde-3-phosphate Dehydrogenase (GAPDH) - EC 1.2.1.12
///
/// GAP + NAD+ + Pi ⇌ 1,3-BPG + NADH
///
/// Key regulatory point, couples to NADH/NAD+ ratio
/// Reference: Rapoport 1976, Mulquiney-Kuchel 1999
pub struct GAPDH {
    pub vmax_f_mM_per_sec: f64,
    pub vmax_r_mM_per_sec: f64,
    pub km_gap_mM: f64,
    pub km_nad_mM: f64,
    pub km_pi_mM: f64,
    pub km_bpg_mM: f64,
    pub km_nadh_mM: f64,
    stoichiometry: ReactionStoichiometry,
    idx_gap: usize,
    idx_nad: usize,
    idx_pi: usize,
    idx_bpg: usize,
    idx_nadh: usize,
}

impl GAPDH {
    pub fn new(indices: &MetaboliteIndices) -> Self {
        Self {
            vmax_f_mM_per_sec: 1.5,
            vmax_r_mM_per_sec: 0.5,
            km_gap_mM: 0.05,
            km_nad_mM: 0.05,
            km_pi_mM: 1.0,
            km_bpg_mM: 0.001,
            km_nadh_mM: 0.01,
            stoichiometry: ReactionStoichiometry::new(
                vec![
                    (indices.glyceraldehyde_3_phosphate, 1.0),
                    (indices.nad, 1.0),
                    (indices.pi, 1.0),
                ],
                vec![
                    (indices.bisphosphoglycerate_1_3, 1.0),
                    (indices.nadh, 1.0),
                ],
            ),
            idx_gap: indices.glyceraldehyde_3_phosphate,
            idx_nad: indices.nad,
            idx_pi: indices.pi,
            idx_bpg: indices.bisphosphoglycerate_1_3,
            idx_nadh: indices.nadh,
        }
    }
}

impl Enzyme for GAPDH {
    fn name(&self) -> &'static str { "GAPDH" }
    fn ec_number(&self) -> &'static str { "1.2.1.12" }

    fn rate(&self, metabolites: &MetabolitePool) -> f64 {
        let gap = metabolites.get(self.idx_gap);
        let nad = metabolites.get(self.idx_nad);
        let pi = metabolites.get(self.idx_pi);
        let bpg = metabolites.get(self.idx_bpg);
        let nadh = metabolites.get(self.idx_nadh);

        // Simplified reversible terreactant kinetics
        let forward = self.vmax_f_mM_per_sec * (gap / self.km_gap_mM)
            * (nad / self.km_nad_mM)
            * (pi / self.km_pi_mM);
        let reverse = self.vmax_r_mM_per_sec * (bpg / self.km_bpg_mM)
            * (nadh / self.km_nadh_mM);

        let denom = 1.0
            + gap / self.km_gap_mM
            + nad / self.km_nad_mM
            + pi / self.km_pi_mM
            + bpg / self.km_bpg_mM
            + nadh / self.km_nadh_mM;

        (forward - reverse) / denom
    }

    fn stoichiometry(&self) -> &ReactionStoichiometry {
        &self.stoichiometry
    }
}

/// Phosphoglycerate Kinase (PGK) - EC 2.7.2.3
///
/// 1,3-BPG + ADP ⇌ 3-PG + ATP
///
/// Major ATP-producing step, near equilibrium
/// Reference: Rapoport 1976
pub struct PhosphoglycerateKinase {
    pub vmax_f_mM_per_sec: f64,
    pub vmax_r_mM_per_sec: f64,
    pub km_bpg_mM: f64,
    pub km_adp_mM: f64,
    pub km_3pg_mM: f64,
    pub km_atp_mM: f64,
    stoichiometry: ReactionStoichiometry,
    idx_bpg: usize,
    idx_adp: usize,
    idx_3pg: usize,
    idx_atp: usize,
}

impl PhosphoglycerateKinase {
    pub fn new(indices: &MetaboliteIndices) -> Self {
        Self {
            vmax_f_mM_per_sec: 5.0,
            vmax_r_mM_per_sec: 0.1,  // Strongly favors forward
            km_bpg_mM: 0.002,
            km_adp_mM: 0.2,
            km_3pg_mM: 0.5,
            km_atp_mM: 0.5,
            stoichiometry: ReactionStoichiometry::new(
                vec![(indices.bisphosphoglycerate_1_3, 1.0), (indices.adp, 1.0)],
                vec![(indices.phosphoglycerate_3, 1.0), (indices.atp, 1.0)],
            ),
            idx_bpg: indices.bisphosphoglycerate_1_3,
            idx_adp: indices.adp,
            idx_3pg: indices.phosphoglycerate_3,
            idx_atp: indices.atp,
        }
    }
}

impl Enzyme for PhosphoglycerateKinase {
    fn name(&self) -> &'static str { "Phosphoglycerate Kinase" }
    fn ec_number(&self) -> &'static str { "2.7.2.3" }

    fn rate(&self, metabolites: &MetabolitePool) -> f64 {
        let bpg = metabolites.get(self.idx_bpg);
        let adp = metabolites.get(self.idx_adp);
        let pg3 = metabolites.get(self.idx_3pg);
        let atp = metabolites.get(self.idx_atp);

        reversible_ordered_bi_bi(
            self.vmax_f_mM_per_sec,
            self.vmax_r_mM_per_sec,
            self.km_bpg_mM,
            self.km_adp_mM,
            self.km_3pg_mM,
            self.km_atp_mM,
            self.km_bpg_mM,  // Ki approximations
            self.km_atp_mM,
            bpg,
            adp,
            pg3,
            atp,
        )
    }

    fn stoichiometry(&self) -> &ReactionStoichiometry {
        &self.stoichiometry
    }
}

/// Phosphoglycerate Mutase (PGM) - EC 5.4.2.11
///
/// 3-PG ⇌ 2-PG
///
/// Requires 2,3-BPG as cofactor, near equilibrium
/// Keq ≈ 0.18
/// Reference: Rapoport 1976
pub struct PhosphoglycerateMutase {
    pub vmax_f_mM_per_sec: f64,
    pub vmax_r_mM_per_sec: f64,
    pub km_3pg_mM: f64,
    pub km_2pg_mM: f64,
    stoichiometry: ReactionStoichiometry,
    idx_3pg: usize,
    idx_2pg: usize,
}

impl PhosphoglycerateMutase {
    pub fn new(indices: &MetaboliteIndices) -> Self {
        Self {
            vmax_f_mM_per_sec: 2.0,
            vmax_r_mM_per_sec: 11.0,  // Keq ≈ 0.18
            km_3pg_mM: 0.2,
            km_2pg_mM: 0.03,
            stoichiometry: ReactionStoichiometry::new(
                vec![(indices.phosphoglycerate_3, 1.0)],
                vec![(indices.phosphoglycerate_2, 1.0)],
            ),
            idx_3pg: indices.phosphoglycerate_3,
            idx_2pg: indices.phosphoglycerate_2,
        }
    }
}

impl Enzyme for PhosphoglycerateMutase {
    fn name(&self) -> &'static str { "Phosphoglycerate Mutase" }
    fn ec_number(&self) -> &'static str { "5.4.2.11" }

    fn rate(&self, metabolites: &MetabolitePool) -> f64 {
        let pg3 = metabolites.get(self.idx_3pg);
        let pg2 = metabolites.get(self.idx_2pg);

        reversible_mm(
            self.vmax_f_mM_per_sec,
            self.km_3pg_mM,
            pg3,
            self.vmax_r_mM_per_sec,
            self.km_2pg_mM,
            pg2,
        )
    }

    fn stoichiometry(&self) -> &ReactionStoichiometry {
        &self.stoichiometry
    }
}

/// Enolase - EC 4.2.1.11
///
/// 2-PG ⇌ PEP + H2O
///
/// Near equilibrium, Keq ≈ 0.5 (for [PEP]/[2-PG])
/// Reference: Rapoport 1976
pub struct Enolase {
    pub vmax_f_mM_per_sec: f64,
    pub vmax_r_mM_per_sec: f64,
    pub km_2pg_mM: f64,
    pub km_pep_mM: f64,
    stoichiometry: ReactionStoichiometry,
    idx_2pg: usize,
    idx_pep: usize,
}

impl Enolase {
    pub fn new(indices: &MetaboliteIndices) -> Self {
        Self {
            vmax_f_mM_per_sec: 1.5,
            vmax_r_mM_per_sec: 3.0,  // Keq ≈ 0.5
            km_2pg_mM: 0.05,
            km_pep_mM: 0.1,
            stoichiometry: ReactionStoichiometry::new(
                vec![(indices.phosphoglycerate_2, 1.0)],
                vec![(indices.phosphoenolpyruvate, 1.0)],
            ),
            idx_2pg: indices.phosphoglycerate_2,
            idx_pep: indices.phosphoenolpyruvate,
        }
    }
}

impl Enzyme for Enolase {
    fn name(&self) -> &'static str { "Enolase" }
    fn ec_number(&self) -> &'static str { "4.2.1.11" }

    fn rate(&self, metabolites: &MetabolitePool) -> f64 {
        let pg2 = metabolites.get(self.idx_2pg);
        let pep = metabolites.get(self.idx_pep);

        reversible_mm(
            self.vmax_f_mM_per_sec,
            self.km_2pg_mM,
            pg2,
            self.vmax_r_mM_per_sec,
            self.km_pep_mM,
            pep,
        )
    }

    fn stoichiometry(&self) -> &ReactionStoichiometry {
        &self.stoichiometry
    }
}

/// Pyruvate Kinase (PK) - EC 2.7.1.40
///
/// PEP + ADP → Pyruvate + ATP
///
/// Major regulatory enzyme with cooperative kinetics
/// Strongly favors forward direction (essentially irreversible)
/// Reference: Rapoport 1976, Joshi-Palsson 1989
pub struct PyruvateKinase {
    pub vmax_mM_per_sec: f64,
    pub k_half_pep_mM: f64,
    pub km_adp_mM: f64,
    pub hill_coefficient: f64,
    stoichiometry: ReactionStoichiometry,
    idx_pep: usize,
    idx_adp: usize,
}

impl PyruvateKinase {
    pub fn new(indices: &MetaboliteIndices) -> Self {
        Self {
            // Vmax increased from 0.8 to 1.5 mM/s to support ATP homeostasis
            vmax_mM_per_sec: 1.5,
            k_half_pep_mM: 0.2,
            km_adp_mM: 0.3,
            hill_coefficient: 2.0,
            stoichiometry: ReactionStoichiometry::new(
                vec![(indices.phosphoenolpyruvate, 1.0), (indices.adp, 1.0)],
                vec![(indices.pyruvate, 1.0), (indices.atp, 1.0)],
            ),
            idx_pep: indices.phosphoenolpyruvate,
            idx_adp: indices.adp,
        }
    }
}

impl Enzyme for PyruvateKinase {
    fn name(&self) -> &'static str { "Pyruvate Kinase" }
    fn ec_number(&self) -> &'static str { "2.7.1.40" }

    fn rate(&self, metabolites: &MetabolitePool) -> f64 {
        let pep = metabolites.get(self.idx_pep);
        let adp = metabolites.get(self.idx_adp);

        // Hill kinetics for PEP, MM for ADP
        let pep_term = hill_kinetics(1.0, self.k_half_pep_mM, pep, self.hill_coefficient);
        let adp_term = michaelis_menten(1.0, self.km_adp_mM, adp);

        self.vmax_mM_per_sec * pep_term * adp_term
    }

    fn stoichiometry(&self) -> &ReactionStoichiometry {
        &self.stoichiometry
    }
}

/// Lactate Dehydrogenase (LDH) - EC 1.1.1.27
///
/// Pyruvate + NADH ⇌ Lactate + NAD+
///
/// Regenerates NAD+ for GAPDH, strongly favors lactate formation
/// Reference: Rapoport 1976
pub struct LactateDehydrogenase {
    pub vmax_f_mM_per_sec: f64,
    pub vmax_r_mM_per_sec: f64,
    pub km_pyruvate_mM: f64,
    pub km_nadh_mM: f64,
    pub km_lactate_mM: f64,
    pub km_nad_mM: f64,
    stoichiometry: ReactionStoichiometry,
    idx_pyruvate: usize,
    idx_nadh: usize,
    idx_lactate: usize,
    idx_nad: usize,
}

impl LactateDehydrogenase {
    pub fn new(indices: &MetaboliteIndices) -> Self {
        Self {
            // Increased Vmax to ensure NAD+ regeneration for GAPDH
            vmax_f_mM_per_sec: 10.0,
            vmax_r_mM_per_sec: 0.01,  // Keq >> 1
            km_pyruvate_mM: 0.2,
            km_nadh_mM: 0.01,
            km_lactate_mM: 10.0,
            km_nad_mM: 0.1,
            stoichiometry: ReactionStoichiometry::new(
                vec![(indices.pyruvate, 1.0), (indices.nadh, 1.0)],
                vec![(indices.lactate, 1.0), (indices.nad, 1.0)],
            ),
            idx_pyruvate: indices.pyruvate,
            idx_nadh: indices.nadh,
            idx_lactate: indices.lactate,
            idx_nad: indices.nad,
        }
    }
}

impl Enzyme for LactateDehydrogenase {
    fn name(&self) -> &'static str { "Lactate Dehydrogenase" }
    fn ec_number(&self) -> &'static str { "1.1.1.27" }

    fn rate(&self, metabolites: &MetabolitePool) -> f64 {
        let pyr = metabolites.get(self.idx_pyruvate);
        let nadh = metabolites.get(self.idx_nadh);
        let lac = metabolites.get(self.idx_lactate);
        let nad = metabolites.get(self.idx_nad);

        reversible_ordered_bi_bi(
            self.vmax_f_mM_per_sec,
            self.vmax_r_mM_per_sec,
            self.km_pyruvate_mM,
            self.km_nadh_mM,
            self.km_lactate_mM,
            self.km_nad_mM,
            self.km_pyruvate_mM,
            self.km_nad_mM,
            pyr,
            nadh,
            lac,
            nad,
        )
    }

    fn stoichiometry(&self) -> &ReactionStoichiometry {
        &self.stoichiometry
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_test_metabolites() -> MetabolitePool {
        let indices = MetaboliteIndices::default();
        let mut pool = MetabolitePool::new(17);

        // Set physiological concentrations (mM)
        pool.set(indices.glucose, 5.0);
        pool.set(indices.glucose_6_phosphate, 0.05);
        pool.set(indices.fructose_6_phosphate, 0.02);
        pool.set(indices.fructose_1_6_bisphosphate, 0.01);
        pool.set(indices.dihydroxyacetone_phosphate, 0.1);
        pool.set(indices.glyceraldehyde_3_phosphate, 0.005);
        pool.set(indices.bisphosphoglycerate_1_3, 0.001);
        pool.set(indices.phosphoglycerate_3, 0.1);
        pool.set(indices.phosphoglycerate_2, 0.02);
        pool.set(indices.phosphoenolpyruvate, 0.02);
        pool.set(indices.pyruvate, 0.1);
        pool.set(indices.lactate, 1.5);
        pool.set(indices.atp, 2.0);
        pool.set(indices.adp, 0.25);
        pool.set(indices.nad, 0.07);
        pool.set(indices.nadh, 0.03);
        pool.set(indices.pi, 1.0);

        pool
    }

    #[test]
    fn test_hexokinase_rate() {
        let indices = MetaboliteIndices::default();
        let hk = Hexokinase::new(&indices);
        let metabolites = create_test_metabolites();

        let rate = hk.rate(&metabolites);
        assert!(rate > 0.0, "Hexokinase should have positive rate");
        assert!(rate < hk.vmax_mM_per_sec, "Rate should be below Vmax");
    }

    #[test]
    fn test_pfk_cooperativity() {
        let indices = MetaboliteIndices::default();
        let pfk = Phosphofructokinase::new(&indices);
        let mut metabolites = create_test_metabolites();

        // Rate at low F6P
        metabolites.set(indices.fructose_6_phosphate, 0.01);
        let rate_low = pfk.rate(&metabolites);

        // Rate at K_half
        metabolites.set(indices.fructose_6_phosphate, pfk.k_half_f6p_mM);
        let rate_half = pfk.rate(&metabolites);

        // Rate at high F6P
        metabolites.set(indices.fructose_6_phosphate, 1.0);
        let rate_high = pfk.rate(&metabolites);

        assert!(rate_low < rate_half);
        assert!(rate_half < rate_high);

        // Cooperative: rate increase should be sigmoidal
        // Ratio from low to half should be greater than linear
    }

    #[test]
    fn test_glycolysis_solver() {
        let solver = GlycolysisSolver::new();
        let metabolites = create_test_metabolites();

        let mut dydt = vec![0.0; 17];
        solver.compute_derivatives(&metabolites, &mut dydt);

        // ATP should have both production and consumption
        // Net depends on balance of kinases vs ATP-consuming reactions
        assert!(dydt[solver.indices.atp].is_finite());

        // Glucose should be consumed
        assert!(dydt[solver.indices.glucose] <= 0.0);

        // Lactate should be produced
        assert!(dydt[solver.indices.lactate] >= 0.0);
    }

    #[test]
    fn test_all_enzymes_return_finite_rates() {
        let solver = GlycolysisSolver::new();
        let metabolites = create_test_metabolites();

        let rates = solver.get_rates(&metabolites);
        for (name, rate) in rates {
            assert!(rate.is_finite(), "Enzyme {} returned non-finite rate: {}", name, rate);
        }
    }
}
