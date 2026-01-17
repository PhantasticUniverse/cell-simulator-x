//! Biochemistry state data structures.
//!
//! Represents the metabolic and chemical state of the red blood cell,
//! including ion concentrations, metabolites, and hemoglobin states.

/// Complete biochemistry state of the cell
#[derive(Debug, Clone)]
pub struct BiochemistryState {
    /// Ion concentrations
    pub ions: IonState,
    /// Metabolite concentrations
    pub metabolites: MetaboliteState,
    /// Hemoglobin state
    pub hemoglobin: HemoglobinState,
    /// Intracellular pH
    /// Reference: Normal RBC pH ~7.2
    /// Source: Jacobs & Stewart, J Gen Physiol 1947
    pub ph: f32,
    /// 2,3-DPG concentration (mM)
    /// Reference: 4.5-5.5 mM in normal RBCs
    /// Source: Benesch & Benesch, Nature 1969
    pub dpg_concentration_mM: f32,
}

impl Default for BiochemistryState {
    fn default() -> Self {
        Self {
            ions: IonState::default(),
            metabolites: MetaboliteState::default(),
            hemoglobin: HemoglobinState::default(),
            ph: 7.2,
            dpg_concentration_mM: 5.0,
        }
    }
}

/// Ion concentrations within the cell (mM)
#[derive(Debug, Clone)]
pub struct IonState {
    /// Potassium concentration (mM)
    /// Reference: ~140 mM intracellular
    /// Source: Bernstein, Physiol Rev 1954
    pub potassium_mM: f32,
    /// Sodium concentration (mM)
    /// Reference: ~10 mM intracellular
    /// Source: Bernstein, Physiol Rev 1954
    pub sodium_mM: f32,
    /// Calcium concentration (mM)
    /// Reference: ~0.0001 mM (100 nM) intracellular free Ca²⁺
    /// Source: Schatzmann, J Physiol 1966
    pub calcium_mM: f32,
    /// Chloride concentration (mM)
    /// Reference: ~80 mM intracellular
    /// Source: Dalmark, J Gen Physiol 1975
    pub chloride_mM: f32,
    /// Magnesium concentration (mM)
    /// Reference: ~2-3 mM intracellular
    /// Source: Flatman & Lew, J Physiol 1980
    pub magnesium_mM: f32,
}

impl Default for IonState {
    fn default() -> Self {
        Self {
            potassium_mM: 140.0,
            sodium_mM: 10.0,
            calcium_mM: 0.0001,
            chloride_mM: 80.0,
            magnesium_mM: 2.5,
        }
    }
}

/// Key metabolite concentrations (mM)
#[derive(Debug, Clone)]
pub struct MetaboliteState {
    /// ATP concentration (mM)
    /// Reference: 1.5-2.5 mM in normal RBCs
    /// Source: Minakami & Yoshikawa, Biochem Biophys Res Commun 1966
    pub atp_mM: f32,
    /// ADP concentration (mM)
    /// Reference: ~0.2-0.3 mM
    /// Source: Minakami & Yoshikawa, 1966
    pub adp_mM: f32,
    /// NADH concentration (mM)
    /// Reference: ~0.03 mM
    /// Source: Zerez et al., Blood 1987
    pub nadh_mM: f32,
    /// NAD+ concentration (mM)
    /// Reference: ~0.07 mM
    /// Source: Zerez et al., Blood 1987
    pub nad_mM: f32,
    /// Glucose concentration (mM)
    /// Reference: ~4-5 mM (equilibrated with plasma)
    /// Source: Widdas, J Physiol 1954
    pub glucose_mM: f32,
    /// Lactate concentration (mM)
    /// Reference: ~1-2 mM
    /// Source: Jacobasch & Rapoport, Mol Aspects Med 1996
    pub lactate_mM: f32,
}

impl Default for MetaboliteState {
    fn default() -> Self {
        Self {
            atp_mM: 2.0,
            adp_mM: 0.25,
            nadh_mM: 0.03,
            nad_mM: 0.07,
            glucose_mM: 5.0,
            lactate_mM: 1.5,
        }
    }
}

/// Hemoglobin oxygenation and conformational state
#[derive(Debug, Clone)]
pub struct HemoglobinState {
    /// Total hemoglobin concentration (mM of tetramer)
    /// Reference: ~5 mM (MCHC ~33 g/dL = ~5 mM tetramer)
    /// Source: Bunn & Forget, Hemoglobin: Molecular, Genetic and Clinical Aspects, 1986
    pub total_concentration_mM: f32,
    /// Fractional saturation (0.0 - 1.0)
    /// Reference: Variable, depends on pO2
    pub saturation: f32,
    /// Fraction in T (tense/deoxy) state
    pub t_state_fraction: f32,
    /// Fraction in R (relaxed/oxy) state
    pub r_state_fraction: f32,
    /// Methemoglobin fraction
    /// Reference: <1% in normal RBCs
    /// Source: Jaffe, Clin Haematol 1981
    pub met_hb_fraction: f32,
}

impl Default for HemoglobinState {
    fn default() -> Self {
        Self {
            total_concentration_mM: 5.0,
            saturation: 0.97, // Arterial blood
            t_state_fraction: 0.03,
            r_state_fraction: 0.97,
            met_hb_fraction: 0.005,
        }
    }
}
