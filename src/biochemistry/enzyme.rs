//! Enzyme kinetics framework for metabolic modeling.
//!
//! Provides kinetic rate equations commonly used in systems biology:
//! - Michaelis-Menten kinetics
//! - Hill kinetics (cooperative binding)
//! - Reversible Michaelis-Menten
//! - Ordered bi-bi mechanism
//!
//! References:
//! - Cornish-Bowden A. Fundamentals of Enzyme Kinetics. 4th ed. Wiley-Blackwell, 2012
//! - Segel IH. Enzyme Kinetics. Wiley-Interscience, 1993

use super::MetabolitePool;

/// Stoichiometry of a metabolic reaction
#[derive(Debug, Clone)]
pub struct ReactionStoichiometry {
    /// Substrates consumed (metabolite index, stoichiometric coefficient)
    pub substrates: Vec<(usize, f64)>,
    /// Products produced (metabolite index, stoichiometric coefficient)
    pub products: Vec<(usize, f64)>,
}

impl ReactionStoichiometry {
    /// Create a new stoichiometry definition
    pub fn new(substrates: Vec<(usize, f64)>, products: Vec<(usize, f64)>) -> Self {
        Self { substrates, products }
    }

    /// Apply the reaction to a derivatives vector given a reaction rate
    pub fn apply(&self, dydt: &mut [f64], rate_mM_per_sec: f64) {
        // Consume substrates
        for &(idx, coeff) in &self.substrates {
            if idx < dydt.len() {
                dydt[idx] -= coeff * rate_mM_per_sec;
            }
        }
        // Produce products
        for &(idx, coeff) in &self.products {
            if idx < dydt.len() {
                dydt[idx] += coeff * rate_mM_per_sec;
            }
        }
    }
}

/// Enzyme trait for metabolic reactions
pub trait Enzyme {
    /// Enzyme name (e.g., "Hexokinase")
    fn name(&self) -> &'static str;

    /// EC number (e.g., "2.7.1.1" for hexokinase)
    fn ec_number(&self) -> &'static str;

    /// Compute reaction rate in mM/s given current metabolite concentrations
    fn rate(&self, metabolites: &MetabolitePool) -> f64;

    /// Get the stoichiometry of the reaction
    fn stoichiometry(&self) -> &ReactionStoichiometry;

    /// Apply enzyme to derivatives vector
    fn apply_to_derivatives(&self, metabolites: &MetabolitePool, dydt: &mut [f64]) {
        let rate = self.rate(metabolites);
        self.stoichiometry().apply(dydt, rate);
    }
}

// ============================================================================
// Standard Kinetic Rate Equations
// ============================================================================

/// Simple Michaelis-Menten kinetics
///
/// v = Vmax * [S] / (Km + [S])
///
/// # Arguments
/// * `vmax_mM_per_sec` - Maximum reaction velocity (mM/s)
/// * `km_mM` - Michaelis constant (mM)
/// * `s_mM` - Substrate concentration (mM)
///
/// # Reference
/// Michaelis L, Menten ML. Biochemische Zeitschrift. 1913;49:333-369
#[inline]
pub fn michaelis_menten(vmax_mM_per_sec: f64, km_mM: f64, s_mM: f64) -> f64 {
    if s_mM <= 0.0 {
        return 0.0;
    }
    vmax_mM_per_sec * s_mM / (km_mM + s_mM)
}

/// Hill kinetics for cooperative enzymes
///
/// v = Vmax * [S]^n / (K0.5^n + [S]^n)
///
/// # Arguments
/// * `vmax_mM_per_sec` - Maximum reaction velocity (mM/s)
/// * `k_half_mM` - Half-saturation constant (mM)
/// * `s_mM` - Substrate concentration (mM)
/// * `n` - Hill coefficient (n > 1 for positive cooperativity)
///
/// # Reference
/// Hill AV. Journal of Physiology. 1910;40:iv-vii
#[inline]
pub fn hill_kinetics(vmax_mM_per_sec: f64, k_half_mM: f64, s_mM: f64, n: f64) -> f64 {
    if s_mM <= 0.0 {
        return 0.0;
    }
    let s_n = s_mM.powf(n);
    let k_n = k_half_mM.powf(n);
    vmax_mM_per_sec * s_n / (k_n + s_n)
}

/// Reversible Michaelis-Menten kinetics
///
/// v = (Vmax_f * [S]/Km_s - Vmax_r * [P]/Km_p) / (1 + [S]/Km_s + [P]/Km_p)
///
/// This equation satisfies the Haldane relationship at equilibrium.
///
/// # Arguments
/// * `vmax_f_mM_per_sec` - Forward maximum velocity (mM/s)
/// * `km_s_mM` - Michaelis constant for substrate (mM)
/// * `s_mM` - Substrate concentration (mM)
/// * `vmax_r_mM_per_sec` - Reverse maximum velocity (mM/s)
/// * `km_p_mM` - Michaelis constant for product (mM)
/// * `p_mM` - Product concentration (mM)
///
/// # Reference
/// Cornish-Bowden A. Fundamentals of Enzyme Kinetics. 4th ed. Chapter 8
#[inline]
pub fn reversible_mm(
    vmax_f_mM_per_sec: f64,
    km_s_mM: f64,
    s_mM: f64,
    vmax_r_mM_per_sec: f64,
    km_p_mM: f64,
    p_mM: f64,
) -> f64 {
    let s_term = s_mM / km_s_mM;
    let p_term = p_mM / km_p_mM;
    let numerator = vmax_f_mM_per_sec * s_term - vmax_r_mM_per_sec * p_term;
    let denominator = 1.0 + s_term + p_term;

    if denominator <= 0.0 {
        return 0.0;
    }
    numerator / denominator
}

/// Ordered bi-bi mechanism (two substrates, two products)
///
/// For reactions like: A + B → P + Q
/// where substrate A binds first, then B
///
/// v = Vmax * [A] * [B] / (Ki_a * Km_b + Km_b * [A] + Km_a * [B] + [A] * [B])
///
/// # Arguments
/// * `vmax_mM_per_sec` - Maximum velocity (mM/s)
/// * `km_a_mM` - Michaelis constant for substrate A (mM)
/// * `km_b_mM` - Michaelis constant for substrate B (mM)
/// * `ki_a_mM` - Inhibition constant for A (mM)
/// * `a_mM` - Concentration of substrate A (mM)
/// * `b_mM` - Concentration of substrate B (mM)
///
/// # Reference
/// Cleland WW. Biochim Biophys Acta. 1963;67:104-137
#[inline]
pub fn ordered_bi_bi(
    vmax_mM_per_sec: f64,
    km_a_mM: f64,
    km_b_mM: f64,
    ki_a_mM: f64,
    a_mM: f64,
    b_mM: f64,
) -> f64 {
    if a_mM <= 0.0 || b_mM <= 0.0 {
        return 0.0;
    }
    let numerator = vmax_mM_per_sec * a_mM * b_mM;
    let denominator = ki_a_mM * km_b_mM + km_b_mM * a_mM + km_a_mM * b_mM + a_mM * b_mM;

    if denominator <= 0.0 {
        return 0.0;
    }
    numerator / denominator
}

/// Reversible ordered bi-bi mechanism
///
/// Full reversible rate equation for reactions like: A + B ⇌ P + Q
///
/// # Arguments
/// * `vmax_f` - Forward maximum velocity (mM/s)
/// * `vmax_r` - Reverse maximum velocity (mM/s)
/// * `km_a`, `km_b`, `km_p`, `km_q` - Michaelis constants (mM)
/// * `ki_a`, `ki_q` - Inhibition constants (mM)
/// * `a`, `b`, `p`, `q` - Substrate/product concentrations (mM)
#[inline]
#[allow(clippy::too_many_arguments)]
pub fn reversible_ordered_bi_bi(
    vmax_f: f64,
    vmax_r: f64,
    _km_a: f64,
    km_b: f64,
    km_p: f64,
    _km_q: f64,
    ki_a: f64,
    ki_q: f64,
    a: f64,
    b: f64,
    p: f64,
    q: f64,
) -> f64 {
    // Numerator: forward - reverse
    let numerator = vmax_f * a * b / (ki_a * km_b) - vmax_r * p * q / (ki_q * km_p);

    // Denominator (simplified - full expression is complex)
    let denom = 1.0
        + a / ki_a
        + b / km_b
        + a * b / (ki_a * km_b)
        + p / km_p
        + q / ki_q
        + p * q / (ki_q * km_p);

    if denom <= 0.0 {
        return 0.0;
    }
    numerator / denom
}

/// Michaelis-Menten with competitive inhibition
///
/// v = Vmax * [S] / (Km * (1 + [I]/Ki) + [S])
///
/// # Arguments
/// * `vmax_mM_per_sec` - Maximum velocity (mM/s)
/// * `km_mM` - Michaelis constant (mM)
/// * `s_mM` - Substrate concentration (mM)
/// * `ki_mM` - Inhibition constant (mM)
/// * `i_mM` - Inhibitor concentration (mM)
#[inline]
pub fn competitive_inhibition(
    vmax_mM_per_sec: f64,
    km_mM: f64,
    s_mM: f64,
    ki_mM: f64,
    i_mM: f64,
) -> f64 {
    if s_mM <= 0.0 {
        return 0.0;
    }
    let km_apparent = km_mM * (1.0 + i_mM / ki_mM);
    vmax_mM_per_sec * s_mM / (km_apparent + s_mM)
}

/// Michaelis-Menten with product inhibition
///
/// Common pattern in metabolic pathways where product accumulation
/// slows the reaction.
///
/// v = Vmax * [S] / (Km + [S] + [S]*[P]/Kip)
#[inline]
pub fn product_inhibition(
    vmax_mM_per_sec: f64,
    km_mM: f64,
    s_mM: f64,
    kip_mM: f64,
    p_mM: f64,
) -> f64 {
    if s_mM <= 0.0 {
        return 0.0;
    }
    let denominator = km_mM + s_mM * (1.0 + p_mM / kip_mM);
    vmax_mM_per_sec * s_mM / denominator
}

/// Calculate equilibrium constant from forward and reverse Vmax and Km
///
/// Keq = (Vmax_f * Km_p) / (Vmax_r * Km_s)
///
/// This is the Haldane relationship.
#[inline]
pub fn haldane_keq(vmax_f: f64, km_s: f64, vmax_r: f64, km_p: f64) -> f64 {
    (vmax_f * km_p) / (vmax_r * km_s)
}

/// Mass action kinetics for near-equilibrium reactions
///
/// v = k_f * [S] - k_r * [P]
///
/// Used for fast equilibrating reactions like triose phosphate isomerase
#[inline]
pub fn mass_action(k_f: f64, s_mM: f64, k_r: f64, p_mM: f64) -> f64 {
    k_f * s_mM.max(0.0) - k_r * p_mM.max(0.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_michaelis_menten() {
        // At [S] = Km, rate should be Vmax/2
        let vmax = 1.0;
        let km = 0.1;
        let rate = michaelis_menten(vmax, km, km);
        assert!((rate - 0.5).abs() < 1e-10);

        // At very high [S], rate approaches Vmax
        let rate_high = michaelis_menten(vmax, km, 100.0 * km);
        assert!((rate_high - vmax).abs() < 0.01);

        // At zero substrate, rate is zero
        let rate_zero = michaelis_menten(vmax, km, 0.0);
        assert_eq!(rate_zero, 0.0);
    }

    #[test]
    fn test_hill_kinetics() {
        let vmax = 1.0;
        let k_half = 0.1;

        // At [S] = K_half, rate should be Vmax/2 regardless of n
        let rate_n1 = hill_kinetics(vmax, k_half, k_half, 1.0);
        let rate_n2 = hill_kinetics(vmax, k_half, k_half, 2.0);
        let rate_n4 = hill_kinetics(vmax, k_half, k_half, 4.0);

        assert!((rate_n1 - 0.5).abs() < 1e-10);
        assert!((rate_n2 - 0.5).abs() < 1e-10);
        assert!((rate_n4 - 0.5).abs() < 1e-10);

        // With n=1, should equal Michaelis-Menten
        let mm_rate = michaelis_menten(vmax, k_half, 0.05);
        let hill_rate = hill_kinetics(vmax, k_half, 0.05, 1.0);
        assert!((mm_rate - hill_rate).abs() < 1e-10);
    }

    #[test]
    fn test_reversible_mm_equilibrium() {
        // At equilibrium, forward rate = reverse rate, so net rate = 0
        // Keq = [P]/[S] = (Vmax_f * Km_p) / (Vmax_r * Km_s)
        let vmax_f = 1.0;
        let vmax_r = 0.5;
        let km_s = 0.1;
        let km_p = 0.2;

        let keq = haldane_keq(vmax_f, km_s, vmax_r, km_p);
        // Keq = (1.0 * 0.2) / (0.5 * 0.1) = 0.2 / 0.05 = 4.0

        // At equilibrium: [P]/[S] = Keq
        let s = 0.1;
        let p = s * keq; // p = 0.4

        let rate = reversible_mm(vmax_f, km_s, s, vmax_r, km_p, p);
        assert!(rate.abs() < 1e-10, "Rate at equilibrium should be ~0, got {}", rate);
    }

    #[test]
    fn test_ordered_bi_bi() {
        let vmax = 1.0;
        let km_a = 0.1;
        let km_b = 0.2;
        let ki_a = 0.05;

        // Zero substrates gives zero rate
        assert_eq!(ordered_bi_bi(vmax, km_a, km_b, ki_a, 0.0, 1.0), 0.0);
        assert_eq!(ordered_bi_bi(vmax, km_a, km_b, ki_a, 1.0, 0.0), 0.0);

        // Rate should be positive with positive substrates
        let rate = ordered_bi_bi(vmax, km_a, km_b, ki_a, 1.0, 1.0);
        assert!(rate > 0.0);
        assert!(rate <= vmax);
    }

    #[test]
    fn test_competitive_inhibition() {
        let vmax = 1.0;
        let km = 0.1;
        let ki = 0.1;
        let s = 0.1;

        // No inhibitor: standard MM
        let rate_no_inhib = competitive_inhibition(vmax, km, s, ki, 0.0);
        let rate_mm = michaelis_menten(vmax, km, s);
        assert!((rate_no_inhib - rate_mm).abs() < 1e-10);

        // With inhibitor, rate should decrease
        let rate_with_inhib = competitive_inhibition(vmax, km, s, ki, 0.1);
        assert!(rate_with_inhib < rate_no_inhib);
    }

    #[test]
    fn test_stoichiometry_apply() {
        // Simple reaction: A → B
        let stoich = ReactionStoichiometry::new(
            vec![(0, 1.0)],  // Consume 1 A
            vec![(1, 1.0)],  // Produce 1 B
        );

        let mut dydt = vec![0.0, 0.0];
        stoich.apply(&mut dydt, 0.1); // Rate = 0.1 mM/s

        assert!((dydt[0] - (-0.1)).abs() < 1e-10); // A decreases
        assert!((dydt[1] - 0.1).abs() < 1e-10);    // B increases
    }
}
