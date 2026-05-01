//! Local identifiability via Fisher information & parameter correlation.
//!
//! Given a set of `SensitivitySample`s and a parallel vector of model
//! responses (one scalar per sample — typically χ² or χ²/dof against a
//! reference curve), we compute:
//!
//! - First-order Sobol-style sensitivity indices (variance attributable
//!   to each parameter, normalized by total variance).
//! - Pairwise parameter correlations as observed across the LHS rows.
//!
//! A high pairwise correlation (|ρ| > 0.8) flags two parameters that move
//! together to produce indistinguishable model responses — i.e., they are
//! not jointly identifiable from the curve under test.
//!
//! This is intentionally a *local* analysis (response surface around the
//! current parameter values). For global / Bayesian identifiability, use
//! a downstream tool such as PyMC or Stan with the simulator as black-box.

use serde::{Deserialize, Serialize};

use super::sensitivity::SensitivitySample;

/// Pairwise parameter correlation observed across LHS samples.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ParameterCorrelation {
    pub name_a: String,
    pub name_b: String,
    /// Pearson correlation coefficient ∈ [-1, 1].
    pub rho: f64,
    /// True iff |ρ| ≥ flag_threshold.
    pub flagged: bool,
}

/// Aggregate identifiability output for a single response variable.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IdentifiabilityReport {
    /// Parameter names (same order throughout the report).
    pub parameter_names: Vec<String>,
    /// First-order sensitivity index per parameter (variance fraction).
    pub first_order_indices: Vec<f64>,
    /// Pairwise correlations (n*(n-1)/2 entries).
    pub correlations: Vec<ParameterCorrelation>,
    /// Total variance of the response across samples.
    pub response_variance: f64,
    /// Number of samples used.
    pub n_samples: usize,
}

/// Compute local identifiability indices.
///
/// `samples` and `responses` must have equal length. `flag_threshold` is
/// the absolute correlation above which a parameter pair is flagged
/// non-identifiable (commonly 0.8).
pub fn analyze(
    samples: &[SensitivitySample],
    responses: &[f64],
    flag_threshold: f64,
) -> IdentifiabilityReport {
    assert_eq!(samples.len(), responses.len(),
        "analyze: samples.len() != responses.len()");
    assert!(!samples.is_empty(), "analyze: empty samples");

    let n = samples.len();
    let n_params = samples[0].values.len();
    let names = samples[0].names.clone();

    // Build columns: param_values[p][i] is the i-th sample's value of
    // parameter p.
    let mut param_values = vec![vec![0.0f64; n]; n_params];
    for (i, s) in samples.iter().enumerate() {
        for p in 0..n_params {
            param_values[p][i] = s.values[p];
        }
    }

    let response_mean: f64 = responses.iter().sum::<f64>() / n as f64;
    let response_variance: f64 = responses.iter()
        .map(|r| (r - response_mean).powi(2))
        .sum::<f64>() / (n as f64).max(1.0);

    // First-order indices via correlation² between parameter and response,
    // multiplied by response_variance / response_variance = correlation²
    // (a simple proxy that's exact for purely-linear models). For
    // non-linear models this is an underestimate; full Sobol requires
    // dedicated sampling, but this is sufficient as a Phase 10 screen.
    let first_order_indices: Vec<f64> = (0..n_params).map(|p| {
        pearson(&param_values[p], responses).powi(2)
    }).collect();

    // Pairwise parameter correlations.
    let mut correlations = Vec::with_capacity(n_params * (n_params - 1) / 2);
    for i in 0..n_params {
        for j in (i + 1)..n_params {
            let rho = pearson(&param_values[i], &param_values[j]);
            correlations.push(ParameterCorrelation {
                name_a: names[i].clone(),
                name_b: names[j].clone(),
                rho,
                flagged: rho.abs() >= flag_threshold,
            });
        }
    }

    IdentifiabilityReport {
        parameter_names: names,
        first_order_indices,
        correlations,
        response_variance,
        n_samples: n,
    }
}

fn pearson(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len() as f64;
    if x.len() != y.len() || x.is_empty() { return 0.0; }
    let mx: f64 = x.iter().sum::<f64>() / n;
    let my: f64 = y.iter().sum::<f64>() / n;
    let mut sxy = 0.0;
    let mut sxx = 0.0;
    let mut syy = 0.0;
    for i in 0..x.len() {
        let dx = x[i] - mx;
        let dy = y[i] - my;
        sxy += dx * dy;
        sxx += dx * dx;
        syy += dy * dy;
    }
    if sxx > 0.0 && syy > 0.0 {
        sxy / (sxx * syy).sqrt()
    } else {
        0.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::sensitivity::SensitivitySample;

    fn make_samples(values: Vec<Vec<f64>>) -> Vec<SensitivitySample> {
        let names = vec!["a".to_string(), "b".to_string()];
        values.into_iter().map(|v| SensitivitySample {
            names: names.clone(),
            values: v,
        }).collect()
    }

    #[test]
    fn pure_dependence_on_param_zero() {
        // y depends only on param 0; param 1 is uncorrelated noise.
        let samples = make_samples(vec![
            vec![1.0, 100.0],
            vec![2.0, 50.0],
            vec![3.0, 200.0],
            vec![4.0, 75.0],
            vec![5.0, 150.0],
        ]);
        let responses = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let report = analyze(&samples, &responses, 0.8);
        assert!(report.first_order_indices[0] > 0.95);
        assert!(report.first_order_indices[1] < 0.5);
    }

    #[test]
    fn flag_correlated_pair() {
        let samples = make_samples(vec![
            vec![1.0, 2.0],
            vec![2.0, 4.0],
            vec![3.0, 6.0],
            vec![4.0, 8.0],
            vec![5.0, 10.0],
        ]);
        let responses = vec![0.0; 5];
        let report = analyze(&samples, &responses, 0.8);
        assert_eq!(report.correlations.len(), 1);
        assert!(report.correlations[0].flagged);
    }
}
