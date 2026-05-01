//! Quantitative fit metrics: R², RMSE, NRMSE, χ²/dof, AIC.
//!
//! These metrics are computed from a `ValidationCurve` (observed) and a
//! parallel vector of model predictions at the same x values.
//!
//! χ²/dof is the primary metric for pass/fail decisions because it correctly
//! handles heteroscedastic measurement uncertainty. For points with very small
//! `sigma`, χ² grows quickly when the model misses; for noisy points, it
//! tolerates more deviation.

use serde::{Deserialize, Serialize};

use super::reference_curve::ValidationCurve;

/// Fit metrics computed from observed and predicted vectors.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FitMetrics {
    /// Number of data points.
    pub n_points: usize,
    /// Number of free parameters used to compute the prediction (for dof).
    pub n_free_params: usize,
    /// Coefficient of determination R² (1.0 = perfect, < 0 = worse than mean).
    pub r_squared: f64,
    /// Root mean square error in original y units.
    pub rmse: f64,
    /// Normalized RMSE (RMSE / range of observed y).
    pub nrmse: f64,
    /// χ² statistic = Σ ((y_obs - y_pred) / sigma)².
    pub chi_squared: f64,
    /// χ² per degree of freedom = χ² / (n_points - n_free_params).
    /// 1.0 = matches quoted error bars; > 2 = poor fit; < 0.5 = overfit.
    pub chi2_per_dof: f64,
    /// Akaike Information Criterion (lower is better).
    pub aic: f64,
    /// Mean residual (signed bias).
    pub mean_residual: f64,
    /// Maximum absolute residual.
    pub max_abs_residual: f64,
}

/// Compute fit metrics from a reference curve and a parallel prediction vector.
///
/// `n_free_params` is the number of model parameters that were fit *to this
/// curve*. For pure prediction (no fitting), pass 0.
pub fn compute_metrics(
    curve: &ValidationCurve,
    predicted: &[f64],
    n_free_params: usize,
) -> FitMetrics {
    assert_eq!(curve.len(), predicted.len(),
        "compute_metrics: predicted.len() != curve.len()");

    let n = curve.len();
    let y_mean: f64 = curve.y.iter().sum::<f64>() / n as f64;

    let mut ss_res: f64 = 0.0;
    let mut ss_tot: f64 = 0.0;
    let mut chi_squared: f64 = 0.0;
    let mut sum_residual: f64 = 0.0;
    let mut max_abs_residual: f64 = 0.0;

    for i in 0..n {
        let r = curve.y[i] - predicted[i];
        ss_res += r * r;
        ss_tot += (curve.y[i] - y_mean).powi(2);
        let s = curve.sigma[i].max(1e-12);
        chi_squared += (r / s).powi(2);
        sum_residual += r;
        max_abs_residual = max_abs_residual.max(r.abs());
    }

    let r_squared = if ss_tot > 0.0 { 1.0 - ss_res / ss_tot } else { f64::NAN };
    let rmse = (ss_res / n as f64).sqrt();
    let y_range = curve.y.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
                - curve.y.iter().cloned().fold(f64::INFINITY, f64::min);
    let nrmse = if y_range > 0.0 { rmse / y_range } else { f64::NAN };
    let dof = (n.saturating_sub(n_free_params)).max(1) as f64;
    let chi2_per_dof = chi_squared / dof;
    // AIC for Gaussian likelihood: 2k + n*ln(RSS/n) — useful for model comparison.
    let aic = if ss_res > 0.0 {
        2.0 * n_free_params as f64 + n as f64 * (ss_res / n as f64).ln()
    } else {
        f64::NEG_INFINITY
    };
    let mean_residual = sum_residual / n as f64;

    FitMetrics {
        n_points: n,
        n_free_params,
        r_squared,
        rmse,
        nrmse,
        chi_squared,
        chi2_per_dof,
        aic,
        mean_residual,
        max_abs_residual,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::validation::reference_curve::CurveMetadata;

    fn metadata() -> CurveMetadata {
        CurveMetadata {
            id: "test".into(),
            citation: "Test 2026".into(),
            doi: "".into(),
            figure: "".into(),
            digitization: "synthetic".into(),
            x_label: "x".into(),
            y_label: "y".into(),
            notes: "".into(),
        }
    }

    #[test]
    fn perfect_fit_gives_unit_chi2_with_unit_sigma() {
        let curve = ValidationCurve::new(
            metadata(),
            vec![1.0, 2.0, 3.0, 4.0],
            vec![10.0, 20.0, 30.0, 40.0],
            vec![1.0, 1.0, 1.0, 1.0],
        );
        let predicted = vec![10.0, 20.0, 30.0, 40.0];
        let m = compute_metrics(&curve, &predicted, 0);
        assert!(m.r_squared > 0.999);
        assert!(m.rmse < 1e-9);
        assert!(m.chi_squared < 1e-9);
        assert!(m.chi2_per_dof < 1e-9);
    }

    #[test]
    fn one_sigma_off_gives_chi2_around_one() {
        // 4 points, each off by exactly 1 sigma. χ² = 4, dof = 4, χ²/dof = 1.
        let curve = ValidationCurve::new(
            metadata(),
            vec![1.0, 2.0, 3.0, 4.0],
            vec![10.0, 20.0, 30.0, 40.0],
            vec![1.0, 1.0, 1.0, 1.0],
        );
        let predicted = vec![11.0, 21.0, 31.0, 41.0];
        let m = compute_metrics(&curve, &predicted, 0);
        assert!((m.chi_squared - 4.0).abs() < 1e-9);
        assert!((m.chi2_per_dof - 1.0).abs() < 1e-9);
    }

    #[test]
    fn r_squared_zero_when_predicting_mean() {
        let curve = ValidationCurve::new(
            metadata(),
            vec![1.0, 2.0, 3.0, 4.0],
            vec![10.0, 20.0, 30.0, 40.0],
            vec![1.0, 1.0, 1.0, 1.0],
        );
        let mean = 25.0;
        let predicted = vec![mean; 4];
        let m = compute_metrics(&curve, &predicted, 0);
        assert!(m.r_squared.abs() < 1e-9);
    }
}
