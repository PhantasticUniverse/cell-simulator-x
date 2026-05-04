//! Latin-hypercube parameter sensitivity sweep.
//!
//! Given a set of named parameter ranges, this generates a Latin-hypercube
//! sample of size N and returns a vector of `SensitivitySample`s with one
//! entry per LHS row. Downstream code evaluates the model at each sample
//! and collects responses for variance-based sensitivity analysis.
//!
//! LHS construction follows the standard recipe (McKay, Beckman, Conover
//! 1979): partition each dimension into N equal-probability intervals,
//! permute the intervals independently per dimension, then jitter within
//! each interval.

use rand::seq::SliceRandom;
use rand::{Rng, SeedableRng};
use rand::rngs::StdRng;
use serde::{Deserialize, Serialize};

/// A named parameter with a [low, high] range. Ranges are uniform unless
/// `log_scale` is true (then sampled in log10 space and exponentiated).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ParameterRange {
    pub name: String,
    pub low: f64,
    pub high: f64,
    pub log_scale: bool,
    /// Optional source for the range (e.g. literature span).
    pub source: String,
}

/// One Latin-hypercube row: the value drawn for each parameter.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SensitivitySample {
    /// Names parallel to `values` (same order as the input ranges).
    pub names: Vec<String>,
    /// Sampled parameter values.
    pub values: Vec<f64>,
}

/// A complete LHS sweep.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SensitivitySweep {
    pub samples: Vec<SensitivitySample>,
    pub ranges: Vec<ParameterRange>,
    pub seed: u64,
}

impl SensitivitySweep {
    /// Generate an LHS of size `n` over the given parameter ranges.
    /// Reproducible via `seed`.
    pub fn generate(ranges: &[ParameterRange], n: usize, seed: u64) -> Self {
        assert!(n > 0, "SensitivitySweep: n must be > 0");
        assert!(!ranges.is_empty(), "SensitivitySweep: ranges must be non-empty");

        let mut rng = StdRng::seed_from_u64(seed);
        let n_dims = ranges.len();

        // For each dimension, build a permutation of [0..n] and jitter within.
        let mut permutations: Vec<Vec<usize>> = (0..n_dims)
            .map(|_| {
                let mut v: Vec<usize> = (0..n).collect();
                v.shuffle(&mut rng);
                v
            })
            .collect();

        let names: Vec<String> = ranges.iter().map(|r| r.name.clone()).collect();
        let mut samples = Vec::with_capacity(n);

        for row in 0..n {
            let mut values = Vec::with_capacity(n_dims);
            for dim in 0..n_dims {
                let bin = permutations[dim][row];
                let jitter: f64 = rng.gen_range(0.0..1.0);
                let u = (bin as f64 + jitter) / n as f64;
                let r = &ranges[dim];
                let v = if r.log_scale {
                    let lo = r.low.max(1e-30).ln();
                    let hi = r.high.ln();
                    (lo + u * (hi - lo)).exp()
                } else {
                    r.low + u * (r.high - r.low)
                };
                values.push(v);
            }
            samples.push(SensitivitySample { names: names.clone(), values });
        }

        // Drop permutations to silence unused-mut warning if any. (No-op.)
        permutations.clear();

        Self {
            samples,
            ranges: ranges.to_vec(),
            seed,
        }
    }

    /// Number of samples drawn.
    pub fn len(&self) -> usize { self.samples.len() }

    /// True iff no samples.
    pub fn is_empty(&self) -> bool { self.samples.is_empty() }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dummy_ranges() -> Vec<ParameterRange> {
        vec![
            ParameterRange {
                name: "vmax".into(), low: 0.01, high: 0.1, log_scale: true,
                source: "test".into(),
            },
            ParameterRange {
                name: "km".into(), low: 0.1, high: 10.0, log_scale: false,
                source: "test".into(),
            },
        ]
    }

    #[test]
    fn lhs_size() {
        let s = SensitivitySweep::generate(&dummy_ranges(), 64, 42);
        assert_eq!(s.len(), 64);
        assert_eq!(s.samples[0].values.len(), 2);
    }

    #[test]
    fn lhs_in_range() {
        let s = SensitivitySweep::generate(&dummy_ranges(), 32, 1);
        for sample in &s.samples {
            // vmax in [0.01, 0.1]
            assert!(sample.values[0] >= 0.01 - 1e-9);
            assert!(sample.values[0] <= 0.1 + 1e-9);
            // km in [0.1, 10.0]
            assert!(sample.values[1] >= 0.1 - 1e-9);
            assert!(sample.values[1] <= 10.0 + 1e-9);
        }
    }

    #[test]
    fn lhs_reproducible() {
        let a = SensitivitySweep::generate(&dummy_ranges(), 16, 7);
        let b = SensitivitySweep::generate(&dummy_ranges(), 16, 7);
        for (sa, sb) in a.samples.iter().zip(b.samples.iter()) {
            for (va, vb) in sa.values.iter().zip(sb.values.iter()) {
                assert!((va - vb).abs() < 1e-12);
            }
        }
    }
}
