//! Empirical validation against published experimental reference curves.
//!
//! This module is the Phase 10 forcing function: it compares the simulator's
//! outputs against quantitative published data (Imai 1981, Mulquiney 1999,
//! Rief 1999, Waugh & Evans 1979, Dao 2003) using χ²/dof as the primary metric.
//!
//! Phase 17.1: carved out into its own workspace member crate so the
//! validation harness is a separately citable artifact for the preprint.

// Allow non-snake-case for unit suffixes in field/variable names (μN, μm, pN, etc.).
// Inherited convention from `cell-simulator-x` ("Units in names": `force_uN`,
// `velocity_um_per_sec`, etc.).
#![allow(non_snake_case)]
//!
//! Unlike the existing `tests/` suite — which checks that outputs fall within
//! physiological ranges — this module evaluates how well the model reproduces
//! measured curves with their published uncertainty bars.
//!
//! ## Layout
//!
//! - [`reference_curve`] — `ValidationCurve` type with x, y, sigma; CSV loader
//! - [`metrics`] — R², RMSE, NRMSE, χ²/dof
//! - [`experiments`] — one submodule per published experiment
//! - [`fitting`] — sensitivity sweep, identifiability analysis
//!
//! ## Usage
//!
//! ```rust,ignore
//! use cell_simulator_x::validation::{ValidationSuite, run_full_suite};
//!
//! let report = run_full_suite();
//! report.print_summary();
//! report.write_json("target/validation/report.json")?;
//! ```
//!
//! ## Why χ²/dof not just R²
//!
//! Published reference curves carry heteroscedastic uncertainty (different
//! sigma per point). R² treats all points equally; χ²/dof correctly weights
//! by reported sigma. A χ²/dof of 1.0 indicates the model fits within the
//! quoted experimental error; > 2 indicates poor fit; < 0.5 indicates
//! over-fitting or over-quoted error bars.

pub mod experiments;
pub mod fitting;
pub mod metrics;
pub mod reference_curve;

pub use metrics::{FitMetrics, compute_metrics};
pub use reference_curve::{ValidationCurve, CurveMetadata};

use serde::{Deserialize, Serialize};
use std::path::Path;

/// Per-experiment fit result.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExperimentResult {
    /// Identifier (e.g. `"imai_1981_oec_standard"`).
    pub name: String,
    /// Citation (e.g. `"Imai 1981, Methods Enzymol 76:438"`).
    pub citation: String,
    /// What the experiment measures.
    pub description: String,
    /// Computed fit metrics.
    pub metrics: FitMetrics,
    /// Pass/fail against the χ²/dof threshold for this experiment.
    pub passed: bool,
    /// Threshold that was applied (default 2.0).
    pub chi2_dof_threshold: f64,
    /// Free-form notes (e.g. why a specific fit fails).
    pub notes: Vec<String>,
}

/// Aggregate validation report across all experiments.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ValidationReport {
    /// Git commit SHA at time of run (or `"uncommitted"`).
    pub commit_sha: String,
    /// ISO 8601 timestamp.
    pub timestamp: String,
    /// Per-experiment results.
    pub experiments: Vec<ExperimentResult>,
    /// Summary: how many experiments passed.
    pub passed_count: usize,
    /// Summary: how many experiments ran.
    pub total_count: usize,
}

impl ValidationReport {
    /// True iff every experiment in the suite passed its threshold.
    pub fn all_passed(&self) -> bool {
        self.passed_count == self.total_count
    }

    /// Format a one-line-per-experiment summary suitable for CLI output.
    pub fn print_summary(&self) {
        println!("=== Validation Report ===");
        println!("Commit:    {}", self.commit_sha);
        println!("Timestamp: {}", self.timestamp);
        println!("Passed:    {}/{}", self.passed_count, self.total_count);
        println!();
        println!(
            "{:<40} {:>10} {:>10} {:>10} {:>8}",
            "Experiment", "χ²/dof", "RMSE", "R²", "Status"
        );
        println!("{}", "-".repeat(82));
        for exp in &self.experiments {
            let status = if exp.passed { "PASS" } else { "FAIL" };
            println!(
                "{:<40} {:>10.3} {:>10.4} {:>10.4} {:>8}",
                truncate(&exp.name, 40),
                exp.metrics.chi2_per_dof,
                exp.metrics.rmse,
                exp.metrics.r_squared,
                status
            );
        }
        if !self.all_passed() {
            println!();
            println!("Failures:");
            for exp in &self.experiments {
                if !exp.passed {
                    println!("  • {} — χ²/dof = {:.3} > {:.1}",
                        exp.name, exp.metrics.chi2_per_dof, exp.chi2_dof_threshold);
                    for note in &exp.notes {
                        println!("    note: {}", note);
                    }
                }
            }
        }
    }

    /// Serialize to JSON file at `path`. Creates parent directories.
    pub fn write_json<P: AsRef<Path>>(&self, path: P) -> std::io::Result<()> {
        if let Some(parent) = path.as_ref().parent() {
            std::fs::create_dir_all(parent)?;
        }
        let json = serde_json::to_string_pretty(self)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        std::fs::write(path, json)
    }
}

fn truncate(s: &str, max: usize) -> String {
    if s.len() <= max { s.to_string() } else { format!("{}…", &s[..max - 1]) }
}

/// Run the full validation suite end-to-end and return an aggregate report.
///
/// Each experiment is run independently; one failure does not abort the suite.
pub fn run_full_suite() -> ValidationReport {
    let experiments = vec![
        experiments::imai_1981::run_standard_oec(),
        experiments::imai_1981::run_bohr_shift(),
        experiments::imai_1981::run_dpg_shift(),
        experiments::mulquiney_1999::run_steady_state_metabolites(),
        experiments::mulquiney_1999::run_ppp_flux_fraction(),
        experiments::rief_1999::run_spectrin_force_extension(),
        experiments::waugh_evans_1979::run_shear_modulus(),
        experiments::dao_2003::run_axial_extension_curve(),
        experiments::fischer_2007::run_tank_treading(),
        experiments::skalak_1973::run_parachute_shape(),
    ];

    let passed_count = experiments.iter().filter(|e| e.passed).count();
    let total_count = experiments.len();

    ValidationReport {
        commit_sha: detect_commit_sha(),
        timestamp: current_timestamp(),
        experiments,
        passed_count,
        total_count,
    }
}

fn detect_commit_sha() -> String {
    std::process::Command::new("git")
        .args(["rev-parse", "--short", "HEAD"])
        .output()
        .ok()
        .filter(|out| out.status.success())
        .map(|out| String::from_utf8_lossy(&out.stdout).trim().to_string())
        .unwrap_or_else(|| "uncommitted".to_string())
}

fn current_timestamp() -> String {
    chrono::Utc::now().to_rfc3339()
}
