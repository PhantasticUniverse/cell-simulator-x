//! Phase 17.1 CLI runner for the empirical validation suite.
//!
//! Replaces `cargo run --features validation -- --validate` (which has been
//! removed in favour of this dedicated binary). Run with:
//!
//!     cargo run -p rbc-validation-suite --bin validate --release
//!
//! The behaviour is identical to the previous `run_validation()` in
//! `cell-simulator-x`'s `main.rs`: run the full suite, print a summary
//! table, persist a JSON report under `target/validation/<sha>.json`, and
//! exit non-zero iff any experiment failed.

use rbc_validation_suite::run_full_suite;
use std::process::ExitCode;

fn main() -> ExitCode {
    println!("=== Cell Simulator X — Phase 10 Validation ===\n");
    let report = run_full_suite();
    report.print_summary();

    let path = format!("target/validation/{}.json", report.commit_sha);
    if let Err(e) = report.write_json(&path) {
        eprintln!("warning: failed to write report to {path}: {e}");
    } else {
        println!("\nFull report written to {}", path);
    }

    if report.all_passed() {
        ExitCode::SUCCESS
    } else {
        eprintln!(
            "{} of {} experiments failed (see notes above)",
            report.total_count - report.passed_count,
            report.total_count,
        );
        ExitCode::FAILURE
    }
}
