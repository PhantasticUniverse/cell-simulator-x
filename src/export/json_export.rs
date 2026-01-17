//! JSON state export for simulation snapshots.

use std::path::PathBuf;

use anyhow::Result;
use chrono::Local;
use serde::Serialize;

use crate::state::SimulationMetrics;

/// Full state export structure
#[derive(Debug, Clone, Serialize)]
pub struct StateExport {
    /// Export timestamp
    pub exported_at: String,
    /// Export version for compatibility
    pub version: &'static str,
    /// Simulation metrics snapshot
    pub metrics: SimulationMetrics,
}

/// Export current simulation state to JSON
///
/// Creates the exports directory if it doesn't exist.
/// Filename is auto-generated with timestamp: `state_YYYYMMDD_HHMMSS.json`
///
/// Returns the path to the saved JSON file.
pub fn export_state_json(metrics: &SimulationMetrics) -> Result<PathBuf> {
    // Create exports directory
    let dir = PathBuf::from("exports");
    std::fs::create_dir_all(&dir)?;

    // Generate filename with timestamp
    let timestamp = Local::now();
    let filename = format!("state_{}.json", timestamp.format("%Y%m%d_%H%M%S"));
    let path = dir.join(&filename);

    // Create export structure
    let export = StateExport {
        exported_at: timestamp.to_rfc3339(),
        version: "1.0.0",
        metrics: metrics.clone(),
    };

    // Write JSON
    let file = std::fs::File::create(&path)?;
    serde_json::to_writer_pretty(file, &export)?;

    log::info!("JSON state exported: {}", path.display());
    Ok(path)
}

/// Export state to a specific file
pub fn export_state_json_to(metrics: &SimulationMetrics, path: &PathBuf) -> Result<()> {
    let timestamp = Local::now();

    let export = StateExport {
        exported_at: timestamp.to_rfc3339(),
        version: "1.0.0",
        metrics: metrics.clone(),
    };

    let file = std::fs::File::create(path)?;
    serde_json::to_writer_pretty(file, &export)?;

    log::info!("JSON state exported: {}", path.display());
    Ok(())
}
