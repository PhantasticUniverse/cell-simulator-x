//! CSV time-series export for simulation metrics.

use std::fs::File;
use std::path::PathBuf;

use anyhow::Result;
use chrono::Local;
use serde::Serialize;

use crate::state::SimulationMetrics;

/// Record for CSV time-series export
#[derive(Debug, Clone, Serialize)]
pub struct TimeSeriesRecord {
    /// Simulation time (seconds)
    pub time_sec: f64,
    /// ATP concentration (mM)
    pub atp_mM: f64,
    /// 2,3-DPG concentration (mM)
    pub dpg_2_3_mM: f64,
    /// Glucose concentration (mM)
    pub glucose_mM: f64,
    /// Lactate concentration (mM)
    pub lactate_mM: f64,
    /// O2 saturation (fraction)
    pub o2_saturation: f64,
    /// pH
    pub ph: f64,
    /// NADPH/NADP+ ratio
    pub nadph_nadp_ratio: f64,
    /// GSH/GSSG ratio
    pub gsh_gssg_ratio: f64,
    /// H2O2 (uM)
    pub h2o2_uM: f64,
    /// Na+ (mM)
    pub na_plus_mM: f64,
    /// K+ (mM)
    pub k_plus_mM: f64,
    /// Ca2+ (nM)
    pub ca_nM: f64,
    /// Membrane tension (pN/nm)
    pub membrane_tension: f64,
}

impl From<&SimulationMetrics> for TimeSeriesRecord {
    fn from(m: &SimulationMetrics) -> Self {
        Self {
            time_sec: m.simulation_time_sec,
            atp_mM: m.atp_mM,
            dpg_2_3_mM: m.dpg_2_3_mM,
            glucose_mM: m.glucose_mM,
            lactate_mM: m.lactate_mM,
            o2_saturation: m.o2_saturation,
            ph: m.ph,
            nadph_nadp_ratio: m.nadph_nadp_ratio,
            gsh_gssg_ratio: m.gsh_gssg_ratio,
            h2o2_uM: m.h2o2_uM,
            na_plus_mM: m.na_plus_mM,
            k_plus_mM: m.k_plus_mM,
            ca_nM: m.ca_nM,
            membrane_tension: m.membrane_tension_pN_per_nm,
        }
    }
}

/// CSV exporter for time-series data
pub struct CsvExporter {
    writer: csv::Writer<File>,
    /// Sample interval in seconds
    sample_interval_sec: f64,
    /// Last sample time
    last_sample_time: f64,
    /// Path to output file
    path: PathBuf,
}

impl CsvExporter {
    /// Create a new CSV exporter with the given sample interval
    ///
    /// Creates the exports directory if it doesn't exist.
    /// Filename is auto-generated with timestamp.
    pub fn new(sample_interval_sec: f64) -> Result<Self> {
        // Create exports directory
        let dir = PathBuf::from("exports");
        std::fs::create_dir_all(&dir)?;

        // Generate filename with timestamp
        let timestamp = Local::now().format("%Y%m%d_%H%M%S");
        let filename = format!("timeseries_{}.csv", timestamp);
        let path = dir.join(&filename);

        // Create writer
        let file = File::create(&path)?;
        let writer = csv::Writer::from_writer(file);

        log::info!("CSV export started: {}", path.display());

        Ok(Self {
            writer,
            sample_interval_sec,
            last_sample_time: -sample_interval_sec, // Ensure first sample is recorded
            path,
        })
    }

    /// Record a sample if the interval has elapsed
    pub fn maybe_record(&mut self, metrics: &SimulationMetrics) -> Result<bool> {
        let time = metrics.simulation_time_sec;

        if time - self.last_sample_time >= self.sample_interval_sec {
            let record = TimeSeriesRecord::from(metrics);
            self.writer.serialize(&record)?;
            self.last_sample_time = time;
            Ok(true)
        } else {
            Ok(false)
        }
    }

    /// Force record a sample regardless of interval
    pub fn record(&mut self, metrics: &SimulationMetrics) -> Result<()> {
        let record = TimeSeriesRecord::from(metrics);
        self.writer.serialize(&record)?;
        self.last_sample_time = metrics.simulation_time_sec;
        Ok(())
    }

    /// Finish writing and return the output path
    pub fn finish(mut self) -> Result<PathBuf> {
        self.writer.flush()?;
        log::info!("CSV export completed: {}", self.path.display());
        Ok(self.path)
    }

    /// Get the output path
    pub fn path(&self) -> &PathBuf {
        &self.path
    }
}
