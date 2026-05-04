//! Reference curve type with measurement uncertainty.
//!
//! Each `ValidationCurve` represents a digitized published figure or table:
//! independent variable `x`, measured response `y`, and the published
//! 1-sigma uncertainty per point. CSV files in `data/reference_curves/`
//! deserialize into this type.

use serde::{Deserialize, Serialize};
use std::path::Path;

/// Metadata for a reference curve, used for provenance and citation tracking.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CurveMetadata {
    /// Short identifier, e.g. `"imai_1981_oec_standard"`.
    pub id: String,
    /// Full citation including author, year, journal.
    pub citation: String,
    /// DOI of the source paper (or `""` if unavailable).
    pub doi: String,
    /// Figure/table number (e.g. `"Fig. 3A"`, `"Table 2"`).
    pub figure: String,
    /// Method used to extract data ("WebPlotDigitizer", "tabulated values", "regenerated from published Adair fit").
    pub digitization: String,
    /// Independent-variable label and unit, e.g. `"pO2 [mmHg]"`.
    pub x_label: String,
    /// Dependent-variable label and unit, e.g. `"saturation [fraction]"`.
    pub y_label: String,
    /// Free-form notes (experimental conditions, caveats).
    pub notes: String,
}

/// A reference curve: x, y, sigma with metadata.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ValidationCurve {
    pub metadata: CurveMetadata,
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    /// Per-point 1-sigma uncertainty in `y`. Same length as `x`/`y`.
    pub sigma: Vec<f64>,
}

impl ValidationCurve {
    /// Construct from raw vectors. Panics if lengths disagree.
    pub fn new(metadata: CurveMetadata, x: Vec<f64>, y: Vec<f64>, sigma: Vec<f64>) -> Self {
        assert_eq!(x.len(), y.len(), "ValidationCurve: x.len() != y.len()");
        assert_eq!(x.len(), sigma.len(), "ValidationCurve: x.len() != sigma.len()");
        assert!(!x.is_empty(), "ValidationCurve: must have at least one point");
        Self { metadata, x, y, sigma }
    }

    /// Number of data points.
    pub fn len(&self) -> usize { self.x.len() }

    /// True iff there are no points.
    pub fn is_empty(&self) -> bool { self.x.is_empty() }

    /// Load from a CSV file. The CSV must have columns x, y, sigma in any
    /// order, plus optional metadata in a leading `# key=value` block.
    /// Provenance metadata is supplied separately because CSV is fragile.
    pub fn from_csv<P: AsRef<Path>>(
        path: P,
        metadata: CurveMetadata,
    ) -> Result<Self, csv::Error> {
        let mut rdr = csv::ReaderBuilder::new()
            .comment(Some(b'#'))
            .from_path(path)?;

        let headers = rdr.headers()?.clone();
        let x_idx = column_index(&headers, "x").ok_or_else(|| {
            csv::Error::from(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "missing 'x' column",
            ))
        })?;
        let y_idx = column_index(&headers, "y").ok_or_else(|| {
            csv::Error::from(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "missing 'y' column",
            ))
        })?;
        let sigma_idx = column_index(&headers, "sigma").ok_or_else(|| {
            csv::Error::from(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "missing 'sigma' column",
            ))
        })?;

        let mut x = Vec::new();
        let mut y = Vec::new();
        let mut sigma = Vec::new();
        for record in rdr.records() {
            let record = record?;
            x.push(parse_field(&record, x_idx)?);
            y.push(parse_field(&record, y_idx)?);
            sigma.push(parse_field(&record, sigma_idx)?);
        }

        Ok(Self::new(metadata, x, y, sigma))
    }
}

fn column_index(headers: &csv::StringRecord, name: &str) -> Option<usize> {
    headers.iter().position(|h| h.trim().eq_ignore_ascii_case(name))
}

fn parse_field(record: &csv::StringRecord, idx: usize) -> Result<f64, csv::Error> {
    let s = record.get(idx).ok_or_else(|| {
        csv::Error::from(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "missing field",
        ))
    })?;
    s.trim().parse::<f64>().map_err(|e| {
        csv::Error::from(std::io::Error::new(std::io::ErrorKind::InvalidData, e))
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dummy_metadata() -> CurveMetadata {
        CurveMetadata {
            id: "test".into(),
            citation: "Test 2026".into(),
            doi: "".into(),
            figure: "Fig. 1".into(),
            digitization: "synthetic".into(),
            x_label: "x".into(),
            y_label: "y".into(),
            notes: "".into(),
        }
    }

    #[test]
    fn construction() {
        let curve = ValidationCurve::new(
            dummy_metadata(),
            vec![1.0, 2.0, 3.0],
            vec![10.0, 20.0, 30.0],
            vec![0.5, 0.5, 0.5],
        );
        assert_eq!(curve.len(), 3);
    }

    #[test]
    #[should_panic]
    fn rejects_mismatched_lengths() {
        ValidationCurve::new(
            dummy_metadata(),
            vec![1.0, 2.0],
            vec![10.0],
            vec![0.5, 0.5],
        );
    }
}
