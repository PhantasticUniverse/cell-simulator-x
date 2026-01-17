//! Parameter structures with citation metadata.
//!
//! All biological parameters must include their source citation.

use serde::{Deserialize, Serialize};
use std::path::Path;

/// Top-level parameters container
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Parameters {
    /// Geometry parameters (Fung-Tong coefficients, etc.)
    pub geometry: GeometryParameters,
    /// Membrane mechanical parameters
    pub membrane: MembraneParameters,
}

impl Parameters {
    /// Load parameters from JSON files, or use defaults if files don't exist
    pub fn load_or_default() -> Self {
        let geometry = GeometryParameters::load_or_default("data/parameters/geometry.json");
        let membrane = MembraneParameters::load_or_default("data/parameters/membrane.json");

        Self { geometry, membrane }
    }

    /// Load parameters from specific directory
    pub fn load_from_dir<P: AsRef<Path>>(dir: P) -> Self {
        let dir = dir.as_ref();
        let geometry = GeometryParameters::load_or_default(dir.join("geometry.json"));
        let membrane = MembraneParameters::load_or_default(dir.join("membrane.json"));

        Self { geometry, membrane }
    }
}

impl Default for Parameters {
    fn default() -> Self {
        Self {
            geometry: GeometryParameters::default(),
            membrane: MembraneParameters::default(),
        }
    }
}

/// Geometry parameters for RBC shape
///
/// Based on Fung-Tong parametric equations.
/// Reference: Fung YC, Tong P, Patitucci P. Biorheology, 1981.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeometryParameters {
    /// Cell radius (μm)
    /// Reference: Mean value for normal human RBC
    /// Source: Evans & Fung, Microvasc Res 1972
    pub cell_radius_um: f32,

    /// Fung-Tong C₀ coefficient (μm)
    /// Determines center thickness
    /// Source: Fung et al., Biorheology 1981
    pub fung_tong_c0_um: f32,

    /// Fung-Tong C₂ coefficient (μm)
    /// Source: Fung et al., Biorheology 1981
    pub fung_tong_c2_um: f32,

    /// Fung-Tong C₄ coefficient (μm)
    /// Source: Fung et al., Biorheology 1981
    pub fung_tong_c4_um: f32,

    /// Mesh resolution (number of radial divisions)
    pub mesh_resolution: usize,

    /// Target number of spectrin tetramers
    /// Reference: ~33,000 per cell
    /// Source: Lux, Blood 2016
    pub spectrin_target_count: usize,
}

impl GeometryParameters {
    /// Load from JSON file or return defaults
    pub fn load_or_default<P: AsRef<Path>>(path: P) -> Self {
        match std::fs::read_to_string(path.as_ref()) {
            Ok(contents) => match serde_json::from_str(&contents) {
                Ok(params) => {
                    log::info!("Loaded geometry parameters from {:?}", path.as_ref());
                    params
                }
                Err(e) => {
                    log::warn!("Failed to parse geometry parameters: {}, using defaults", e);
                    Self::default()
                }
            },
            Err(_) => {
                log::info!("Geometry parameters file not found, using defaults");
                Self::default()
            }
        }
    }
}

impl Default for GeometryParameters {
    fn default() -> Self {
        Self {
            // Evans & Fung, Microvasc Res 1972
            cell_radius_um: 3.91,

            // Fung et al., Biorheology 1981
            fung_tong_c0_um: 0.81,
            fung_tong_c2_um: 7.83,
            fung_tong_c4_um: -4.39,

            // Simulation settings
            mesh_resolution: 50,
            spectrin_target_count: 33000,
        }
    }
}

/// Membrane mechanical parameters
///
/// Properties of the RBC membrane composite (lipid bilayer + spectrin network).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MembraneParameters {
    /// Shear modulus (μN/m)
    /// Reference: 5.5 ± 1.8 μN/m
    /// Source: Evans & Waugh, Biophys J 1977; Henon et al., Biophys J 1999
    pub shear_modulus_uN_per_m: f32,

    /// Bending modulus (pN·μm = 10⁻¹⁸ J)
    /// Reference: 1.8 × 10⁻¹⁹ J ≈ 44 kT at 37°C
    /// Source: Evans, Biophys J 1983
    pub bending_modulus_pN_um: f32,

    /// Area expansion modulus (mN/m)
    /// Reference: ~450 mN/m
    /// Source: Evans et al., Biophys J 1976
    pub area_modulus_mN_per_m: f32,

    /// Spectrin contour length (μm)
    /// Reference: ~200 nm fully extended
    /// Source: Rief et al., J Mol Biol 1999
    pub spectrin_contour_length_um: f32,

    /// Spectrin persistence length (μm)
    /// Reference: ~20 nm
    /// Source: Rief et al., J Mol Biol 1999
    pub spectrin_persistence_length_um: f32,

    /// Spectrin rest length (μm)
    /// Reference: ~75 nm in situ
    /// Source: Lux, Blood 2016
    pub spectrin_rest_length_um: f32,

    /// Membrane viscosity (pN·s/μm)
    /// Reference: Surface viscosity
    /// Source: Hochmuth et al., Biophys J 1979
    pub membrane_viscosity_pN_s_per_um: f32,
}

impl MembraneParameters {
    /// Load from JSON file or return defaults
    pub fn load_or_default<P: AsRef<Path>>(path: P) -> Self {
        match std::fs::read_to_string(path.as_ref()) {
            Ok(contents) => match serde_json::from_str(&contents) {
                Ok(params) => {
                    log::info!("Loaded membrane parameters from {:?}", path.as_ref());
                    params
                }
                Err(e) => {
                    log::warn!("Failed to parse membrane parameters: {}, using defaults", e);
                    Self::default()
                }
            },
            Err(_) => {
                log::info!("Membrane parameters file not found, using defaults");
                Self::default()
            }
        }
    }
}

impl Default for MembraneParameters {
    fn default() -> Self {
        Self {
            // Evans & Waugh 1977, Henon et al. 1999
            shear_modulus_uN_per_m: 5.5,

            // Evans 1983
            bending_modulus_pN_um: 0.18,

            // Evans et al. 1976
            area_modulus_mN_per_m: 450.0,

            // Rief et al. 1999
            spectrin_contour_length_um: 0.200,
            spectrin_persistence_length_um: 0.020,

            // Lux 2016
            spectrin_rest_length_um: 0.075,

            // Hochmuth et al. 1979
            membrane_viscosity_pN_s_per_um: 0.6,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_geometry_params() {
        let params = GeometryParameters::default();
        assert!((params.cell_radius_um - 3.91).abs() < 0.01);
    }

    #[test]
    fn test_default_membrane_params() {
        let params = MembraneParameters::default();
        assert!((params.shear_modulus_uN_per_m - 5.5).abs() < 0.1);
    }

    #[test]
    fn test_serialization() {
        let params = Parameters::default();
        let json = serde_json::to_string_pretty(&params).unwrap();
        let parsed: Parameters = serde_json::from_str(&json).unwrap();
        assert!((parsed.geometry.cell_radius_um - params.geometry.cell_radius_um).abs() < 0.001);
    }
}
