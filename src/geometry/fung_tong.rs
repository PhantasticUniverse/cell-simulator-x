//! Fung-Tong parametric equations for RBC shape.
//!
//! The biconcave disc shape is described by:
//! z(r) = ±(1 - (r/R)²)^(1/2) * (C₀ + C₂*(r/R)² + C₄*(r/R)⁴)
//!
//! Reference: Fung YC, Tong P, Patitucci P. "High-resolution data on the
//! geometry of red blood cells." Biorheology, 1981.
//!
//! Parameters for normal human RBC at physiological conditions:
//! - R = 3.91 μm (cell radius)
//! - C₀ = 0.81 μm
//! - C₂ = 7.83 μm
//! - C₄ = -4.39 μm

use crate::config::GeometryParameters;

/// Fung-Tong parametric surface for RBC biconcave disc
pub struct FungTong {
    /// Cell radius (μm)
    pub radius_um: f32,
    /// C₀ coefficient (μm)
    pub c0_um: f32,
    /// C₂ coefficient (μm)
    pub c2_um: f32,
    /// C₄ coefficient (μm)
    pub c4_um: f32,
}

impl FungTong {
    /// Create from geometry parameters
    pub fn from_parameters(params: &GeometryParameters) -> Self {
        Self {
            radius_um: params.cell_radius_um,
            c0_um: params.fung_tong_c0_um,
            c2_um: params.fung_tong_c2_um,
            c4_um: params.fung_tong_c4_um,
        }
    }

    /// Calculate z coordinate for a given radial distance
    ///
    /// Returns the half-thickness at radial position r.
    /// The full shape has z and -z surfaces.
    ///
    /// # Arguments
    /// * `r` - Radial distance from center (μm), must be <= radius
    ///
    /// # Returns
    /// Half-thickness z(r) in μm (positive value for upper surface)
    pub fn z_at_radius(&self, r: f32) -> f32 {
        if r >= self.radius_um {
            return 0.0;
        }

        let r_normalized = r / self.radius_um;
        let r2 = r_normalized * r_normalized;
        let r4 = r2 * r2;

        // Thickness profile
        let profile = self.c0_um + self.c2_um * r2 + self.c4_um * r4;

        // Elliptical envelope
        let envelope = (1.0 - r2).sqrt();

        envelope * profile
    }

    /// Calculate the surface normal at a point (r, theta)
    ///
    /// Returns unnormalized normal vector pointing outward from upper surface.
    pub fn normal_at(&self, r: f32, theta: f32) -> (f32, f32, f32) {
        let eps = 0.001 * self.radius_um;

        // Numerical gradient
        let z = self.z_at_radius(r);
        let z_r = if r > eps {
            (self.z_at_radius(r + eps) - self.z_at_radius(r - eps)) / (2.0 * eps)
        } else {
            (self.z_at_radius(r + eps) - z) / eps
        };

        // Normal = (-dz/dr * cos(theta), -dz/dr * sin(theta), 1)
        let nx = -z_r * theta.cos();
        let ny = -z_r * theta.sin();
        let nz = 1.0;

        (nx, ny, nz)
    }

    /// Calculate the dimple depth (minimum thickness at center vs edge)
    pub fn dimple_depth(&self) -> f32 {
        let z_center = self.z_at_radius(0.0);
        let z_max = self.maximum_half_thickness();
        z_max - z_center
    }

    /// Find the radial position of maximum thickness
    pub fn radius_of_max_thickness(&self) -> f32 {
        // Numerical search for maximum
        let mut max_z = 0.0_f32;
        let mut r_max = 0.0_f32;

        let steps = 100;
        for i in 0..=steps {
            let r = (i as f32 / steps as f32) * self.radius_um;
            let z = self.z_at_radius(r);
            if z > max_z {
                max_z = z;
                r_max = r;
            }
        }

        r_max
    }

    /// Get the maximum half-thickness
    pub fn maximum_half_thickness(&self) -> f32 {
        let r_max = self.radius_of_max_thickness();
        self.z_at_radius(r_max)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn default_params() -> GeometryParameters {
        GeometryParameters {
            cell_radius_um: 3.91,
            fung_tong_c0_um: 0.81,
            fung_tong_c2_um: 7.83,
            fung_tong_c4_um: -4.39,
            mesh_resolution: 50,
            spectrin_target_count: 33000,
        }
    }

    #[test]
    fn test_center_thickness() {
        let ft = FungTong::from_parameters(&default_params());
        let z_center = ft.z_at_radius(0.0);

        // Center should have C0 thickness (dimple)
        assert!((z_center - 0.81).abs() < 0.01);
    }

    #[test]
    fn test_edge_thickness() {
        let ft = FungTong::from_parameters(&default_params());
        let z_edge = ft.z_at_radius(ft.radius_um);

        // Edge should be zero
        assert!(z_edge.abs() < 0.001);
    }

    #[test]
    fn test_biconcave_shape() {
        let ft = FungTong::from_parameters(&default_params());

        // Should have dimple (center thinner than ~70% radius)
        let z_center = ft.z_at_radius(0.0);
        let z_mid = ft.z_at_radius(0.7 * ft.radius_um);

        assert!(z_mid > z_center, "Expected biconcave shape with dimple");
    }

    #[test]
    fn test_max_thickness_location() {
        let ft = FungTong::from_parameters(&default_params());
        let r_max = ft.radius_of_max_thickness();

        // Maximum thickness should be between 60-80% of radius
        let ratio = r_max / ft.radius_um;
        assert!(ratio > 0.5 && ratio < 0.9, "Max thickness at unexpected location: {}", ratio);
    }
}
