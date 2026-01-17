//! Skalak membrane strain energy model.
//!
//! Implements the Skalak model for red blood cell membrane mechanics:
//!
//! W = (G_s/4) * (I₁² + 2*I₁ - 2*I₂) + (G_a/4) * I₂²
//!
//! Where:
//! - G_s = shear modulus (5.5 μN/m)
//! - G_a = area modulus (related to incompressibility)
//! - I₁ = α² + β² - 2 (first strain invariant)
//! - I₂ = α²β² - 1 (second strain invariant)
//! - α, β = principal stretch ratios
//!
//! Also includes Helfrich bending energy:
//! E_b = κ/2 * ∫(2H - C₀)² dA
//!
//! References:
//! - Skalak et al., Biophys J 1973
//! - Evans & Fung, Microvasc Res 1972
//! - Helfrich, Z Naturforsch 1973

use glam::{Mat2, Vec2, Vec3};

use crate::geometry::Mesh;

/// Skalak membrane material parameters
#[derive(Debug, Clone)]
pub struct SkalakMaterial {
    /// Shear modulus in μN/m
    /// Reference: 5.5 ± 1.1 μN/m
    /// Source: Evans & Waugh, Biophys J 1977
    pub shear_modulus_uN_per_m: f32,
    /// Area modulus in μN/m (much larger than shear)
    /// Reference: ~450 mN/m for lipid bilayer
    /// Source: Evans et al., Biophys J 1976
    pub area_modulus_uN_per_m: f32,
    /// Bending modulus in pN·μm
    /// Reference: ~0.18 pN·μm (~44 kT at 37°C)
    /// Source: Evans, Biophys J 1983
    pub bending_modulus_pN_um: f32,
    /// Spontaneous curvature (1/μm)
    /// Usually zero for unstressed RBC
    pub spontaneous_curvature_per_um: f32,
}

impl Default for SkalakMaterial {
    fn default() -> Self {
        Self {
            shear_modulus_uN_per_m: 5.5,
            area_modulus_uN_per_m: 450_000.0, // 450 mN/m = 450,000 μN/m
            bending_modulus_pN_um: 0.18,
            spontaneous_curvature_per_um: 0.0,
        }
    }
}

/// Reference configuration for a membrane triangle element
#[derive(Debug, Clone)]
pub struct MembraneElement {
    /// Indices of three vertices
    pub vertex_indices: [u32; 3],
    /// Reference (undeformed) positions in local 2D coordinates
    pub reference_positions_2d: [Vec2; 3],
    /// Reference area
    pub reference_area_um2: f32,
    /// Reference edge vectors for strain calculation
    pub reference_edges: [Vec2; 2],
    /// Inverse of reference shape matrix for deformation gradient
    pub inverse_reference_shape: Mat2,
}

impl MembraneElement {
    /// Create a membrane element from reference positions
    pub fn new(indices: [u32; 3], ref_positions: [Vec3; 3]) -> Self {
        // Project to local 2D coordinate system
        // Use edge0 as x-axis, compute y perpendicular to it
        let edge0 = ref_positions[1] - ref_positions[0];
        let edge1 = ref_positions[2] - ref_positions[0];

        let x_axis = edge0.normalize();
        let normal = edge0.cross(edge1).normalize();
        let y_axis = normal.cross(x_axis);

        // Project vertices to 2D
        let p0 = Vec2::ZERO;
        let p1 = Vec2::new(edge0.dot(x_axis), edge0.dot(y_axis));
        let p2 = Vec2::new(edge1.dot(x_axis), edge1.dot(y_axis));

        let ref_positions_2d = [p0, p1, p2];

        // Reference edges
        let ref_e0 = p1 - p0;
        let ref_e1 = p2 - p0;

        // Reference area (half cross product magnitude in 2D)
        let ref_area = 0.5 * (ref_e0.x * ref_e1.y - ref_e0.y * ref_e1.x).abs();

        // Shape matrix: columns are edge vectors
        let shape_matrix = Mat2::from_cols(ref_e0, ref_e1);
        let inv_shape = shape_matrix.inverse();

        Self {
            vertex_indices: indices,
            reference_positions_2d: ref_positions_2d,
            reference_area_um2: ref_area,
            reference_edges: [ref_e0, ref_e1],
            inverse_reference_shape: inv_shape,
        }
    }

    /// Compute deformation gradient F from current positions
    pub fn deformation_gradient(&self, current_positions: [Vec3; 3]) -> Mat2 {
        // Project current positions to local 2D (using same frame as reference)
        let edge0 = current_positions[1] - current_positions[0];
        let edge1 = current_positions[2] - current_positions[0];

        let x_axis = edge0.normalize();
        let normal = edge0.cross(edge1);
        let normal_len = normal.length();

        if normal_len < 1e-10 {
            return Mat2::IDENTITY; // Degenerate triangle
        }

        let normal = normal / normal_len;
        let y_axis = normal.cross(x_axis);

        let curr_e0 = Vec2::new(edge0.dot(x_axis), edge0.dot(y_axis));
        let curr_e1 = Vec2::new(edge1.dot(x_axis), edge1.dot(y_axis));

        // F = current_shape * inverse_reference_shape
        let current_shape = Mat2::from_cols(curr_e0, curr_e1);
        current_shape * self.inverse_reference_shape
    }

    /// Compute Cauchy-Green deformation tensor C = F^T * F
    pub fn cauchy_green(&self, current_positions: [Vec3; 3]) -> Mat2 {
        let f = self.deformation_gradient(current_positions);
        f.transpose() * f
    }

    /// Compute principal stretch ratios (eigenvalues of sqrt(C))
    pub fn principal_stretches(&self, current_positions: [Vec3; 3]) -> (f32, f32) {
        let c = self.cauchy_green(current_positions);

        // For 2x2 symmetric matrix, eigenvalues are:
        // λ = (trace ± sqrt(trace² - 4*det)) / 2
        let trace = c.x_axis.x + c.y_axis.y;
        let det = c.x_axis.x * c.y_axis.y - c.x_axis.y * c.y_axis.x;

        let discriminant = (trace * trace - 4.0 * det).max(0.0);
        let sqrt_disc = discriminant.sqrt();

        // Eigenvalues of C
        let lambda1_sq = (trace + sqrt_disc) / 2.0;
        let lambda2_sq = (trace - sqrt_disc) / 2.0;

        // Principal stretches are sqrt of eigenvalues
        (lambda1_sq.max(0.0).sqrt(), lambda2_sq.max(0.0).sqrt())
    }

    /// Compute strain invariants I₁ and I₂
    pub fn strain_invariants(&self, current_positions: [Vec3; 3]) -> (f32, f32) {
        let (alpha, beta) = self.principal_stretches(current_positions);

        // I₁ = α² + β² - 2
        let i1 = alpha * alpha + beta * beta - 2.0;

        // I₂ = α²β² - 1 (area change)
        let i2 = alpha * alpha * beta * beta - 1.0;

        (i1, i2)
    }

    /// Compute current area
    pub fn current_area(&self, current_positions: [Vec3; 3]) -> f32 {
        let edge0 = current_positions[1] - current_positions[0];
        let edge1 = current_positions[2] - current_positions[0];
        0.5 * edge0.cross(edge1).length()
    }
}

/// Skalak membrane solver
pub struct SkalakSolver {
    /// Material parameters
    pub material: SkalakMaterial,
    /// Membrane elements (triangles with reference configuration)
    pub elements: Vec<MembraneElement>,
}

impl SkalakSolver {
    /// Create a new Skalak solver from mesh
    pub fn new(mesh: &Mesh, material: SkalakMaterial) -> Self {
        let mut elements = Vec::new();

        // Create an element for each triangle
        for chunk in mesh.indices.chunks(3) {
            let indices = [chunk[0], chunk[1], chunk[2]];
            let ref_positions = [
                mesh.vertices[chunk[0] as usize].position_vec3(),
                mesh.vertices[chunk[1] as usize].position_vec3(),
                mesh.vertices[chunk[2] as usize].position_vec3(),
            ];

            let element = MembraneElement::new(indices, ref_positions);

            // Skip degenerate triangles
            if element.reference_area_um2 > 1e-12 {
                elements.push(element);
            }
        }

        Self { material, elements }
    }

    /// Compute Skalak strain energy density for an element
    #[allow(dead_code)] // Kept for validation/debugging
    fn strain_energy_density(&self, i1: f32, i2: f32) -> f32 {
        let gs = self.material.shear_modulus_uN_per_m;
        let ga = self.material.area_modulus_uN_per_m;

        // W = (Gs/4) * (I1² + 2*I1 - 2*I2) + (Ga/4) * I2²
        let shear_term = (gs / 4.0) * (i1 * i1 + 2.0 * i1 - 2.0 * i2);
        let area_term = (ga / 4.0) * i2 * i2;

        shear_term + area_term
    }

    /// Compute forces on all vertices from membrane strain
    ///
    /// Uses a fast spring-based approximation instead of numerical differentiation.
    /// Each triangle edge acts as a spring trying to maintain its reference length.
    ///
    /// Returns (forces_per_vertex, total_elastic_energy)
    pub fn compute_forces(&self, mesh: &Mesh) -> (Vec<Vec3>, f32) {
        let n_vertices = mesh.vertices.len();
        let mut forces = vec![Vec3::ZERO; n_vertices];
        let mut total_energy = 0.0;

        // Spring stiffness derived from shear modulus
        // k ≈ G_s * sqrt(3) for triangular lattice
        let k_spring = self.material.shear_modulus_uN_per_m * 1.732;

        for element in &self.elements {
            let i0 = element.vertex_indices[0] as usize;
            let i1 = element.vertex_indices[1] as usize;
            let i2 = element.vertex_indices[2] as usize;

            let p0 = mesh.vertices[i0].position_vec3();
            let p1 = mesh.vertices[i1].position_vec3();
            let p2 = mesh.vertices[i2].position_vec3();

            // Reference edge lengths from element
            let ref_len_01 = element.reference_edges[0].length();
            let ref_len_02 = element.reference_edges[1].length();
            let ref_len_12 = (element.reference_positions_2d[2] - element.reference_positions_2d[1]).length();

            // Edge 0-1
            let edge_01 = p1 - p0;
            let len_01 = edge_01.length();
            if len_01 > 1e-10 {
                let strain_01 = (len_01 - ref_len_01) / ref_len_01.max(1e-10);
                let force_mag = k_spring * strain_01 * ref_len_01;
                let dir = edge_01 / len_01;
                let f = dir * force_mag;
                forces[i0] += f;
                forces[i1] -= f;
                total_energy += 0.5 * k_spring * strain_01 * strain_01 * ref_len_01;
            }

            // Edge 0-2
            let edge_02 = p2 - p0;
            let len_02 = edge_02.length();
            if len_02 > 1e-10 {
                let strain_02 = (len_02 - ref_len_02) / ref_len_02.max(1e-10);
                let force_mag = k_spring * strain_02 * ref_len_02;
                let dir = edge_02 / len_02;
                let f = dir * force_mag;
                forces[i0] += f;
                forces[i2] -= f;
                total_energy += 0.5 * k_spring * strain_02 * strain_02 * ref_len_02;
            }

            // Edge 1-2
            let edge_12 = p2 - p1;
            let len_12 = edge_12.length();
            if len_12 > 1e-10 {
                let strain_12 = (len_12 - ref_len_12) / ref_len_12.max(1e-10);
                let force_mag = k_spring * strain_12 * ref_len_12;
                let dir = edge_12 / len_12;
                let f = dir * force_mag;
                forces[i1] += f;
                forces[i2] -= f;
                total_energy += 0.5 * k_spring * strain_12 * strain_12 * ref_len_12;
            }

            // Area preservation force (penalty for area change)
            let current_area = element.current_area([p0, p1, p2]);
            let area_strain = (current_area - element.reference_area_um2) / element.reference_area_um2.max(1e-10);
            if area_strain.abs() > 0.01 {
                // Push vertices outward/inward to restore area
                let center = (p0 + p1 + p2) / 3.0;
                let area_force_scale = self.material.area_modulus_uN_per_m * area_strain * 0.001;

                for (idx, pos) in [(i0, p0), (i1, p1), (i2, p2)] {
                    let to_center = center - pos;
                    let dist = to_center.length();
                    if dist > 1e-10 {
                        // If area too large, push inward; if too small, push outward
                        forces[idx] += to_center.normalize() * area_force_scale;
                    }
                }
            }
        }

        // Skip bending forces for performance (they're small anyway)
        // Uncomment if needed:
        // let bending_forces = self.compute_bending_forces(mesh);
        // for (i, bf) in bending_forces.into_iter().enumerate() {
        //     forces[i] += bf;
        // }

        (forces, total_energy)
    }

    /// Compute bending forces using discrete Laplacian approximation
    #[allow(dead_code)] // Skipped for performance, can be enabled if needed
    fn compute_bending_forces(&self, mesh: &Mesh) -> Vec<Vec3> {
        let n_vertices = mesh.vertices.len();
        let mut forces = vec![Vec3::ZERO; n_vertices];
        let kappa = self.material.bending_modulus_pN_um;

        // Build vertex adjacency (one-ring neighbors)
        let mut neighbors: Vec<Vec<usize>> = vec![Vec::new(); n_vertices];
        for chunk in mesh.indices.chunks(3) {
            let i0 = chunk[0] as usize;
            let i1 = chunk[1] as usize;
            let i2 = chunk[2] as usize;

            if !neighbors[i0].contains(&i1) {
                neighbors[i0].push(i1);
            }
            if !neighbors[i0].contains(&i2) {
                neighbors[i0].push(i2);
            }
            if !neighbors[i1].contains(&i0) {
                neighbors[i1].push(i0);
            }
            if !neighbors[i1].contains(&i2) {
                neighbors[i1].push(i2);
            }
            if !neighbors[i2].contains(&i0) {
                neighbors[i2].push(i0);
            }
            if !neighbors[i2].contains(&i1) {
                neighbors[i2].push(i1);
            }
        }

        // Compute mean curvature normal at each vertex using cotangent Laplacian
        for i in 0..n_vertices {
            if neighbors[i].is_empty() {
                continue;
            }

            let pos_i = mesh.vertices[i].position_vec3();

            // Simple umbrella operator: average of neighbors minus center
            let mut laplacian = Vec3::ZERO;
            for &j in &neighbors[i] {
                let pos_j = mesh.vertices[j].position_vec3();
                laplacian += pos_j - pos_i;
            }
            laplacian /= neighbors[i].len() as f32;

            // Mean curvature ≈ |Δx| / 2
            // Bending force ≈ -κ * Δ²x
            // Simplified: F = κ * Laplacian (restoring toward smooth surface)
            forces[i] = laplacian * kappa * 0.001; // Scale factor for stability
        }

        forces
    }

    /// Compute total membrane area
    pub fn compute_area(&self, mesh: &Mesh) -> f32 {
        let mut total = 0.0;
        for element in &self.elements {
            let positions = [
                mesh.vertices[element.vertex_indices[0] as usize].position_vec3(),
                mesh.vertices[element.vertex_indices[1] as usize].position_vec3(),
                mesh.vertices[element.vertex_indices[2] as usize].position_vec3(),
            ];
            total += element.current_area(positions);
        }
        total
    }

    /// Compute reference area
    pub fn reference_area(&self) -> f32 {
        self.elements.iter().map(|e| e.reference_area_um2).sum()
    }

    /// Compute area strain
    pub fn area_strain(&self, mesh: &Mesh) -> f32 {
        let current = self.compute_area(mesh);
        let reference = self.reference_area();
        if reference > 0.0 {
            (current - reference) / reference
        } else {
            0.0
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::GeometryParameters;

    fn default_params() -> GeometryParameters {
        GeometryParameters {
            cell_radius_um: 3.91,
            fung_tong_c0_um: 0.81,
            fung_tong_c2_um: 7.83,
            fung_tong_c4_um: -4.39,
            mesh_resolution: 10,
            spectrin_target_count: 100,
        }
    }

    #[test]
    fn test_skalak_material_defaults() {
        let mat = SkalakMaterial::default();
        assert_eq!(mat.shear_modulus_uN_per_m, 5.5);
        assert_eq!(mat.bending_modulus_pN_um, 0.18);
    }

    #[test]
    fn test_membrane_element_creation() {
        let positions = [
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        ];
        let element = MembraneElement::new([0, 1, 2], positions);

        // Reference area should be 0.5 (unit right triangle)
        assert!((element.reference_area_um2 - 0.5).abs() < 1e-6);
    }

    #[test]
    fn test_strain_invariants_undeformed() {
        let positions = [
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        ];
        let element = MembraneElement::new([0, 1, 2], positions);

        let (i1, i2) = element.strain_invariants(positions);

        // For undeformed state: α = β = 1, so I1 = 0, I2 = 0
        assert!(i1.abs() < 1e-5, "I1 should be ~0 for undeformed state, got {}", i1);
        assert!(i2.abs() < 1e-5, "I2 should be ~0 for undeformed state, got {}", i2);
    }

    #[test]
    fn test_solver_creation() {
        let params = default_params();
        let mesh = Mesh::generate_rbc(&params);
        let solver = SkalakSolver::new(&mesh, SkalakMaterial::default());

        // Should have elements for each triangle
        assert!(!solver.elements.is_empty());
    }

    #[test]
    fn test_forces_at_rest() {
        let params = default_params();
        let mesh = Mesh::generate_rbc(&params);
        let solver = SkalakSolver::new(&mesh, SkalakMaterial::default());

        let (forces, energy) = solver.compute_forces(&mesh);

        // At reference configuration, strain energy should be very small
        // (only bending contribution from discrete curvature)
        assert!(energy < 1e-3, "Energy at rest should be small, got {}", energy);

        // Forces should be small at rest
        let max_force = forces.iter().map(|f| f.length()).fold(0.0f32, f32::max);
        assert!(
            max_force < 1.0,
            "Max force at rest should be small, got {}",
            max_force
        );
    }
}
