//! Spectrin-actin cytoskeleton network.
//!
//! The RBC membrane skeleton consists of:
//! - ~33,000 spectrin tetramers (αβ heterodimers that associate head-to-head)
//! - Junctional complexes at network nodes (short actin filaments + associated proteins)
//! - Network forms a quasi-hexagonal lattice underlying the lipid bilayer
//!
//! Reference: Mohandas & Gallagher, Blood 2008
//! Reference: Lux, Blood 2016

use bytemuck::{Pod, Zeroable};
use glam::Vec3;

use super::Mesh;
use crate::config::GeometryParameters;

/// A node in the spectrin network (junctional complex)
///
/// Each junction contains:
/// - Short actin protofilament (~13-14 monomers)
/// - Tropomyosin
/// - Protein 4.1
/// - Adducin
/// - Dematin
#[repr(C)]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
pub struct SpectrinNode {
    /// Position in μm
    pub position: [f32; 3],
    /// Index of nearest mesh vertex
    pub mesh_vertex_idx: u32,
    /// Number of connected spectrin tetramers (typically 5-7, usually 6)
    pub degree: u32,
    /// Padding for GPU alignment
    _padding: [u32; 3],
}

impl SpectrinNode {
    pub fn new(position: Vec3, mesh_vertex_idx: u32, degree: u32) -> Self {
        Self {
            position: position.to_array(),
            mesh_vertex_idx,
            degree,
            _padding: [0; 3],
        }
    }

    pub fn position_vec3(&self) -> Vec3 {
        Vec3::from_array(self.position)
    }
}

/// A spectrin tetramer connecting two junctional complexes
///
/// Properties:
/// - Contour length: ~200 nm (fully extended)
/// - Rest length: ~75 nm (in situ)
/// - Behaves as worm-like chain (WLC)
/// - Persistence length: ~20 nm
///
/// Reference: Rief et al., J Mol Biol 1999
#[repr(C)]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
pub struct SpectrinEdge {
    /// Index of first node
    pub node_a: u32,
    /// Index of second node
    pub node_b: u32,
    /// Rest length in μm (typically ~0.075 μm = 75 nm)
    pub rest_length_um: f32,
    /// Contour length in μm (typically ~0.2 μm = 200 nm)
    pub contour_length_um: f32,
}

impl SpectrinEdge {
    pub fn new(node_a: u32, node_b: u32) -> Self {
        Self {
            node_a,
            node_b,
            rest_length_um: 0.075,    // 75 nm
            contour_length_um: 0.200, // 200 nm
        }
    }
}

/// The complete spectrin-actin cytoskeleton network
pub struct SpectrinNetwork {
    /// Junctional complexes
    pub nodes: Vec<SpectrinNode>,
    /// Spectrin tetramers
    pub edges: Vec<SpectrinEdge>,
}

impl SpectrinNetwork {
    /// Generate spectrin network on the RBC surface
    ///
    /// Creates a quasi-hexagonal lattice of junctional complexes
    /// connected by spectrin tetramers.
    pub fn generate(mesh: &Mesh, params: &GeometryParameters) -> Self {
        // Estimate number of nodes needed for target spectrin count
        // For a hexagonal network: edges ≈ 3 * nodes
        // So for ~33,000 spectrin: ~11,000 nodes
        let target_spectrin = params.spectrin_target_count;
        let target_nodes = target_spectrin / 3;

        // Generate nodes by sampling mesh vertices with subdivision
        let (nodes, edges) = Self::generate_hexagonal_network(mesh, target_nodes);

        Self { nodes, edges }
    }

    /// Generate hexagonal network by placing nodes on mesh surface
    fn generate_hexagonal_network(mesh: &Mesh, target_nodes: usize) -> (Vec<SpectrinNode>, Vec<SpectrinEdge>) {
        let mut nodes = Vec::with_capacity(target_nodes);
        let mut edges = Vec::new();

        // Calculate node spacing based on surface area and target count
        let surface_area = mesh.calculate_surface_area();
        let area_per_node = surface_area / target_nodes as f32;
        let spacing = (area_per_node * 2.0 / 3.0_f32.sqrt()).sqrt(); // Hexagonal packing

        // Use a simpler approach: sample vertices uniformly and connect neighbors
        // This creates an approximate hexagonal lattice

        // First, collect unique positions from mesh at appropriate density
        let mut sampled_positions: Vec<(Vec3, u32)> = Vec::new();

        for (idx, vertex) in mesh.vertices.iter().enumerate() {
            let pos = vertex.position_vec3();

            // Check if this position is far enough from existing samples
            let too_close = sampled_positions.iter().any(|(p, _)| {
                pos.distance(*p) < spacing * 0.7
            });

            if !too_close {
                sampled_positions.push((pos, idx as u32));
            }

            if sampled_positions.len() >= target_nodes {
                break;
            }
        }

        // If we don't have enough nodes, add more by interpolation
        if sampled_positions.len() < target_nodes / 2 {
            // Subdivide: add midpoints of triangles
            for chunk in mesh.indices.chunks(3) {
                if sampled_positions.len() >= target_nodes {
                    break;
                }

                let v1 = mesh.vertices[chunk[0] as usize].position_vec3();
                let v2 = mesh.vertices[chunk[1] as usize].position_vec3();
                let v3 = mesh.vertices[chunk[2] as usize].position_vec3();
                let center = (v1 + v2 + v3) / 3.0;

                let too_close = sampled_positions.iter().any(|(p, _)| {
                    center.distance(*p) < spacing * 0.7
                });

                if !too_close {
                    sampled_positions.push((center, chunk[0])); // Use first vertex as reference
                }
            }
        }

        // Create nodes
        for (pos, mesh_idx) in &sampled_positions {
            nodes.push(SpectrinNode::new(*pos, *mesh_idx, 0));
        }

        // Create edges by connecting nearby nodes (within ~1.5x spacing)
        let connect_distance = spacing * 1.5;
        let mut node_degrees: Vec<u32> = vec![0; nodes.len()];

        for i in 0..nodes.len() {
            for j in (i + 1)..nodes.len() {
                let pos_i = nodes[i].position_vec3();
                let pos_j = nodes[j].position_vec3();
                let dist = pos_i.distance(pos_j);

                // Connect if within range and both nodes have < 7 connections
                if dist < connect_distance && node_degrees[i] < 7 && node_degrees[j] < 7 {
                    edges.push(SpectrinEdge::new(i as u32, j as u32));
                    node_degrees[i] += 1;
                    node_degrees[j] += 1;
                }
            }
        }

        // Update node degrees
        for (i, degree) in node_degrees.into_iter().enumerate() {
            nodes[i].degree = degree;
        }

        (nodes, edges)
    }

    /// Get the total length of all spectrin tetramers
    pub fn total_spectrin_length(&self) -> f32 {
        self.edges.iter().map(|e| {
            let pos_a = self.nodes[e.node_a as usize].position_vec3();
            let pos_b = self.nodes[e.node_b as usize].position_vec3();
            pos_a.distance(pos_b)
        }).sum()
    }

    /// Get average connectivity (degree) of nodes
    pub fn average_degree(&self) -> f32 {
        if self.nodes.is_empty() {
            return 0.0;
        }
        let total_degree: u32 = self.nodes.iter().map(|n| n.degree).sum();
        total_degree as f32 / self.nodes.len() as f32
    }

    /// Get edge positions as pairs of Vec3 (for line rendering)
    pub fn edge_positions(&self) -> Vec<(Vec3, Vec3)> {
        self.edges.iter().map(|e| {
            let pos_a = self.nodes[e.node_a as usize].position_vec3();
            let pos_b = self.nodes[e.node_b as usize].position_vec3();
            (pos_a, pos_b)
        }).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::Mesh;

    fn default_params() -> GeometryParameters {
        GeometryParameters {
            cell_radius_um: 3.91,
            fung_tong_c0_um: 0.81,
            fung_tong_c2_um: 7.83,
            fung_tong_c4_um: -4.39,
            mesh_resolution: 30,
            spectrin_target_count: 1000, // Reduced for testing
        }
    }

    #[test]
    fn test_network_generation() {
        let params = default_params();
        let mesh = Mesh::generate_rbc(&params);
        let network = SpectrinNetwork::generate(&mesh, &params);

        assert!(!network.nodes.is_empty());
        assert!(!network.edges.is_empty());
    }

    #[test]
    fn test_connectivity() {
        let params = default_params();
        let mesh = Mesh::generate_rbc(&params);
        let network = SpectrinNetwork::generate(&mesh, &params);

        // Average degree should be around 5-6 for hexagonal network
        let avg_degree = network.average_degree();
        assert!(avg_degree > 2.0 && avg_degree < 10.0,
                "Average degree {} out of expected range", avg_degree);
    }
}
