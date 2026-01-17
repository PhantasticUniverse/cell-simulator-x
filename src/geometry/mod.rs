//! Geometry module for RBC mesh and cytoskeleton representation.
//!
//! Contains the Fung-Tong parametric surface for biconcave disc shape
//! and spectrin-actin network topology.

mod fung_tong;
mod mesh;
mod spectrin;

pub use fung_tong::FungTong;
pub use mesh::{Mesh, Vertex};
pub use spectrin::{SpectrinEdge, SpectrinNetwork, SpectrinNode};
