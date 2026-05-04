//! Splenic interendothelial slit geometry (Phase C-hybrid.2).
//!
//! In the spleen, RBCs traverse 0.5–1.0 μm-wide endothelial slits as part
//! of normal clearance. This is the principal physiological filter that
//! separates deformable cells (which pass) from rigid or aged cells
//! (which are sequestered and removed).
//!
//! ## Geometry
//!
//! A slit is a rectangular slot:
//! - `width_um`     — narrow direction (0.5–1.0 μm physiologically)
//! - `height_um`    — wider transverse direction (~3 μm)
//! - `length_um`    — axial (flow direction) length (~5 μm)
//! - `axis`         — flow direction (typically +Z)
//! - `narrow_axis`  — narrow-direction axis (typically +X)
//! - `entry`        — center of the slit entrance face
//!
//! ## Velocity profile
//!
//! For a slit of width `w` ≪ height `h`, the steady-state Poiseuille
//! velocity is well approximated by a parabolic-in-narrow-direction
//! profile:
//!
//!   v_axial(x) = v_max · (1 - (2x/w)²)  for |x| < w/2
//!
//! where x is the perpendicular distance from the slit centerline in the
//! narrow direction. The wall shear rate is `4 v_max / w`, so even
//! modest centerline velocities produce enormous wall shears at slit
//! widths approaching 0.5 μm — the source of the splenic filter's
//! mechanical selectivity.
//!
//! ## References
//!
//! - Pivkin IV, Peng Z, Karniadakis GE, et al. PNAS 2016;113:7804-7809
//!   (splenic-slit RBC transit).
//! - Picas L, Rico F, Scheuring S. Biophys J 2013;103:51-60 (membrane
//!   bending stiffness vs splenic clearance).
//! - Drasler WJ, Smith CM, Keller KH. Biophys J 1987;52:357-365 (slit
//!   geometry and RBC residence times).

use glam::Vec3;

use crate::geometry::Mesh;
use crate::state::PhysicsState;

/// A splenic interendothelial slit. Models a finite rectangular slot of
/// length `length_um` along the flow axis, width `width_um` in the
/// narrow direction, and height `height_um` in the wide transverse
/// direction. The entry face is centered at `entry`; flow is along
/// `axis`.
#[derive(Debug, Clone, Copy)]
pub struct SplenicSlit {
    pub width_um: f32,
    pub height_um: f32,
    pub length_um: f32,
    pub axis: Vec3,
    pub narrow_axis: Vec3,
    pub entry: Vec3,
}

impl Default for SplenicSlit {
    fn default() -> Self {
        // Canonical splenic slit (Pivkin et al. 2016 typical values).
        Self {
            width_um: 0.7,
            height_um: 3.0,
            length_um: 5.0,
            axis: Vec3::Z,
            narrow_axis: Vec3::X,
            entry: Vec3::ZERO,
        }
    }
}

impl SplenicSlit {
    /// Build a slit with given width and the canonical height + length.
    pub fn with_width(width_um: f32) -> Self {
        Self { width_um, ..Self::default() }
    }

    /// True if the point lies inside the slit volume (between entry and
    /// exit planes, within width and height).
    pub fn contains(&self, p: Vec3) -> bool {
        let to_p = p - self.entry;
        let along = to_p.dot(self.axis);
        if along < 0.0 || along > self.length_um {
            return false;
        }
        let narrow = to_p.dot(self.narrow_axis);
        if narrow.abs() > self.width_um * 0.5 {
            return false;
        }
        // Wide direction: assumed perpendicular to both axis and narrow_axis.
        let wide_axis = self.axis.cross(self.narrow_axis).normalize_or_zero();
        let wide = to_p.dot(wide_axis);
        if wide.abs() > self.height_um * 0.5 {
            return false;
        }
        true
    }
}

/// Steady-state slit flow with parabolic-in-narrow-direction profile.
#[derive(Debug, Clone, Copy)]
pub struct SlitFlow {
    pub slit: SplenicSlit,
    pub max_velocity_um_per_sec: f32,
}

impl SlitFlow {
    pub fn new(slit: SplenicSlit, max_velocity_um_per_sec: f32) -> Self {
        Self { slit, max_velocity_um_per_sec }
    }

    /// Velocity vector at a position. Inside the slit, parabolic profile
    /// in the narrow direction along the flow axis. Outside the slit,
    /// zero (no-slip outside; ignores upstream/downstream chamber
    /// hydrodynamics).
    pub fn velocity_at(&self, p: Vec3) -> Vec3 {
        if !self.slit.contains(p) {
            return Vec3::ZERO;
        }
        let to_p = p - self.slit.entry;
        let narrow = to_p.dot(self.slit.narrow_axis);
        let xi = 2.0 * narrow / self.slit.width_um;
        let speed = self.max_velocity_um_per_sec * (1.0 - xi * xi).max(0.0);
        self.slit.axis * speed
    }

    /// Wall shear rate `4 v_max / w` (1/s) — peak shear at the slit
    /// walls, where deformability selection happens.
    pub fn wall_shear_rate_per_sec(&self) -> f32 {
        4.0 * self.max_velocity_um_per_sec / self.slit.width_um
    }

    /// Mean transit time for a passive tracer at the slit centerline:
    /// `length / v_max`.
    pub fn centerline_transit_time_sec(&self) -> f32 {
        self.slit.length_um / self.max_velocity_um_per_sec
    }
}

/// Stokes-style drag force on a vertex inside a SlitFlow.
///
/// `F_drag = drag_coeff * (v_fluid - v_vertex)`
///
/// Outside the slit, drag is zero (the cell is in the slack chamber).
pub fn slit_drag_force_uN(
    flow: &SlitFlow,
    vertex_pos_um: Vec3,
    vertex_vel_um_per_sec: Vec3,
    drag_coeff: f32,
) -> Vec3 {
    let v_fluid = flow.velocity_at(vertex_pos_um);
    drag_coeff * (v_fluid - vertex_vel_um_per_sec)
}

/// Apply slit-flow drag to a mesh's vertices, writing into
/// `state.external_forces_uN`. Vertices outside the slit get zero
/// external force.
pub fn apply_slit_drag_to_external_forces(
    state: &mut PhysicsState,
    mesh: &Mesh,
    flow: &SlitFlow,
    drag_coeff: f32,
) {
    let n = mesh.vertices.len();
    if state.external_forces_uN.len() != n {
        state.external_forces_uN = vec![Vec3::ZERO; n];
    }
    if state.vertex_velocities_um_per_sec.len() != n {
        state.vertex_velocities_um_per_sec = vec![Vec3::ZERO; n];
    }
    for i in 0..n {
        let pos = mesh.vertices[i].position_vec3();
        let vel = state.vertex_velocities_um_per_sec[i];
        state.external_forces_uN[i] = slit_drag_force_uN(flow, pos, vel, drag_coeff);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn default_flow() -> SlitFlow {
        SlitFlow::new(SplenicSlit::default(), 200.0)
    }

    #[test]
    fn velocity_zero_outside_slit() {
        let flow = default_flow();
        // Far from slit entrance.
        let p = Vec3::new(0.0, 0.0, 100.0);
        assert_eq!(flow.velocity_at(p), Vec3::ZERO);
    }

    #[test]
    fn velocity_max_at_center() {
        let flow = default_flow();
        // Centerline at midpoint of slit length.
        let p = Vec3::new(0.0, 0.0, 2.5);
        let v = flow.velocity_at(p);
        assert!((v.length() - flow.max_velocity_um_per_sec).abs() < 1e-3);
        assert!(v.z > 0.0);
    }

    #[test]
    fn velocity_zero_at_wall() {
        let flow = default_flow();
        // At narrow-direction wall (x = +width/2).
        let p = Vec3::new(flow.slit.width_um * 0.5, 0.0, 2.5);
        let v = flow.velocity_at(p);
        assert!(v.length() < 1e-3, "wall velocity: {:?}", v);
    }

    #[test]
    fn parabolic_profile_in_narrow_direction() {
        let flow = default_flow();
        // At x = w/4 (half-way to wall), speed = v_max·(1 - 0.5²) = 0.75 v_max.
        let p = Vec3::new(flow.slit.width_um * 0.25, 0.0, 2.5);
        let v = flow.velocity_at(p);
        let expected = flow.max_velocity_um_per_sec * 0.75;
        assert!((v.length() - expected).abs() < 1e-3, "v at x=w/4: {} (expected {})", v.length(), expected);
    }

    #[test]
    fn wall_shear_high_at_narrow_slit() {
        let narrow = SlitFlow::new(SplenicSlit::with_width(0.5), 100.0);
        let wide = SlitFlow::new(SplenicSlit::with_width(1.0), 100.0);
        // Narrower slit → higher wall shear.
        assert!(narrow.wall_shear_rate_per_sec() > wide.wall_shear_rate_per_sec());
        // For w = 0.5 μm, γ̇ = 4 · 100 / 0.5 = 800 /s.
        assert!((narrow.wall_shear_rate_per_sec() - 800.0).abs() < 1e-3);
    }

    #[test]
    fn transit_time_canonical() {
        let flow = default_flow();
        // length = 5 μm, v_max = 200 μm/s → transit = 0.025 s.
        assert!((flow.centerline_transit_time_sec() - 0.025).abs() < 1e-5);
    }

    #[test]
    fn contains_inside_and_outside() {
        let slit = SplenicSlit::default();
        // Inside.
        assert!(slit.contains(Vec3::new(0.0, 0.0, 2.5)));
        // Outside (past exit plane).
        assert!(!slit.contains(Vec3::new(0.0, 0.0, 10.0)));
        // Outside (narrow direction).
        assert!(!slit.contains(Vec3::new(slit.width_um, 0.0, 2.5)));
    }
}
