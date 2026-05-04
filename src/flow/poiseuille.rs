//! Analytic Poiseuille flow in a cylindrical channel.
//!
//! For a cylindrical pipe of radius `R` aligned along the +Z axis, the
//! steady-state laminar (Poiseuille) velocity profile is
//!
//!   v_z(r) = v_max * (1 - (r/R)²),    v_x = v_y = 0
//!
//! where `r = sqrt(x² + y²)` and `v_max` is the centerline speed. The
//! mean velocity is `v_max / 2` and the wall shear rate is `2 v_max / R`.
//!
//! Per-vertex drag is approximated as Stokes drag on a small effective
//! sphere:
//!
//!   F_drag = drag_coeff * (v_fluid - v_vertex)
//!
//! The drag coefficient absorbs the viscosity, hydrodynamic radius, and
//! local geometry — for an RBC vertex with effective radius ~0.05 µm in
//! plasma (μ ≈ 1.2 mPa·s), `6πμa ≈ 1.1 × 10⁻⁹ N·s/m = 1.1 µN·s/µm`. The
//! default below uses a smaller value because vertex velocities in this
//! simulator are in µm/s while forces are in µN; calibration is left to
//! validation in Phase 12.C.

use glam::Vec3;

use crate::geometry::Mesh;
use crate::state::PhysicsState;

/// A right-circular cylindrical channel aligned along the axis. The
/// vessel inscribes the radius `R` and is infinite along the axis — the
/// Poiseuille profile applies wherever a query point sits (no
/// inlet/outlet boundary conditions yet).
#[derive(Debug, Clone, Copy)]
pub struct CylindricalChannel {
    /// Radius of the channel (μm).
    pub radius_um: f32,
    /// Axis direction, normalized.
    pub axis: Vec3,
    /// A point on the centerline (μm).
    pub center: Vec3,
}

impl Default for CylindricalChannel {
    fn default() -> Self {
        // 5 µm radius — comparable to a small arteriole.
        Self {
            radius_um: 5.0,
            axis: Vec3::Z,
            center: Vec3::ZERO,
        }
    }
}

impl CylindricalChannel {
    /// Radial distance from a point to the centerline.
    pub fn radial_distance_um(&self, p: Vec3) -> f32 {
        let to_p = p - self.center;
        let along = to_p.dot(self.axis);
        let radial = to_p - self.axis * along;
        radial.length()
    }
}

/// Steady-state Poiseuille flow inside a [`CylindricalChannel`].
#[derive(Debug, Clone, Copy)]
pub struct Poiseuille {
    /// Channel geometry.
    pub channel: CylindricalChannel,
    /// Centerline (peak) velocity in μm/s.
    pub max_velocity_um_per_sec: f32,
}

impl Poiseuille {
    /// Construct a Poiseuille flow with given channel and centerline speed.
    pub fn new(channel: CylindricalChannel, max_velocity_um_per_sec: f32) -> Self {
        Self {
            channel,
            max_velocity_um_per_sec,
        }
    }

    /// Mean velocity (`v_max / 2`).
    pub fn mean_velocity_um_per_sec(&self) -> f32 {
        self.max_velocity_um_per_sec * 0.5
    }

    /// Wall shear rate: `dv/dr|_{r=R}` = `2 v_max / R` (1/s).
    pub fn wall_shear_rate_per_sec(&self) -> f32 {
        2.0 * self.max_velocity_um_per_sec / self.channel.radius_um
    }

    /// Velocity vector at a position (μm/s). Inside the channel, parabolic
    /// profile along the axis. Outside the channel, zero (no-slip outside).
    pub fn velocity_at(&self, p: Vec3) -> Vec3 {
        let r = self.channel.radial_distance_um(p);
        if r >= self.channel.radius_um {
            return Vec3::ZERO;
        }
        let xi = r / self.channel.radius_um;
        let speed = self.max_velocity_um_per_sec * (1.0 - xi * xi);
        self.channel.axis * speed
    }

    /// Shear rate at radial position `r`: `dv/dr = -2 v_max r / R²` (1/s,
    /// signed; magnitude grows linearly toward the wall).
    pub fn shear_rate_at(&self, p: Vec3) -> f32 {
        let r = self.channel.radial_distance_um(p);
        if r >= self.channel.radius_um {
            return 0.0;
        }
        -2.0 * self.max_velocity_um_per_sec * r
            / (self.channel.radius_um * self.channel.radius_um)
    }
}

/// Stokes-style drag force on a vertex, given its position and velocity.
///
/// `F_drag = drag_coeff * (v_fluid - v_vertex)`
///
/// Returns force in μN.
pub fn drag_force_uN(
    flow: &Poiseuille,
    vertex_pos_um: Vec3,
    vertex_vel_um_per_sec: Vec3,
    drag_coeff: f32,
) -> Vec3 {
    let v_fluid = flow.velocity_at(vertex_pos_um);
    drag_coeff * (v_fluid - vertex_vel_um_per_sec)
}

/// Compute drag forces on every vertex of `mesh` and write them into
/// `state.external_forces_uN`. Existing external forces are overwritten.
///
/// Used as a pre-step before `PhysicsSolver::step` (CPU path) so the
/// integrator picks up the drag as part of its force accumulation. For
/// the GPU path (Phase 12.B), this gets replaced by a kernel.
pub fn apply_drag_to_external_forces(
    state: &mut PhysicsState,
    mesh: &Mesh,
    flow: &Poiseuille,
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
        state.external_forces_uN[i] = drag_force_uN(flow, pos, vel, drag_coeff);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn default_flow() -> Poiseuille {
        Poiseuille::new(CylindricalChannel::default(), 100.0)
    }

    #[test]
    fn velocity_zero_outside_channel() {
        let flow = default_flow();
        let outside = Vec3::new(10.0, 0.0, 0.0); // r = 10 > R = 5
        assert_eq!(flow.velocity_at(outside), Vec3::ZERO);
    }

    #[test]
    fn velocity_max_at_center() {
        let flow = default_flow();
        let v = flow.velocity_at(Vec3::ZERO);
        assert!((v.length() - flow.max_velocity_um_per_sec).abs() < 1e-5);
        assert!(v.dot(flow.channel.axis) > 0.0); // points along axis
    }

    #[test]
    fn velocity_zero_at_wall() {
        let flow = default_flow();
        let p = Vec3::new(flow.channel.radius_um, 0.0, 5.0);
        let v = flow.velocity_at(p);
        // Right at the wall, parabolic gives 0.
        assert!(v.length() < 1e-5);
    }

    #[test]
    fn parabolic_profile() {
        let flow = default_flow();
        // At r = R/2, v = v_max * 3/4.
        let p = Vec3::new(flow.channel.radius_um * 0.5, 0.0, 0.0);
        let v = flow.velocity_at(p);
        let expected = flow.max_velocity_um_per_sec * 0.75;
        assert!((v.length() - expected).abs() < 1e-4);
    }

    #[test]
    fn mean_and_wall_shear() {
        let flow = default_flow();
        assert!((flow.mean_velocity_um_per_sec() - 50.0).abs() < 1e-5);
        // Wall shear = 2 * 100 / 5 = 40 /s
        assert!((flow.wall_shear_rate_per_sec() - 40.0).abs() < 1e-5);
    }

    #[test]
    fn drag_force_zero_when_vertex_matches_fluid() {
        let flow = default_flow();
        let pos = Vec3::ZERO;
        let v_fluid = flow.velocity_at(pos);
        let f = drag_force_uN(&flow, pos, v_fluid, 0.5);
        assert!(f.length() < 1e-5);
    }

    #[test]
    fn drag_force_pulls_stationary_vertex_along_flow() {
        let flow = default_flow();
        let pos = Vec3::ZERO;
        let f = drag_force_uN(&flow, pos, Vec3::ZERO, 0.5);
        // At center, fluid moves along +Z, so drag on stationary vertex
        // should also push +Z.
        assert!(f.z > 0.0);
        assert!(f.x.abs() < 1e-5);
        assert!(f.y.abs() < 1e-5);
    }

    #[test]
    fn external_forces_populated_for_all_vertices() {
        use crate::config::GeometryParameters;
        let params = GeometryParameters {
            cell_radius_um: 3.91,
            fung_tong_c0_um: 0.81,
            fung_tong_c2_um: 7.83,
            fung_tong_c4_um: -4.39,
            mesh_resolution: 6,
            spectrin_target_count: 30,
        };
        let mesh = Mesh::generate_rbc(&params);
        let mut state = PhysicsState::new(mesh.vertices.len());
        let flow = default_flow();
        apply_drag_to_external_forces(&mut state, &mesh, &flow, 0.5);
        assert_eq!(state.external_forces_uN.len(), mesh.vertices.len());
        // Some forces should be non-zero (cells of radius ~4 µm sit inside
        // the 5 µm channel).
        let any_nonzero = state.external_forces_uN.iter().any(|f| f.length() > 1e-6);
        assert!(any_nonzero, "expected some non-zero drag forces");
    }
}
