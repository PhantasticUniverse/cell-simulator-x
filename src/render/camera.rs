//! Orbit camera for 3D visualization.
//!
//! Provides a camera that orbits around the RBC, allowing rotation
//! with mouse drag and zoom with scroll.

use bytemuck::{Pod, Zeroable};
use glam::{Mat4, Vec3};

/// Camera uniform data sent to GPU
#[repr(C)]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
pub struct CameraUniform {
    /// Combined view-projection matrix
    pub view_proj: [[f32; 4]; 4],
    /// Camera position for lighting calculations
    pub camera_pos: [f32; 4],
}

/// Orbit camera that rotates around a target point
pub struct Camera {
    /// Target point to orbit around (μm)
    pub target: Vec3,
    /// Distance from target (μm)
    pub distance: f32,
    /// Azimuthal angle (radians, rotation around Y axis)
    pub azimuth: f32,
    /// Polar angle (radians, elevation from horizontal)
    pub elevation: f32,
    /// Field of view (radians)
    pub fov: f32,
    /// Near clipping plane (μm)
    pub near: f32,
    /// Far clipping plane (μm)
    pub far: f32,
    /// Aspect ratio (width / height)
    pub aspect: f32,
    /// Auto-rotation speed (radians per second), 0 to disable
    pub auto_rotate_speed: f32,
}

impl Camera {
    /// Create a new camera with default settings for RBC viewing
    pub fn new(aspect: f32) -> Self {
        Self {
            target: Vec3::ZERO,
            distance: 15.0, // μm, gives good view of ~8 μm diameter cell
            azimuth: 0.0,
            elevation: 0.3, // Slight top-down view
            fov: std::f32::consts::FRAC_PI_4, // 45 degrees
            near: 0.1,
            far: 100.0,
            aspect,
            auto_rotate_speed: 0.3, // Slow rotation
        }
    }

    /// Calculate camera position from orbital parameters
    pub fn position(&self) -> Vec3 {
        let x = self.distance * self.elevation.cos() * self.azimuth.sin();
        let y = self.distance * self.elevation.sin();
        let z = self.distance * self.elevation.cos() * self.azimuth.cos();
        self.target + Vec3::new(x, y, z)
    }

    /// Calculate view matrix
    pub fn view_matrix(&self) -> Mat4 {
        Mat4::look_at_rh(self.position(), self.target, Vec3::Y)
    }

    /// Calculate projection matrix
    pub fn projection_matrix(&self) -> Mat4 {
        Mat4::perspective_rh(self.fov, self.aspect, self.near, self.far)
    }

    /// Get combined view-projection matrix
    pub fn view_projection_matrix(&self) -> Mat4 {
        self.projection_matrix() * self.view_matrix()
    }

    /// Get camera uniform data for GPU
    pub fn to_uniform(&self) -> CameraUniform {
        CameraUniform {
            view_proj: self.view_projection_matrix().to_cols_array_2d(),
            camera_pos: [self.position().x, self.position().y, self.position().z, 1.0],
        }
    }

    /// Orbit the camera by delta angles
    pub fn orbit(&mut self, delta_azimuth: f32, delta_elevation: f32) {
        self.azimuth += delta_azimuth;
        self.elevation = (self.elevation + delta_elevation)
            .clamp(-std::f32::consts::FRAC_PI_2 + 0.1, std::f32::consts::FRAC_PI_2 - 0.1);
    }

    /// Zoom by changing distance
    pub fn zoom(&mut self, delta: f32) {
        self.distance = (self.distance * (1.0 - delta * 0.1)).clamp(5.0, 50.0);
    }

    /// Update camera with auto-rotation
    pub fn update(&mut self, delta_time: f32) {
        if self.auto_rotate_speed != 0.0 {
            self.azimuth += self.auto_rotate_speed * delta_time;
        }
    }

    /// Reset camera to default position
    pub fn reset(&mut self) {
        self.azimuth = 0.0;
        self.elevation = 0.3;
        self.distance = 15.0;
    }

    /// Update aspect ratio (call on window resize)
    pub fn set_aspect(&mut self, aspect: f32) {
        self.aspect = aspect;
    }
}
