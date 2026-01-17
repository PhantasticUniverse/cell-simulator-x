//! WebGPU rendering pipeline for RBC visualization.
//!
//! Sets up wgpu device, surface, and render pipelines for
//! drawing the cell mesh and spectrin network overlay.

use std::sync::Arc;
use std::time::Instant;

use anyhow::Result;
use bytemuck::{Pod, Zeroable};
use wgpu::util::DeviceExt;
use winit::{dpi::PhysicalSize, window::Window};

use super::camera::Camera;
use crate::state::CellState;

/// Vertex for spectrin network lines
#[repr(C)]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
struct LineVertex {
    position: [f32; 3],
    color: [f32; 4],
}

/// Render settings uniform
#[repr(C)]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
struct RenderSettings {
    /// Light direction (normalized)
    light_dir: [f32; 4],
    /// Ambient color
    ambient: [f32; 4],
    /// Diffuse color
    diffuse: [f32; 4],
    /// Specular color and shininess
    specular: [f32; 4],
}

impl Default for RenderSettings {
    fn default() -> Self {
        Self {
            light_dir: [0.5, 0.8, 0.3, 0.0],
            ambient: [0.15, 0.08, 0.08, 1.0],   // Dark red ambient
            diffuse: [0.85, 0.2, 0.2, 1.0],     // Red diffuse (hemoglobin color)
            specular: [1.0, 0.9, 0.9, 32.0],    // White specular, shininess in w
        }
    }
}

/// Main render state managing all GPU resources
pub struct RenderState {
    #[allow(dead_code)]
    window: Arc<Window>,
    surface: wgpu::Surface<'static>,
    device: wgpu::Device,
    queue: wgpu::Queue,
    config: wgpu::SurfaceConfiguration,
    pub size: PhysicalSize<u32>,

    // Pipelines
    mesh_pipeline: wgpu::RenderPipeline,
    line_pipeline: wgpu::RenderPipeline,

    // Buffers
    mesh_vertex_buffer: wgpu::Buffer,
    mesh_index_buffer: wgpu::Buffer,
    mesh_index_count: u32,
    line_vertex_buffer: wgpu::Buffer,
    line_vertex_count: u32,

    // Uniforms
    camera_buffer: wgpu::Buffer,
    #[allow(dead_code)]
    settings_buffer: wgpu::Buffer,
    bind_group: wgpu::BindGroup,

    // Depth buffer
    depth_texture: wgpu::TextureView,

    // State
    pub camera: Camera,
    show_spectrin: bool,
    last_frame_time: Instant,
}

impl RenderState {
    /// Create new render state from window and cell state
    pub async fn new(window: Arc<Window>, cell_state: &CellState) -> Result<Self> {
        let size = window.inner_size();

        // Create wgpu instance
        let instance = wgpu::Instance::new(wgpu::InstanceDescriptor {
            backends: wgpu::Backends::all(),
            ..Default::default()
        });

        // Create surface
        // SAFETY: The Arc<Window> is stored in RenderState, ensuring the window
        // outlives the surface. The surface is dropped before the window.
        let surface = instance.create_surface(window.clone())?;

        // Get adapter
        let adapter = instance
            .request_adapter(&wgpu::RequestAdapterOptions {
                power_preference: wgpu::PowerPreference::HighPerformance,
                compatible_surface: Some(&surface),
                force_fallback_adapter: false,
            })
            .await
            .ok_or_else(|| anyhow::anyhow!("Failed to find suitable GPU adapter"))?;

        log::info!("Using adapter: {:?}", adapter.get_info());

        // Create device and queue
        let (device, queue) = adapter
            .request_device(
                &wgpu::DeviceDescriptor {
                    label: Some("RBC Device"),
                    required_features: wgpu::Features::empty(),
                    required_limits: wgpu::Limits::default(),
                },
                None,
            )
            .await?;

        // Configure surface
        let surface_caps = surface.get_capabilities(&adapter);
        let surface_format = surface_caps
            .formats
            .iter()
            .find(|f| f.is_srgb())
            .copied()
            .unwrap_or(surface_caps.formats[0]);

        let config = wgpu::SurfaceConfiguration {
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT,
            format: surface_format,
            width: size.width,
            height: size.height,
            present_mode: wgpu::PresentMode::AutoVsync,
            alpha_mode: surface_caps.alpha_modes[0],
            view_formats: vec![],
            desired_maximum_frame_latency: 2,
        };
        surface.configure(&device, &config);

        // Load shaders
        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("Cell Shader"),
            source: wgpu::ShaderSource::Wgsl(include_str!("../../shaders/cell.wgsl").into()),
        });

        // Create camera
        let camera = Camera::new(size.width as f32 / size.height as f32);

        // Create uniform buffers
        let camera_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("Camera Buffer"),
            contents: bytemuck::cast_slice(&[camera.to_uniform()]),
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        });

        let settings = RenderSettings::default();
        let settings_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("Settings Buffer"),
            contents: bytemuck::cast_slice(&[settings]),
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        });

        // Create bind group layout
        let bind_group_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("Uniform Bind Group Layout"),
            entries: &[
                wgpu::BindGroupLayoutEntry {
                    binding: 0,
                    visibility: wgpu::ShaderStages::VERTEX | wgpu::ShaderStages::FRAGMENT,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Uniform,
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
                wgpu::BindGroupLayoutEntry {
                    binding: 1,
                    visibility: wgpu::ShaderStages::FRAGMENT,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Uniform,
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
            ],
        });

        // Create bind group
        let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("Uniform Bind Group"),
            layout: &bind_group_layout,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: camera_buffer.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: settings_buffer.as_entire_binding(),
                },
            ],
        });

        // Create pipeline layout
        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("Render Pipeline Layout"),
            bind_group_layouts: &[&bind_group_layout],
            push_constant_ranges: &[],
        });

        // Create mesh pipeline
        let mesh_pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some("Mesh Pipeline"),
            layout: Some(&pipeline_layout),
            vertex: wgpu::VertexState {
                module: &shader,
                entry_point: "vs_main",
                buffers: &[wgpu::VertexBufferLayout {
                    array_stride: std::mem::size_of::<crate::geometry::Vertex>() as wgpu::BufferAddress,
                    step_mode: wgpu::VertexStepMode::Vertex,
                    attributes: &[
                        wgpu::VertexAttribute {
                            offset: 0,
                            shader_location: 0,
                            format: wgpu::VertexFormat::Float32x3,
                        },
                        wgpu::VertexAttribute {
                            offset: 12,
                            shader_location: 1,
                            format: wgpu::VertexFormat::Float32x3,
                        },
                        wgpu::VertexAttribute {
                            offset: 24,
                            shader_location: 2,
                            format: wgpu::VertexFormat::Float32x2,
                        },
                    ],
                }],
            },
            fragment: Some(wgpu::FragmentState {
                module: &shader,
                entry_point: "fs_main",
                targets: &[Some(wgpu::ColorTargetState {
                    format: config.format,
                    blend: Some(wgpu::BlendState::REPLACE),
                    write_mask: wgpu::ColorWrites::ALL,
                })],
            }),
            primitive: wgpu::PrimitiveState {
                topology: wgpu::PrimitiveTopology::TriangleList,
                strip_index_format: None,
                front_face: wgpu::FrontFace::Ccw,
                cull_mode: Some(wgpu::Face::Back),
                polygon_mode: wgpu::PolygonMode::Fill,
                ..Default::default()
            },
            depth_stencil: Some(wgpu::DepthStencilState {
                format: wgpu::TextureFormat::Depth32Float,
                depth_write_enabled: true,
                depth_compare: wgpu::CompareFunction::Less,
                stencil: wgpu::StencilState::default(),
                bias: wgpu::DepthBiasState::default(),
            }),
            multisample: wgpu::MultisampleState::default(),
            multiview: None,
        });

        // Create line pipeline for spectrin network
        let line_pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some("Line Pipeline"),
            layout: Some(&pipeline_layout),
            vertex: wgpu::VertexState {
                module: &shader,
                entry_point: "vs_line",
                buffers: &[wgpu::VertexBufferLayout {
                    array_stride: std::mem::size_of::<LineVertex>() as wgpu::BufferAddress,
                    step_mode: wgpu::VertexStepMode::Vertex,
                    attributes: &[
                        wgpu::VertexAttribute {
                            offset: 0,
                            shader_location: 0,
                            format: wgpu::VertexFormat::Float32x3,
                        },
                        wgpu::VertexAttribute {
                            offset: 12,
                            shader_location: 3,
                            format: wgpu::VertexFormat::Float32x4,
                        },
                    ],
                }],
            },
            fragment: Some(wgpu::FragmentState {
                module: &shader,
                entry_point: "fs_line",
                targets: &[Some(wgpu::ColorTargetState {
                    format: config.format,
                    blend: Some(wgpu::BlendState::ALPHA_BLENDING),
                    write_mask: wgpu::ColorWrites::ALL,
                })],
            }),
            primitive: wgpu::PrimitiveState {
                topology: wgpu::PrimitiveTopology::LineList,
                strip_index_format: None,
                front_face: wgpu::FrontFace::Ccw,
                cull_mode: None,
                polygon_mode: wgpu::PolygonMode::Fill,
                ..Default::default()
            },
            depth_stencil: Some(wgpu::DepthStencilState {
                format: wgpu::TextureFormat::Depth32Float,
                depth_write_enabled: false, // Lines don't write depth
                depth_compare: wgpu::CompareFunction::Less,
                stencil: wgpu::StencilState::default(),
                bias: wgpu::DepthBiasState::default(),
            }),
            multisample: wgpu::MultisampleState::default(),
            multiview: None,
        });

        // Create mesh buffers
        // Use VERTEX | COPY_DST to allow dynamic updates for physics simulation
        let mesh = &cell_state.geometry.mesh;
        let mesh_vertex_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("Mesh Vertex Buffer"),
            contents: bytemuck::cast_slice(&mesh.vertices),
            usage: wgpu::BufferUsages::VERTEX | wgpu::BufferUsages::COPY_DST,
        });

        let mesh_index_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("Mesh Index Buffer"),
            contents: bytemuck::cast_slice(&mesh.indices),
            usage: wgpu::BufferUsages::INDEX,
        });

        // Create spectrin line buffer
        let spectrin = &cell_state.geometry.spectrin_network;
        let mut line_vertices: Vec<LineVertex> = Vec::new();
        let line_color = [0.2, 0.8, 0.3, 0.7]; // Green, semi-transparent

        for edge in &spectrin.edges {
            let pos_a = spectrin.nodes[edge.node_a as usize].position;
            let pos_b = spectrin.nodes[edge.node_b as usize].position;

            // Offset slightly outward from surface for visibility
            let offset = 0.02; // μm
            let dir_a = glam::Vec3::from_array(pos_a).normalize();
            let dir_b = glam::Vec3::from_array(pos_b).normalize();

            line_vertices.push(LineVertex {
                position: (glam::Vec3::from_array(pos_a) + dir_a * offset).to_array(),
                color: line_color,
            });
            line_vertices.push(LineVertex {
                position: (glam::Vec3::from_array(pos_b) + dir_b * offset).to_array(),
                color: line_color,
            });
        }

        let line_vertex_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("Line Vertex Buffer"),
            contents: bytemuck::cast_slice(&line_vertices),
            usage: wgpu::BufferUsages::VERTEX | wgpu::BufferUsages::COPY_DST,
        });

        // Create depth texture
        let depth_texture = Self::create_depth_texture(&device, &config);

        Ok(Self {
            window,
            surface,
            device,
            queue,
            config,
            size,
            mesh_pipeline,
            line_pipeline,
            mesh_vertex_buffer,
            mesh_index_buffer,
            mesh_index_count: mesh.indices.len() as u32,
            line_vertex_buffer,
            line_vertex_count: line_vertices.len() as u32,
            camera_buffer,
            settings_buffer,
            bind_group,
            depth_texture,
            camera,
            show_spectrin: true,
            last_frame_time: Instant::now(),
        })
    }

    fn create_depth_texture(
        device: &wgpu::Device,
        config: &wgpu::SurfaceConfiguration,
    ) -> wgpu::TextureView {
        let texture = device.create_texture(&wgpu::TextureDescriptor {
            label: Some("Depth Texture"),
            size: wgpu::Extent3d {
                width: config.width,
                height: config.height,
                depth_or_array_layers: 1,
            },
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format: wgpu::TextureFormat::Depth32Float,
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::TEXTURE_BINDING,
            view_formats: &[],
        });
        texture.create_view(&wgpu::TextureViewDescriptor::default())
    }

    /// Resize the render surface
    pub fn resize(&mut self, new_size: PhysicalSize<u32>) {
        if new_size.width > 0 && new_size.height > 0 {
            self.size = new_size;
            self.config.width = new_size.width;
            self.config.height = new_size.height;
            self.surface.configure(&self.device, &self.config);
            self.depth_texture = Self::create_depth_texture(&self.device, &self.config);
            self.camera.set_aspect(new_size.width as f32 / new_size.height as f32);
        }
    }

    /// Toggle spectrin network visibility
    pub fn set_show_spectrin(&mut self, show: bool) {
        self.show_spectrin = show;
    }

    /// Update mesh vertices for dynamic simulation
    ///
    /// This uploads new vertex data to the GPU buffer.
    pub fn update_mesh(&mut self, mesh: &crate::geometry::Mesh) {
        self.queue.write_buffer(
            &self.mesh_vertex_buffer,
            0,
            bytemuck::cast_slice(&mesh.vertices),
        );
    }

    /// Update spectrin network lines after physics simulation
    pub fn update_spectrin(&mut self, spectrin: &crate::geometry::SpectrinNetwork) {
        let mut line_vertices: Vec<LineVertex> = Vec::new();
        let line_color = [0.2, 0.8, 0.3, 0.7]; // Green, semi-transparent

        for edge in &spectrin.edges {
            let pos_a = spectrin.nodes[edge.node_a as usize].position;
            let pos_b = spectrin.nodes[edge.node_b as usize].position;

            // Offset slightly outward from surface for visibility
            let offset = 0.02; // μm
            let dir_a = glam::Vec3::from_array(pos_a).normalize();
            let dir_b = glam::Vec3::from_array(pos_b).normalize();

            line_vertices.push(LineVertex {
                position: (glam::Vec3::from_array(pos_a) + dir_a * offset).to_array(),
                color: line_color,
            });
            line_vertices.push(LineVertex {
                position: (glam::Vec3::from_array(pos_b) + dir_b * offset).to_array(),
                color: line_color,
            });
        }

        // Only update if we have the same number of vertices
        if line_vertices.len() as u32 == self.line_vertex_count {
            self.queue.write_buffer(
                &self.line_vertex_buffer,
                0,
                bytemuck::cast_slice(&line_vertices),
            );
        }
    }

    /// Update state (called each frame before render)
    pub fn update(&mut self) {
        let now = Instant::now();
        let delta_time = (now - self.last_frame_time).as_secs_f32();
        self.last_frame_time = now;

        // Update camera with auto-rotation
        self.camera.update(delta_time);

        // Update camera uniform
        self.queue.write_buffer(
            &self.camera_buffer,
            0,
            bytemuck::cast_slice(&[self.camera.to_uniform()]),
        );
    }

    /// Render a frame
    pub fn render(&mut self) -> Result<(), wgpu::SurfaceError> {
        let output = self.surface.get_current_texture()?;
        let view = output
            .texture
            .create_view(&wgpu::TextureViewDescriptor::default());

        let mut encoder = self
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("Render Encoder"),
            });

        {
            let mut render_pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("Render Pass"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &view,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Clear(wgpu::Color {
                            r: 0.02,
                            g: 0.02,
                            b: 0.05,
                            a: 1.0,
                        }),
                        store: wgpu::StoreOp::Store,
                    },
                })],
                depth_stencil_attachment: Some(wgpu::RenderPassDepthStencilAttachment {
                    view: &self.depth_texture,
                    depth_ops: Some(wgpu::Operations {
                        load: wgpu::LoadOp::Clear(1.0),
                        store: wgpu::StoreOp::Store,
                    }),
                    stencil_ops: None,
                }),
                timestamp_writes: None,
                occlusion_query_set: None,
            });

            // Draw mesh
            render_pass.set_pipeline(&self.mesh_pipeline);
            render_pass.set_bind_group(0, &self.bind_group, &[]);
            render_pass.set_vertex_buffer(0, self.mesh_vertex_buffer.slice(..));
            render_pass.set_index_buffer(self.mesh_index_buffer.slice(..), wgpu::IndexFormat::Uint32);
            render_pass.draw_indexed(0..self.mesh_index_count, 0, 0..1);

            // Draw spectrin network
            if self.show_spectrin && self.line_vertex_count > 0 {
                render_pass.set_pipeline(&self.line_pipeline);
                render_pass.set_bind_group(0, &self.bind_group, &[]);
                render_pass.set_vertex_buffer(0, self.line_vertex_buffer.slice(..));
                render_pass.draw(0..self.line_vertex_count, 0..1);
            }
        }

        self.queue.submit(std::iter::once(encoder.finish()));
        output.present();

        Ok(())
    }
}
