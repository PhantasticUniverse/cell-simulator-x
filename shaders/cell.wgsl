// Cell Simulator X - WGSL Shader
// Renders RBC mesh with Phong shading and spectrin network as lines

// Uniforms
struct CameraUniform {
    view_proj: mat4x4<f32>,
    camera_pos: vec4<f32>,
};

struct RenderSettings {
    light_dir: vec4<f32>,
    ambient: vec4<f32>,
    diffuse: vec4<f32>,
    specular: vec4<f32>,  // w component is shininess
};

@group(0) @binding(0)
var<uniform> camera: CameraUniform;

@group(0) @binding(1)
var<uniform> settings: RenderSettings;

// Vertex input for mesh
struct MeshVertexInput {
    @location(0) position: vec3<f32>,
    @location(1) normal: vec3<f32>,
    @location(2) uv: vec2<f32>,
};

// Vertex output for mesh
struct MeshVertexOutput {
    @builtin(position) clip_position: vec4<f32>,
    @location(0) world_position: vec3<f32>,
    @location(1) world_normal: vec3<f32>,
    @location(2) uv: vec2<f32>,
};

// Mesh vertex shader
@vertex
fn vs_main(in: MeshVertexInput) -> MeshVertexOutput {
    var out: MeshVertexOutput;

    out.world_position = in.position;
    out.world_normal = normalize(in.normal);
    out.uv = in.uv;
    out.clip_position = camera.view_proj * vec4<f32>(in.position, 1.0);

    return out;
}

// Mesh fragment shader with Phong lighting
@fragment
fn fs_main(in: MeshVertexOutput) -> @location(0) vec4<f32> {
    let N = normalize(in.world_normal);
    let L = normalize(settings.light_dir.xyz);
    let V = normalize(camera.camera_pos.xyz - in.world_position);
    let H = normalize(L + V);

    // Ambient
    let ambient = settings.ambient.rgb;

    // Diffuse (Lambertian)
    let NdotL = max(dot(N, L), 0.0);
    let diffuse = settings.diffuse.rgb * NdotL;

    // Specular (Blinn-Phong)
    let NdotH = max(dot(N, H), 0.0);
    let shininess = settings.specular.w;
    let specular = settings.specular.rgb * pow(NdotH, shininess) * NdotL;

    // Subsurface scattering approximation for RBC translucency
    // Light passing through the thin membrane
    let transmit = max(dot(-N, L), 0.0) * 0.15;
    let sss = settings.diffuse.rgb * transmit;

    // Rim lighting for edge visibility
    let rim = pow(1.0 - max(dot(N, V), 0.0), 2.0) * 0.2;
    let rim_color = vec3<f32>(0.9, 0.3, 0.3) * rim;

    // Combine all lighting
    var color = ambient + diffuse + specular + sss + rim_color;

    // Add subtle variation based on position (oxygen gradient simulation)
    let oxygen_factor = 0.5 + 0.5 * sin(in.world_position.x * 2.0 + in.world_position.y * 2.0);
    color = mix(color, color * vec3<f32>(1.0, 0.85, 0.85), oxygen_factor * 0.1);

    return vec4<f32>(color, 1.0);
}

// Line vertex input for spectrin network
struct LineVertexInput {
    @location(0) position: vec3<f32>,
    @location(3) color: vec4<f32>,
};

// Line vertex output
struct LineVertexOutput {
    @builtin(position) clip_position: vec4<f32>,
    @location(0) color: vec4<f32>,
};

// Line vertex shader
@vertex
fn vs_line(in: LineVertexInput) -> LineVertexOutput {
    var out: LineVertexOutput;

    out.clip_position = camera.view_proj * vec4<f32>(in.position, 1.0);
    out.color = in.color;

    return out;
}

// Line fragment shader
@fragment
fn fs_line(in: LineVertexOutput) -> @location(0) vec4<f32> {
    return in.color;
}
