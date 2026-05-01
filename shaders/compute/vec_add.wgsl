// vec_add — Phase 11.0 sentinel kernel.
//
// Adds two f32 input buffers element-wise into an output buffer. Used by
// `--diagnose-gpu` to assert that the wgpu compute pipeline is wired up
// correctly and that CPU↔GPU correctness is bitwise.
//
// Workgroup size 64 is the Phase 11 default for 1-D dispatches: two
// subgroups of 32 on Apple Silicon (Metal), two warps on NVIDIA, two
// wavefronts on AMD.

@group(0) @binding(0) var<storage, read>       a: array<f32>;
@group(0) @binding(1) var<storage, read>       b: array<f32>;
@group(0) @binding(2) var<storage, read_write> out: array<f32>;

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    if (i >= arrayLength(&out)) {
        return;
    }
    out[i] = a[i] + b[i];
}
