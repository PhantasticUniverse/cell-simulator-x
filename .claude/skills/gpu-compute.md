# GPU Compute Domain Knowledge

## Metal (Apple Silicon)

### Architecture Priorities
1. Unified memory - zero-copy between CPU/GPU
2. Tile-based rendering - efficient for deferred
3. Threadgroup memory - 32KB per threadgroup
4. SIMD groups - 32 threads

### Compute Kernel Structure
```metal
kernel void membrane_forces(
    device float3* positions [[buffer(0)]],
    device float3* forces [[buffer(1)]],
    constant SimParams& params [[buffer(2)]],
    uint tid [[thread_position_in_grid]],
    uint lid [[thread_position_in_threadgroup]],
    threadgroup float3* shared_pos [[threadgroup(0)]]
) {
    // Tile loading into threadgroup memory
    // Force calculation
    // Atomic accumulation if needed
}
```

### Memory Optimization
- Use `MTLStorageMode.shared` for dynamic data
- Use `MTLStorageMode.private` for GPU-only
- Prefer `simd_sum()` over atomic adds
- Batch small dispatches

### Recommended Thread Counts
- Membrane vertices: 64 threads/group
- Particle simulation: 256 threads/group
- Reduction operations: 1024 threads/group

## CUDA (NVIDIA)

### Memory Hierarchy
1. Global (slow, large)
2. Shared (fast, 48KB/SM)
3. Registers (fastest, limited)
4. Constant (cached, read-only)

### Kernel Pattern
```cuda
__global__ void dpd_forces(
    float3* positions,
    float3* velocities,
    float3* forces,
    int n_particles,
    SimParams params
) {
    __shared__ float3 s_pos[BLOCK_SIZE];

    int tid = blockIdx.x * blockDim.x + threadIdx.x;

    // Tile-based neighbor interaction
    for (int tile = 0; tile < gridDim.x; tile++) {
        // Load tile to shared memory
        __syncthreads();
        // Compute interactions
        __syncthreads();
    }
}
```

### Occupancy Optimization
- Target 50%+ occupancy
- Balance registers vs threads
- Use `__launch_bounds__` for critical kernels

## WebGPU

### Compute Shader (WGSL)
```wgsl
@group(0) @binding(0) var<storage, read> positions: array<vec3f>;
@group(0) @binding(1) var<storage, read_write> forces: array<vec3f>;

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3u) {
    let idx = gid.x;
    if (idx >= arrayLength(&positions)) { return; }

    // Compute forces
    forces[idx] = compute_force(idx);
}
```

### Limitations vs Native
- No true shared memory atomics (emulate)
- Workgroup size limits (256 typical)
- Buffer size limits (check adapter)
- No subgroup operations (WebGPU 2.0)

## Common Patterns

### Parallel Reduction
```
// Sum N values
// Phase 1: Local reduction in threadgroup
// Phase 2: Global reduction across groups
```

### Neighbor Lists (for DPD)
- Cell list: O(N) build, O(N) query
- Verlet list: Rebuild when max displacement > skin/2
- Spatial hashing for uniform distributions

### Double Buffering
```rust
struct SimState {
    positions: [Buffer; 2],
    current: usize,
}

impl SimState {
    fn swap(&mut self) {
        self.current = 1 - self.current;
    }
    fn read_buffer(&self) -> &Buffer {
        &self.positions[self.current]
    }
    fn write_buffer(&self) -> &Buffer {
        &self.positions[1 - self.current]
    }
}
```

## Performance Targets

| Operation | Target | Notes |
|-----------|--------|-------|
| Membrane force calc | <0.5ms | 10K vertices |
| DPD timestep | <1ms | 100K particles |
| Neighbor rebuild | <2ms | Cell list |
| ODE step | <0.1ms | 200 reactions |

## Profiling Tools

### Metal
- Xcode GPU debugger
- Metal System Trace
- `MTLCaptureManager` for programmatic capture

### CUDA
- Nsight Compute
- nvprof / ncu
- cuda-memcheck

### WebGPU
- Chrome DevTools (GPU panel)
- RenderDoc (limited support)

## Cross-Platform Abstraction

```rust
pub trait ComputeBackend {
    fn dispatch_membrane_forces(&self, state: &mut CellState);
    fn dispatch_dpd_step(&self, state: &mut CellState, dt: f32);
    fn dispatch_metabolic_step(&self, state: &mut BiochemState, dt: f32);
    fn sync(&self);
}
```

Implementations: `MetalBackend`, `CudaBackend`, `WebGpuBackend`, `CpuFallback`

## Key Citations

- Apple Metal Best Practices Guide
- NVIDIA CUDA C++ Programming Guide
- WebGPU Specification (W3C)
- Fedosov GPU-DPD (2010) - GPU blood cell simulation
