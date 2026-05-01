//! GPU compute diagnostics (Phase 11.0).
//!
//! Exposes [`vec_add_gpu`], the CPU-vs-GPU correctness sentinel, and
//! [`run_diagnose_gpu`], the entry point for the `--diagnose-gpu` CLI flag.
//!
//! `vec_add_gpu` dispatches a trivial element-wise add on the GPU and
//! reads the result back. It is intentionally simple: any failure means
//! the wgpu pipeline plumbing is broken, not the simulation.

use std::time::Instant;

use anyhow::Result;
use wgpu::util::DeviceExt;

use super::ComputeContext;

/// Source of the `vec_add` compute kernel. Compiled at runtime by naga.
const VEC_ADD_WGSL: &str = include_str!("../../shaders/compute/vec_add.wgsl");

/// Compute `out[i] = a[i] + b[i]` on the GPU and read results back.
///
/// Returns the output vector on the host. Panics if `a.len() != b.len()`.
/// Used by [`run_diagnose_gpu`] and by `tests/compute_smoke.rs` as a
/// CPU-vs-GPU correctness check.
pub fn vec_add_gpu(ctx: &ComputeContext, a: &[f32], b: &[f32]) -> Result<Vec<f32>> {
    assert_eq!(a.len(), b.len(), "vec_add_gpu: length mismatch");
    let n = a.len();
    let bytes = (n * std::mem::size_of::<f32>()) as u64;

    // Storage buffers for inputs.
    let buf_a = ctx.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("vec_add a"),
        contents: bytemuck::cast_slice(a),
        usage: wgpu::BufferUsages::STORAGE,
    });
    let buf_b = ctx.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("vec_add b"),
        contents: bytemuck::cast_slice(b),
        usage: wgpu::BufferUsages::STORAGE,
    });

    // Output storage buffer (initialized to zero, gets written by the kernel).
    let buf_out = ctx.device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("vec_add out"),
        size: bytes,
        usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
        mapped_at_creation: false,
    });

    // Staging buffer for readback.
    let buf_staging = ctx.device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("vec_add staging"),
        size: bytes,
        usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
        mapped_at_creation: false,
    });

    // Bind group layout matches the WGSL @group(0) bindings.
    let bgl = ctx
        .device
        .create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("vec_add bgl"),
            entries: &[
                read_only_storage_entry(0),
                read_only_storage_entry(1),
                read_write_storage_entry(2),
            ],
        });

    let bind_group = ctx.device.create_bind_group(&wgpu::BindGroupDescriptor {
        label: Some("vec_add bg"),
        layout: &bgl,
        entries: &[
            wgpu::BindGroupEntry { binding: 0, resource: buf_a.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 1, resource: buf_b.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 2, resource: buf_out.as_entire_binding() },
        ],
    });

    let shader = ctx.device.create_shader_module(wgpu::ShaderModuleDescriptor {
        label: Some("vec_add shader"),
        source: wgpu::ShaderSource::Wgsl(VEC_ADD_WGSL.into()),
    });

    let pipeline_layout = ctx
        .device
        .create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("vec_add pl"),
            bind_group_layouts: &[&bgl],
            push_constant_ranges: &[],
        });

    let pipeline = ctx
        .device
        .create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("vec_add pipeline"),
            layout: Some(&pipeline_layout),
            module: &shader,
            entry_point: "main",
            compilation_options: Default::default(),
        });

    let mut encoder = ctx.device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
        label: Some("vec_add encoder"),
    });

    {
        let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
            label: Some("vec_add pass"),
            timestamp_writes: None,
        });
        pass.set_pipeline(&pipeline);
        pass.set_bind_group(0, &bind_group, &[]);
        // Workgroup size in WGSL is 64; dispatch ceil(n / 64) groups.
        let groups = ((n as u32) + 63) / 64;
        pass.dispatch_workgroups(groups, 1, 1);
    }

    encoder.copy_buffer_to_buffer(&buf_out, 0, &buf_staging, 0, bytes);
    ctx.queue.submit(std::iter::once(encoder.finish()));

    // Map the staging buffer and copy back. wgpu requires a poll on the
    // device to drive the mapping completion; on headless contexts we
    // just call `device.poll(Wait)`.
    let slice = buf_staging.slice(..);
    let (tx, rx) = std::sync::mpsc::channel();
    slice.map_async(wgpu::MapMode::Read, move |result| {
        let _ = tx.send(result);
    });
    ctx.device.poll(wgpu::Maintain::Wait);
    rx.recv()
        .map_err(|e| anyhow::anyhow!("vec_add: map_async channel closed: {e}"))?
        .map_err(|e| anyhow::anyhow!("vec_add: map_async failed: {e:?}"))?;

    let mapped = slice.get_mapped_range();
    let result: Vec<f32> = bytemuck::cast_slice(&mapped).to_vec();
    drop(mapped);
    buf_staging.unmap();

    Ok(result)
}

fn read_only_storage_entry(binding: u32) -> wgpu::BindGroupLayoutEntry {
    wgpu::BindGroupLayoutEntry {
        binding,
        visibility: wgpu::ShaderStages::COMPUTE,
        ty: wgpu::BindingType::Buffer {
            ty: wgpu::BufferBindingType::Storage { read_only: true },
            has_dynamic_offset: false,
            min_binding_size: None,
        },
        count: None,
    }
}

fn read_write_storage_entry(binding: u32) -> wgpu::BindGroupLayoutEntry {
    wgpu::BindGroupLayoutEntry {
        binding,
        visibility: wgpu::ShaderStages::COMPUTE,
        ty: wgpu::BindingType::Buffer {
            ty: wgpu::BufferBindingType::Storage { read_only: false },
            has_dynamic_offset: false,
            min_binding_size: None,
        },
        count: None,
    }
}

/// `--diagnose-gpu` entry point.
///
/// Builds a headless `ComputeContext`, dispatches `vec_add` on 1024 f32
/// elements, asserts CPU=GPU bitwise on every element, and prints a
/// short report. Returns `Err` if anything fails (so the binary exits
/// non-zero on regression).
pub fn run_diagnose_gpu() -> Result<()> {
    println!("=== Cell Simulator X — GPU Compute Diagnostic (Phase 11.0) ===\n");

    let setup_start = Instant::now();
    let ctx = ComputeContext::new_headless_blocking()?;
    let setup_elapsed = setup_start.elapsed();
    println!("Adapter:    {}", ctx.adapter_summary());
    println!("Setup:      {:.2?}", setup_elapsed);

    // Reproducible test vectors (no RNG so failures are deterministic).
    let n = 1024;
    let a: Vec<f32> = (0..n).map(|i| i as f32).collect();
    let b: Vec<f32> = (0..n).map(|i| (n - i) as f32 * 0.5).collect();

    let cpu_expected: Vec<f32> = a.iter().zip(b.iter()).map(|(x, y)| x + y).collect();

    let dispatch_start = Instant::now();
    let gpu_result = vec_add_gpu(&ctx, &a, &b)?;
    let dispatch_elapsed = dispatch_start.elapsed();
    println!("Dispatch:   {:.2?} (1024 f32 elements)", dispatch_elapsed);

    // Bitwise equality check. vec_add is purely additive so f32 add is
    // commutative-but-non-associative; for this kernel each output is a
    // single add of two CPU-floats, so the GPU result should be bitwise
    // identical to `a[i] + b[i]` computed on the host.
    let mut mismatches = 0usize;
    let mut first_diff = None;
    for i in 0..n {
        if gpu_result[i].to_bits() != cpu_expected[i].to_bits() {
            if first_diff.is_none() {
                first_diff = Some((i, cpu_expected[i], gpu_result[i]));
            }
            mismatches += 1;
        }
    }

    if mismatches == 0 {
        println!("Result:     CPU = GPU bitwise on all {} elements ✓", n);
        println!("\n=== Phase 11.0 sentinel: PASS ===");
        Ok(())
    } else {
        if let Some((i, cpu, gpu)) = first_diff {
            println!("Result:     {} mismatches; first at i={}: CPU={:.6} GPU={:.6}",
                mismatches, i, cpu, gpu);
        }
        println!("\n=== Phase 11.0 sentinel: FAIL ===");
        Err(anyhow::anyhow!(
            "vec_add CPU/GPU mismatch on {} of {} elements", mismatches, n
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn vec_add_cpu_gpu_bitwise() {
        // This test runs on the developer's machine (skipped on CI without GPU).
        // It dispatches vec_add and asserts bit-exact CPU=GPU correctness.
        let ctx = match ComputeContext::new_headless_blocking() {
            Ok(c) => c,
            Err(e) => {
                eprintln!("skipping GPU test: no adapter available: {e}");
                return;
            }
        };

        let n = 256;
        let a: Vec<f32> = (0..n).map(|i| i as f32 * 1.5).collect();
        let b: Vec<f32> = (0..n).map(|i| (i * i) as f32).collect();

        let gpu = vec_add_gpu(&ctx, &a, &b).expect("dispatch");
        for i in 0..n {
            let expected = a[i] + b[i];
            assert_eq!(
                gpu[i].to_bits(),
                expected.to_bits(),
                "i={}: CPU={} GPU={}",
                i, expected, gpu[i],
            );
        }
    }
}
