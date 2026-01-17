//! Geometry benchmarks

use criterion::{black_box, criterion_group, criterion_main, Criterion};

use cell_simulator_x::config::GeometryParameters;
use cell_simulator_x::geometry::{Mesh, SpectrinNetwork};

fn bench_mesh_generation(c: &mut Criterion) {
    let params = GeometryParameters::default();

    c.bench_function("mesh_generation", |b| {
        b.iter(|| Mesh::generate_rbc(black_box(&params)))
    });
}

fn bench_spectrin_generation(c: &mut Criterion) {
    let params = GeometryParameters {
        spectrin_target_count: 1000, // Smaller for benchmarking
        ..Default::default()
    };
    let mesh = Mesh::generate_rbc(&params);

    c.bench_function("spectrin_generation", |b| {
        b.iter(|| SpectrinNetwork::generate(black_box(&mesh), black_box(&params)))
    });
}

fn bench_volume_calculation(c: &mut Criterion) {
    let params = GeometryParameters::default();
    let mesh = Mesh::generate_rbc(&params);

    c.bench_function("volume_calculation", |b| {
        b.iter(|| black_box(&mesh).calculate_volume())
    });
}

criterion_group!(
    benches,
    bench_mesh_generation,
    bench_spectrin_generation,
    bench_volume_calculation
);
criterion_main!(benches);
