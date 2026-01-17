# /benchmark - Performance Benchmarks

Run performance benchmarks for simulation components.

## Usage

```
/benchmark [target] [--profile] [--iterations N]
```

## Targets

### all (default)
Run complete benchmark suite.

### membrane
Membrane force calculations.
- Target: <0.5ms for 10K vertices

### dpd
DPD fluid solver timestep.
- Target: <1ms for 100K particles

### metabolism
ODE integration step.
- Target: <0.1ms for 200 reactions

### neighbor
Neighbor list rebuild.
- Target: <2ms for cell list

### render
Rendering frame time.
- Target: <16ms (60 fps)

## Options

- `--profile`: Generate detailed profiling data
- `--iterations N`: Number of iterations (default: 100)
- `--warmup N`: Warmup iterations (default: 10)
- `--export`: Export results to JSON

## Output

```
=== Cell Simulator X Benchmarks ===
Hardware: Apple M4 Max, 64GB RAM
Backend: Metal

Component          | Mean    | Std Dev | Target  | Status
-------------------|---------|---------|---------|--------
membrane_forces    | 0.42ms  | 0.03ms  | <0.5ms  | ✓
dpd_timestep       | 0.87ms  | 0.05ms  | <1.0ms  | ✓
metabolic_step     | 0.08ms  | 0.01ms  | <0.1ms  | ✓
neighbor_rebuild   | 1.65ms  | 0.12ms  | <2.0ms  | ✓
render_frame       | 12.3ms  | 0.8ms   | <16ms   | ✓

Total simulation step: 3.02ms (331 fps theoretical)
```

## Profiling Output

With `--profile`, generates:
- `benchmarks/profiles/membrane_profile.json`
- `benchmarks/profiles/gpu_timeline.json`
- Flame graph HTML (if `inferno` installed)

## Implementation

Runs: `cargo bench --features bench`

Uses criterion.rs for statistical analysis.

## Comparison Mode

Compare against previous results:

```
/benchmark --compare benchmarks/baseline.json
```

Output shows regression/improvement percentages.

## Example

```
/benchmark membrane --iterations 1000 --profile
```

```
=== Membrane Forces Benchmark ===
Iterations: 1000 (after 10 warmup)

Percentiles:
  p50: 0.41ms
  p95: 0.48ms
  p99: 0.52ms

Memory bandwidth: 12.4 GB/s
GPU occupancy: 68%

Profile saved to: benchmarks/profiles/membrane_20240115.json
```
