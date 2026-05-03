# Cell Simulator X — Headline Figures

Reproduction recipes for each figure that supports the preprint. Every
figure has a single command that produces a CSV under `target/` and a
follow-on plotting step (the plotting scripts are TBD; CSVs are
provided so any plotting tool can be used).

All commands assume the project is built (`cargo build --release`).

---

## Figure 1: 42-day storage lesion (the headline)

The quantitative metabolomics curve over 42 days of refrigerated
(4°C) storage in CPD baseline. Reproduces Hess 2010 ion gradients
at days 0/14/42 simultaneously, ATP per Hess 21-day half-life, DPG per
Zimrin 2009, deformability per Pivkin 2011.

```bash
cargo test --release --test storage_curve write_csv_succeeds -- --nocapture
# Writes /tmp/storage_curve_test.csv (43 rows: day 0..=42)
# Or run the full curve test for stdout output:
cargo test --release --test storage_curve ion_gradients_match_hess_2010_quantitatively -- --nocapture
```

Or via the headline diagnostic CLI:

```bash
cargo run --release -- --diagnose-disease storage --disease-param 42
```

CSV columns: `day, atp_mM, dpg23_mM, gsh_mM, gssg_mM, h2o2_mM,
nadph_mM, na_cyt_mM, k_cyt_mM, lactate_mM, pump_efficiency,
leak_multiplier, oxidative_stress, deformability_relative`.

Expected day-42 values: ATP ≈ 0.50, DPG ≈ 0.0, Na ≈ 59.7, K ≈ 90.6,
deformability ≈ 0.73.

---

## Figure 2: Storage parameter sweep (Phase 14.D)

Three additive solutions × ATP retention curves over 42 days, plus
supercooled (-4°C / 100-day) extension.

```bash
# Each preset's day-42 ATP within 15% of literature anchor
cargo test --release --test storage_additive_comparator -- --nocapture

# Day-42 ATP retention by additive:
#   CPD:    0.50 mM (25%)
#   SAGM:   1.00 mM (50%)
#   AS-3:   1.40 mM (70%)
#   PAGGSM: 1.70 mM (85%)
```

Sensitivity table (±20% perturbation per envelope parameter):

```bash
cargo test --release --test storage_sensitivity baseline_oat_sweep_emits_csv -- --nocapture
# Writes target/storage_sensitivity.csv (10 rows)
```

ATP half-life dominates day-42 ATP and deformability sensitivity by
~3 orders of magnitude over the other envelope parameters — robustness
finding.

---

## Figure 3: Tank-treading frequency vs shear rate (Phase 12.C.1)

Keller-Skalak 1982 analytic prediction matched against Fischer 2007
reported physiological-viscosity range.

```bash
cargo test --features validation --release \
    validation::experiments::fischer_2007 -- --nocapture
```

Predicted f_TT ≈ 0.075·γ̇ Hz (with γ̇ in 1/s) for canonical RBC
semi-axes 4×1 μm, falling in Fischer's K(λ) range 0.04–0.15.

---

## Figure 4: Parachute aspect ratio under Poiseuille flow (Phase 12.C.2)

Capillary number + empirical AR(Ca) fit.

```bash
cargo test --features validation --release \
    validation::experiments::skalak_1973 -- --nocapture
```

Predicted AR: 1.40 / 1.66 / 2.13 at γ̇_wall ∈ {100, 200, 400} 1/s,
within Skalak 1973's reported 1.5–2.0 band.

---

## Figure 5: Splenic transit at storage day N (Phase C-hybrid — the demo)

The integrated coupling demo: cells at storage days 0, 7, 14, 21, 28,
35, 42 traversing 0.5/0.7/1.0 μm splenic slits.

```bash
cargo test --release --test splenic_transit_sweep transit_sweep_emits_csv -- --nocapture
# Writes target/splenic_transit_storage_curve.csv (21 rows: 7 days × 3 widths)
```

Or for a single (day, width):

```bash
cargo run --release -- --diagnose-splenic-transit --storage-day 21 --slit-width 0.7
```

CSV columns: `storage_day, slit_width_um, wall_shear_rate_per_sec,
transit_time_sec, centroid_displacement_um, completed,
peak_strain_relative, peak_velocity_um_per_sec,
deformability_relative, wall_clock_ms`.

Slit-width effect is monotonic (0.5 μm: 1.27 μm displacement; 1.0 μm:
1.63 μm). Storage-day effect is muted with current
`SpectrinModulator` 0.5x stiffening cap — see
`docs/phase_c_hybrid_notes.md` for the detailed limitation discussion
and follow-on options.

---

## Figure 6: GPU vs CPU parity validation

Per-species relative error after 1 s of simulated time across the
38-metabolite ODE + the 6-dispatch integrated mechanics step. This
isn't a "preprint figure" but is the technical validation that the
GPU port preserves the science.

```bash
cargo test --release --test biochem_gpu_full_solver_parity -- --nocapture
cargo test --release --test physics_backend_parity -- --nocapture
```

Worst-fit relative errors:
- Biochemistry (38 species): 0.0023% (Na+ after 1 s).
- Mechanics (Skalak/WLC/DPD/Verlet integrated): 1.9e-15 μm position
  drift after 10 substeps (CPU≈GPU bit-near-perfect).

---

## Reproducibility

A `Makefile` (Phase 17.2) provides one-line targets per figure:

```bash
make figure-storage         # Figure 1
make figure-additives       # Figure 2 part A
make figure-sensitivity     # Figure 2 part B
make figure-tank-treading   # Figure 3
make figure-parachute       # Figure 4
make figure-splenic-transit # Figure 5
make figures                # all of the above
```

Each target writes a CSV to `target/`. Plotting is left to the consumer
(matplotlib, gnuplot, R, etc.); the plotting scripts will be added in
the preprint repo.

To reproduce the entire paper artifact:

```bash
git clone <repo>
cd cell-simulator-x
cargo build --release
make figures
# All CSVs land in target/
```

For Nix users: `nix develop && make figures`. (`flake.nix` is
deferred to a follow-on phase; `Dockerfile` is the universal
fallback.)
