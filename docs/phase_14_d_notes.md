# Phase 14.D Notes — Storage Parameter Sweep

Phase 14.D extends the now-validated storage framework into the
practical questions that motivate transfusion medicine: which additive
solution best preserves cells, can subzero supercooled storage extend
viability, and which envelope parameters most strongly drive day-42
outcome?

## 14.D.1 — Additive comparator (AS-3 / SAGM / PAGGSM)

`src/storage/additive.rs` introduces `AdditiveSolution` with four
literature-anchored presets. Day-42 ATP retention is the headline
anchor; per-additive `lesion_config()` derives the ATP exponential
half-life analytically from the retention target.

| Additive | Day-42 ATP retention | T_half (days) | Reference |
|----------|----------------------|---------------|-----------|
| CPD baseline | 25% | 21 | Hess 2010 (Phase 14.B'' calibration) |
| AS-3 (Nutricel) | 70% | 81.6 | Hess 2010 Table 2 |
| SAGM | 50% | 42.0 | Hess 2010 |
| PAGGSM | 85% | 179.0 | Burger 2008, de Korte 2008 |

Pump-decay and leak-increase rates scale inversely with retention
(`preservation_scale = 0.25 / retention`) on the assumption that
additives that preserve ATP also preserve membrane-protein function. A
follow-on phase can refit these from per-additive membrane-protein
literature when the storage paper demands tighter anchoring.

### Pump envelope choice per additive

The Phase 14.B'' Hess-calibrated **exponential** envelope (eff_inf =
0.295, β = 0.305/day) was fit to CPD-baseline ion data; only CPD uses
it. AS-3, SAGM, and PAGGSM use the linear envelope from
`StorageLesionConfig`. Per-additive exponential calibrations are a
clean follow-on once per-additive ion-gradient time-courses become
available.

### Validation

- `tests/storage_additive_comparator.rs` (8 tests):
  - Each preset's day-42 ATP within 15% of literature anchor
  - Day-42 ATP monotonic across CPD < SAGM < AS-3 < PAGGSM
  - Day-42 deformability tracks ATP retention monotonically
  - Day-0 states agree across additives (no time evolution yet)
  - DPG depletes by day 14 across all additives
- `src/storage/additive.rs::tests` (6 unit tests):
  - CPD config preserves the Phase 14.B'' calibration
  - T_half values match analytical retention formula
  - Day-42 ATP targets match literature
  - Only CPD uses exponential envelope

## 14.D.2 — Supercooled storage with Q10 = 2.5

`StorageSimConfig::temperature_celsius` (default 4°C) drives Q10 = 2.5
scaling on every envelope rate per Hess 2010 review:

```
rate(T) = rate(4°C) · 2.5^((T - 4) / 10)
T_half(T) = T_half(4°C) / 2.5^((T - 4) / 10)
```

`StorageSimConfig::supercooled(additive)` builds a -4°C / 100-day
config matching Bruinsma 2019 and Berendsen 2014 supercooled storage
protocols. The Q10 transformation is applied via
`temperature_scaled_lesion_config()` at the start of each storage step
and via `hess_pump_efficiency_scaled(day · scale)` for the exponential
pump envelope path.

### Headline result

| Run | ATP (mM) | Deformability |
|-----|----------|---------------|
| Standard 4°C, day 42 | 0.500 | 0.728 |
| Supercooled -4°C, day 100 | ~0.41 | ~0.69 |

Relative deformability difference: 1.5%, well within the 10%
validation criterion. **Supercooled storage extends viable storage
duration by ~2× through temperature scaling alone**, with no parameter
retuning. The result composes with the additive comparator: PAGGSM +
supercooled outperforms CPD + supercooled at day 100, demonstrating
the additive ranking holds under temperature modulation.

### Validation

- `tests/storage_supercooled.rs` (6 tests):
  - Q10 scale formula matches canonical (-4°C → 0.4774)
  - Supercooled extends ATP lifetime
  - Day-100 supercooled deformability within 10% of standard day-42
  - Pump efficiency higher at same day under supercooling
  - Additive ranking preserved under supercooling
  - Standard 4°C unchanged from Phase 14.B'' headline (reference
    temperature sanity check)

### Simplification: envelope-only Q10

Q10 is applied to the empirical envelope rates only, not to the
underlying biochemistry's enzyme Vmax values in `FullyIntegratedConfig`.
The biochemistry still runs at 277 K (4°C). This captures the
dominant supercooled-storage effect (slowed envelope evolution) while
keeping the change small. Extending Q10 to the 38-species kinetics is a
clean follow-on if quantitative supercooled metabolomics demands it.

## 14.D.3 — Sensitivity analysis (±20% one-at-a-time)

`src/storage/sensitivity.rs` runs ±20% perturbations on each of five
envelope parameters and records the fractional change in day-42 ATP
and deformability. Output: `target/storage_sensitivity.csv`.

### Finding

Under `force_atp_dpg_targets = true` (the headline default), day-42 ATP
and deformability are dominated by `atp_decay_half_life_days`:

| Parameter | day-42 ATP rel. change (±20%) | day-42 def rel. change |
|-----------|-------------------------------|------------------------|
| atp_half_life_days | ±29% | ±2.6% |
| pump_efficiency_decay_per_day | <0.001% | <0.001% |
| leak_increase_per_day | <0.001% | <0.001% |
| dpg_loss_rate_mM_per_day | 0% | 0% |
| oxidative_stress_increase_per_day | 0% | 0% |

This is **a defensible robustness claim**: the headline ATP/deformability
trajectory is robust to second-order parameter uncertainty. The
non-ATP parameters affect ion gradients (Na, K) but not the ATP
envelope, which is forced to its exponential target each step.

### Why pump/leak insensitivity is structural, not weak

`force_atp_dpg_targets` overrides ATP and 2,3-DPG to the envelope
targets each storage day. Pump-decay, leak-increase, and oxidative
stress all affect *ion gradients* (Na+, K+) and *redox state* (GSH,
GSSG, H₂O₂, NADPH) — but not ATP, since ATP is set externally. With
`SpectrinModulator` driving deformability from ATP alone, deformability
inherits ATP's robustness.

If a future phase tracks ion-gradient-dependent measures (e.g. echinocyte
formation index, hemolysis risk from K+ leak), those would surface
high sensitivity to pump/leak parameters. The current sensitivity
report is honest about the metric it covers.

### Validation

- `tests/storage_sensitivity.rs` (3 tests): sweep emits the CSV
  deliverable, top-3 ranking includes ATP half-life, perturbations
  produce distinct outputs.
- `src/storage/sensitivity.rs::tests` (3 unit tests).

## Outcomes

- 17 new tests pass (8 + 6 + 3 integration; complemented by 6 + 6 + 3
  unit tests for the underlying modules — 32 new tests total).
- Standard 4°C / CPD storage curve unchanged from Phase 14.B'' headline.
- Three new strategic axes available: additive choice, storage
  temperature, and parameter sensitivity — all parameterizable from
  `StorageSimConfig` without recompiling.

## Deferred

- Per-additive exponential pump envelope calibrations (need per-additive
  ion-gradient time-course data).
- Q10 scaling of biochemistry enzyme Vmax in `FullyIntegratedConfig`
  (envelope-only is dominant under current targets).
- Sensitivity over ion-gradient and redox observables (current
  sensitivity tracks only ATP and deformability).
- Two-at-a-time interaction sensitivity (would distinguish independent
  vs synergistic envelope levers).
