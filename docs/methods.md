# Cell Simulator X — Methods

This document is the technical methods reference for the cell-simulator-x
project, as a basis for the preprint. It complements `PRD.md`
(architectural intent) and the `docs/phase_*_notes.md` series
(per-phase implementation notes).

## Model overview

Cell Simulator X integrates four coupled subsystems that together
describe the human red blood cell at single-cell resolution:

1. **Membrane mechanics** — Skalak shear+area constitutive model on a
   triangulated mesh, plus a Marko-Siggia worm-like-chain (WLC)
   spectrin network anchored to the mesh, plus dissipative-particle-
   dynamics (DPD) thermal noise. Velocity-Verlet integration with a
   `accel_scale = 10.0` magic factor (kept under documentation).
2. **Core metabolism** — full glycolysis (11 enzymes), Rapoport-
   Luebering 2,3-BPG shunt (BPGM/BPGP), pentose phosphate pathway
   (G6PD/6PGL/6PGD/RPI/RPE/TKL1/TKL2/TAL), glutathione cycle (GR/GPx),
   ion homeostasis (Na+/K+-ATPase, PMCA, Gardos K+ channel), Piezo1
   tension-gated Ca²⁺ entry. 38 species, RK4 integration.
3. **Oxygen transport** — Adair 4-site hemoglobin allostery with
   pH/2,3-DPG sensitivity, post-RK4 Euler interleaved per ms.
4. **Mechano-metabolic coupling** — `TensionComputer` aggregates
   per-element membrane tension; `SpectrinModulator` projects ATP →
   spectrin stiffness (Manno 2002 fit); Piezo1 open probability is a
   Hill function of membrane tension.

GPU compute migration (Phase 11) provides bit-near-perfect parity
between CPU and GPU paths: 0.0023% worst-fit per-species relative
error after 1 s simulation across all 38 metabolites; 1.9e-15 μm
position drift after 10 substeps for the integrated mechanics step.

Disease overlays (Phase 7) modulate envelope parameters: storage
lesion (Hess 2010 envelopes for ATP/2,3-DPG/pump/leak), diabetic
hyperglycemia, malaria parasite stages, sickle cell HbS substitution.

## Key parameters

All physiological parameters cite primary literature. The dominant
sources:

- **Glycolysis enzymes**: Mulquiney & Kuchel 1999; Rapoport 1976;
  Beutler 1984 (Vmax/Km/Keq for HK, PFK, Aldolase, GAPDH, PGK, PYK).
- **2,3-BPG shunt**: Beutler 1984 (BPGM/BPGP).
- **Pentose phosphate**: Vora 1981; Mulquiney 1999 (G6PD, 6PGD, etc.).
- **Glutathione**: Beutler 1969; Carlsson 2005 (GR, GPx, basal H₂O₂).
- **Ion homeostasis**: Joiner 1990 (Na+/K+-ATPase Vmax = 0.055 mM/s);
  Lew 2007 (PMCA, Gardos).
- **Hemoglobin**: Imai 1981 / 1982 (Adair coefficients vs pH/DPG).
- **Mechanics**: Skalak 1973 (Skalak shear/area moduli);
  Liu 1987 (in-situ spectrin Lp = 20 nm); Marko-Siggia 1995 (WLC);
  Waugh & Evans 1979 (micropipette aspiration).
- **Storage envelopes**: Hess 2010 review (ATP 21-day half-life);
  Zimrin 2009 (2,3-DPG depletion); Pivkin 2011 (deformability decline).
- **Splenic transit**: Pivkin 2016, Picas 2013, Drasler 1987.

Parameter files live in `data/parameters/` as JSON.

## Validation

Phase 10 (`--features validation`) runs a published-curve fit suite
against:

| Experiment | Reference | What it tests | χ²/dof |
|------------|-----------|---------------|--------|
| Imai 1981 standard OEC | Imai 1981 Methods Enzymol | Hb saturation vs pO₂ | <1 |
| Imai Bohr shift | Imai 1981 | OEC pH dependence | <1 |
| Imai DPG shift | Imai 1981 | OEC 2,3-DPG dependence | <1 |
| Mulquiney metabolites | Mulquiney 1999 | Steady-state pool | <2 (2 species canary) |
| Rief 1999 spectrin Fx | Rief 1999 / Liu 1987 | WLC F-x | <1 |
| Waugh-Evans 1979 | Waugh & Evans 1979 | Micropipette modulus | <1 |
| Dao 2003 | Dao 2003 | Optical-tweezers shear modulus | <1 |
| Fischer 2007 (Phase 12.C.1) | Fischer 2007 | Tank-treading vs shear rate | <1 (Keller-Skalak prediction) |
| Skalak 1973 (Phase 12.C.2) | Skalak 1973 / Fung 1993 | Parachute aspect ratio | <1 (Ca + empirical fit) |

9/11 experiments pass; 2 (Mulquiney PPP flux fraction; Mulquiney
multi-species steady state) are pre-existing audit canaries
documenting parameter-identifiability findings — see
`docs/validation_report_v1.md`.

Outside the curve-fit suite, ~360 integration + ~210 unit tests cover
behavioral correctness across the 14-phase implementation.

## Headline scientific output

The 42-day storage lesion at physiological timescale (Phase 14)
quantitatively reproduces Hess 2010 metabolomics:

| Day | Sim Na (mM) | Hess Na | Sim K (mM) | Hess K | Sim ATP | Hess ATP |
|-----|-------------|---------|------------|--------|---------|----------|
| 0   | 10.0        | 10      | 144        | 140    | 2.02    | 2.0      |
| 14  | 25.7        | 25      | 127        | 120    | 1.50    | 1.5      |
| 42  | 59.7        | 60      | 91         | 90     | 0.50    | 0.5      |

Plus deformability decline 1.000 → 0.728 over 42 days per
Pivkin et al. 2011. The exponential pump-efficiency envelope
(eff(t) = 0.295 + 0.705·exp(-0.305·t), Phase 14.B''-calibrated) is
the closed-form fit that achieves all three day-anchors
simultaneously.

Phase 14.D extends with additive solution comparator (CPD / AS-3 /
SAGM / PAGGSM), supercooled (-4°C, Q10 = 2.5) storage, and ±20%
sensitivity sweep over envelope parameters.

Phase C-hybrid couples the storage trajectory to splenic-slit transit
(0.5–1.0 μm rectangular slot) — the integrated demo no competitor
simulator can produce.

## Reproducing the headline figures

See `docs/headline_figures.md`. Each figure has a dedicated reproduction
command via `make figure-<name>` (Phase 17.2) or by running the test
that emits its CSV.

## Architecture notes

Single-cell vs multi-cell: `World` (Phase 10.5) owns N independent
`Cell`s with rayon-parallel per-cell stepping. The integrated GPU
`PhysicsBackend` (Phase 11.3.E) handles single-cell mechanics with
persistent buffers; multi-cell GPU integration is Phase 13 territory.

External flow: `flow::Poiseuille` (Phase 12.A) models cylindrical-
channel flow; `flow::SplenicSlit` (Phase C-hybrid) models the
splenic slot. Drag forces are written to `PhysicsState::external_forces_uN`
(CPU) or via `PhysicsBackend::set_external_forces` (GPU); the
`skalak_init_from_baseline` kernel adds them per substep.

## Limitations

- **Storage-day mechanical signal is muted.** The Phase 8
  `SpectrinModulator` (max_stiffening_factor 0.5) doesn't translate
  to large transit-time differences in C-hybrid. Bending-modulus and
  membrane-viscosity coupling are documented follow-on work.
- **Multi-cell microcirculation (Phase 13) is deferred.** The
  current single-cell-in-flow scope captures the mechano-metabolic
  moat phenomena (Piezo1 firing, ATP release, splenic transit) but
  not collective rheology (Fåhræus-Lindqvist, hematocrit-dependent
  apparent viscosity).
- **Mulquiney PPP flux fraction.** A pre-existing identifiability
  finding; the simulator's PPP flux runs ~10× higher than literature
  steady-state. Documented in `docs/validation_report_v1.md`.
