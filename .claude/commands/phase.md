# /phase - Development Phase Progress

Show and manage development phase progress based on PRD.md roadmap.

## Usage

```
/phase [action] [args]
```

## Actions

### status (default)
Show current phase status and progress.

```
/phase
/phase status
```

### list
List all phases with completion status.

```
/phase list
```

### detail
Show detailed tasks for a specific phase.

```
/phase detail 1
/phase detail "Foundation"
```

### mark
Mark a task as complete or in-progress.

```
/phase mark 1.3 complete
/phase mark "RBC mesh generation" complete
```

### next
Show next recommended tasks based on dependencies.

```
/phase next
```

## Phases Overview

| Phase | Name | Timeline | Status |
|-------|------|----------|--------|
| 1 | Foundation | Months 1-3 | ✅ Complete |
| 2 | Mechanics | Months 4-6 | ✅ Complete |
| 3 | Core Metabolism | Months 7-9 | Current |
| 4 | Oxygen Transport | Months 10-11 | Pending |
| 5 | Integration | Months 12-14 | Pending |
| 6 | Extended Biochemistry | Months 15-17 | Pending |
| 7 | Disease Models | Months 18-20 | Pending |
| 8 | Polish & Release | Months 21-24 | Pending |

## Example Output

```
> /phase status

=== Cell Simulator X Development Status ===

Current Phase: 3 - Core Metabolism
Progress: 0/6 tasks (0%)

Completed Phases:
  ✅ Phase 1: Foundation (6/6 tasks)
  ✅ Phase 2: Mechanics (5/6 tasks)

Next Tasks:
  ○ ODE solver (adaptive Runge-Kutta)
  ○ Glycolysis implementation (all 10 reactions)
  ○ Rapoport-Luebering shunt (2,3-DPG)
  ○ ATP consumption/production balance
  ○ Steady-state validation

Phase Deliverable: Validated metabolic model reproducing Joshi-Palsson results

Next milestone: Implement adaptive ODE solver for metabolic reactions
```

```
> /phase detail 2

=== Phase 2: Mechanics (Months 4-6) ✅ COMPLETE ===

Goal: Accurate mechanical simulation

Tasks:
  ✓ DPD fluid solver (CPU-based)
  ✓ Membrane mechanics (Skalak model)
  ✓ Spectrin elasticity (WLC model)
  ✓ Velocity-Verlet time integration
  ✓ Dynamic mesh rendering
  ○ Deformation validation vs micropipette data (ongoing)

Implementation:
  - src/physics/wlc.rs - WLC spectrin elasticity
  - src/physics/membrane.rs - Skalak strain energy
  - src/physics/dpd.rs - DPD fluid dynamics
  - src/physics/integrator.rs - Velocity-Verlet

Deliverable: Deformable RBC with physics simulation ✅

Validation Targets:
  - Shear modulus: 5.5 ± 1.1 μN/m ✓
  - Bending modulus: 0.18 pN·μm ✓
```

## Progress Tracking

Progress stored in `.claude/phase-progress.json`:

```json
{
  "current_phase": 1,
  "tasks": {
    "1.1": {"status": "complete", "date": "2024-01-10"},
    "1.2": {"status": "complete", "date": "2024-01-12"},
    "1.3": {"status": "in_progress", "started": "2024-01-13"}
  }
}
```

## Reference

Full phase details in PRD.md Section 7.
