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
| 1 | Foundation | Months 1-3 | Current |
| 2 | Mechanics | Months 4-6 | Pending |
| 3 | Core Metabolism | Months 7-9 | Pending |
| 4 | Oxygen Transport | Months 10-11 | Pending |
| 5 | Integration | Months 12-14 | Pending |
| 6 | Extended Biochemistry | Months 15-17 | Pending |
| 7 | Disease Models | Months 18-20 | Pending |
| 8 | Polish & Release | Months 21-24 | Pending |

## Example Output

```
> /phase status

=== Cell Simulator X Development Status ===

Current Phase: 1 - Foundation
Progress: 2/6 tasks (33%)

Completed:
  ✓ Project setup (Rust/Metal/WebGPU toolchain)
  ✓ Core data structures (CellState)

In Progress:
  → RBC mesh generation with Fung-Tong parametric surface

Remaining:
  ○ Basic rendering pipeline with camera controls
  ○ Spectrin network graph structure
  ○ Configuration system for parameters

Phase Deliverable: Rotating 3D RBC visualization with spectrin network overlay

Next milestone: Complete mesh generation to unblock rendering work
```

```
> /phase detail 2

=== Phase 2: Mechanics (Months 4-6) ===

Goal: Accurate mechanical simulation

Tasks:
  ○ DPD fluid solver (GPU-accelerated)
  ○ Membrane mechanics (Skalak model)
  ○ Spectrin elasticity (WLC model)
  ○ Deformation validation vs micropipette aspiration data
  ○ Shear flow simulation
  ○ Osmotic swelling/shrinking

Dependencies:
  - Requires Phase 1 complete (mesh, rendering)

Deliverable: Deformable RBC that matches experimental mechanics data

Validation Targets:
  - Shear modulus: 5.5 ± 1.5 μN/m
  - Match Dao et al. 2003 force-extension curves
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
