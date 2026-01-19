# Implementation Plan: Phase 3: Atomistic Region COM Tracking for AdResS

**Branch**: `001-atomistic-com-tracking` | **Date**: 2025-01-27 | **Spec**: [spec.md](spec.md)
**Input**: Feature specification from `/specs/001-atomistic-com-tracking/spec.md`

## Summary

Implement center-of-mass (COM) calculation and CG particle tracking for molecules in the atomistic region of AdResS simulations. When molecules are in the atomistic region (lambda > threshold), the CG particle must accurately track the COM of the molecule's atomistic atoms. This enables proper force interpolation and smooth transitions between CG and atomistic regions.

**Technical Approach**: Complete the existing `calculate_com_from_atoms()` method (already has MPI), implement `update_cg_from_com()` to update CG particle from COM, and integrate into `post_force()` for atomistic region molecules. Handle periodic boundaries correctly using LAMMPS domain infrastructure.

## Technical Context

**Language/Version**: C++11 (LAMMPS standard)  
**Primary Dependencies**: LAMMPS core (atom, domain, memory, error, comm), MPI (for distributed molecule handling)  
**Storage**: In-memory arrays in `FixAdResSConstraint` class (x_com, v_com, mass_com)  
**Testing**: LAMMPS input scripts with verification, manual COM calculation comparison, MPI correctness tests  
**Target Platform**: Linux with MPI support (standard LAMMPS build)  
**Project Type**: LAMMPS user package extension (fix style)  
**Performance Goals**: COM calculation overhead < 5% of total simulation time for 1000+ molecules  
**Constraints**: Must maintain 1e-6 numerical precision, handle periodic boundaries correctly, work with existing Phase 2 CG constraints  
**Scale/Scope**: Systems with 1-10000 molecules, molecules with 2-100 atoms, distributed across 1-1000 MPI processors

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

### I. Scientific Accuracy First
- ✅ **PASS**: COM calculations use standard physics formulas (Σ(m_i * x_i) / Σ(m_i))
- ✅ **PASS**: Numerical precision requirement (1e-6) is explicitly stated
- ✅ **PASS**: Periodic boundary handling preserves physical correctness
- **Action**: Validate against manual calculations in testing

### II. Example-Driven Development
- ✅ **PASS**: Specification requires example test cases with verification scripts
- ✅ **PASS**: Success criteria include example test execution
- **Action**: Create example test cases in Phase 1 design

### III. Build Verification (NON-NEGOTIABLE)
- ✅ **PASS**: Success criteria explicitly require build success
- ✅ **PASS**: Constitution requires build verification after each task
- **Action**: Verify build after each implementation task

### IV. LAMMPS Integration Standards
- ✅ **PASS**: Uses existing LAMMPS infrastructure (MPI, domain, memory management)
- ✅ **PASS**: Follows LAMMPS fix style patterns
- ✅ **PASS**: Uses existing `calculate_com_from_atoms()` method structure
- **Action**: Follow LAMMPS coding conventions throughout

### V. Documentation and Analysis
- ✅ **PASS**: Specification requires documentation of design decisions
- ✅ **PASS**: Performance analysis required (overhead < 5%)
- **Action**: Document implementation details and verify performance

**Gate Status**: ✅ **PASS** - All constitution principles satisfied

## Project Structure

### Documentation (this feature)

```text
specs/001-atomistic-com-tracking/
├── plan.md              # This file (/speckit.plan command output)
├── research.md          # Phase 0 output (/speckit.plan command)
├── data-model.md        # Phase 1 output (/speckit.plan command)
├── quickstart.md        # Phase 1 output (/speckit.plan command)
├── contracts/           # Phase 1 output (/speckit.plan command)
└── tasks.md             # Phase 2 output (/speckit.tasks command - NOT created by /speckit.plan)
```

### Source Code (repository root)

```text
src/USER-ADRESS/
├── fix_adress_constraint.cpp    # Main implementation file
├── fix_adress_constraint.h      # Header file
└── IMPLEMENTATION_PLAN.md       # Overall project plan

examples/USER/adres/phase3_atomistic_com/
├── in.dimer_atomistic           # Dimer in atomistic region test
├── data.dimer_atomistic         # Data file
├── verify_com_tracking.py       # Verification script
└── README.md                    # Documentation
```

**Structure Decision**: Single LAMMPS user package extension. Code lives in `src/USER-ADRESS/` following LAMMPS conventions. Examples in `examples/USER/adres/phase3_atomistic_com/` following Phase 2 pattern.

## Complexity Tracking

> **No violations - all complexity justified by AdResS methodology requirements**

---

## Phase 0: Research Complete

**Status**: ✅ Complete

All research questions resolved. Key findings:
- Existing `calculate_com_from_atoms()` is correct and complete
- CG particle update requires direct assignment with image flag handling
- Region detection follows Phase 2 pattern (lambda-based)
- Integration extends existing `post_force()` method

See [research.md](research.md) for detailed findings.

---

## Phase 1: Design Complete

**Status**: ✅ Complete

Design artifacts created:
- **Data Model**: [data-model.md](data-model.md) - Entities, relationships, validation rules
- **Quickstart**: [quickstart.md](quickstart.md) - Usage guide and examples
- **API Contracts**: [contracts/fix-adress-constraint-api.md](contracts/fix-adress-constraint-api.md) - Method contracts and error handling

Key design decisions:
- Reuse existing `calculate_com_from_atoms()` method (no changes needed)
- Implement `update_cg_from_com()` with direct assignment and image flag handling
- Add `is_molecule_in_atomistic_region()` helper method
- Extend `post_force()` to handle atomistic region molecules

---

## Next Steps

Ready for Phase 2: Task breakdown via `/speckit.tasks` command.
