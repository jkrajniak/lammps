# Implementation Plan: Force-based AdResS Verification (Dimer System)

**Branch**: `003-force-based-adress-verification` | **Date**: 2025-01-24 | **Spec**: [spec.md](./spec.md)
**Input**: Feature specification from `/specs/003-force-based-adress-verification/spec.md`

**Note**: This template is filled in by the `/speckit.plan` command. See `.specify/templates/commands/plan.md` for the execution workflow.

## Summary

Implement a Python-based verification pipeline for Force-based AdResS in LAMMPS. The pipeline consists of three components: (1) a system generator (`gen_dimers.py`) that creates LAMMPS data files for dimer fluids, (2) a LAMMPS input script (`in.adress_dimer`) that runs AdResS simulations with thermodynamic force tables, and (3) an iterative verification manager (`run_verification.py`) that uses Iterative Boltzmann Inversion to tune the thermodynamic force until a flat density profile (target 0.844) is achieved across all AdResS zones. The system must run in cloud environments (Kaggle/Colab) using only standard Python libraries.

## Technical Context

**Language/Version**: Python 3.8+ (compatible with Kaggle/Colab default versions)  
**Primary Dependencies**: numpy, matplotlib, subprocess (standard library)  
**Storage**: File-based (LAMMPS data files, table files, dump files, density profile outputs)  
**Testing**: Manual verification via LAMMPS execution and density profile analysis  
**Target Platform**: Linux-based cloud kernels (Kaggle, Google Colab) - no GUI, headless execution  
**Project Type**: Single Python scripts (standalone utilities, no web/mobile components)  
**Performance Goals**: Verification loop should complete 20 iterations in reasonable time (< 1 hour for typical system sizes)  
**Constraints**: Must use only standard Python libraries (numpy, matplotlib, subprocess), no external dependencies. LAMMPS binary must be available as `lmp_serial`. Must work without GUI in cloud environment.  
**Scale/Scope**: Small to medium dimer systems (hundreds to thousands of molecules), 3 Python scripts, 1 LAMMPS input script

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

### I. Scientific Accuracy First
✅ **PASS**: Verification pipeline validates physical correctness of AdResS implementation by measuring density profiles and ensuring convergence to target density. Numerical accuracy maintained through proper LAMMPS input configuration.

### II. Example-Driven Development
✅ **PASS**: This feature IS an example/verification system. The generated scripts serve as both documentation and verification tools for AdResS implementation. Includes complete input scripts and verification procedures.

### III. Build Verification (NON-NEGOTIABLE)
⚠️ **NOTE**: This feature does not modify C++ code - it's a Python verification pipeline. LAMMPS binary (`lmp_serial`) must be pre-compiled. Python scripts must be syntactically correct and executable.

### IV. LAMMPS Integration Standards
✅ **PASS**: LAMMPS input script follows standard LAMMPS conventions (`units lj`, `atom_style molecular`, proper fix/compute usage). Uses standard LAMMPS data file format.

### V. Documentation and Analysis
✅ **PASS**: All design decisions documented in plan/research files. Verification results tracked through density profile plots and convergence analysis.

**Overall Status**: ✅ **PASS** - All gates satisfied. This is a verification/example system that supports AdResS development without modifying core LAMMPS code.

## Project Structure

### Documentation (this feature)

```text
specs/[###-feature]/
├── plan.md              # This file (/speckit.plan command output)
├── research.md          # Phase 0 output (/speckit.plan command)
├── data-model.md        # Phase 1 output (/speckit.plan command)
├── quickstart.md        # Phase 1 output (/speckit.plan command)
├── contracts/           # Phase 1 output (/speckit.plan command)
└── tasks.md             # Phase 2 output (/speckit.tasks command - NOT created by /speckit.plan)
```

### Source Code (repository root)

```text
examples/USER/adres/force_verification/
├── gen_dimers.py          # System generator script
├── in.adress_dimer        # LAMMPS input script
├── run_verification.py    # Iterative verification manager
├── system.data            # Generated LAMMPS data file (output)
├── thermo_force.table     # Thermodynamic force table (generated/updated)
├── density_profile.dat    # Density profile output (generated)
└── plots/                 # Density profile plots (generated)
    └── iteration_*.png
```

**Structure Decision**: Standalone Python scripts in `examples/USER/adres/force_verification/` directory. This follows LAMMPS convention of placing example scripts in `examples/USER/` subdirectories. All scripts are self-contained and can be run independently. Output files are generated in the same directory for simplicity in cloud environments.

## Complexity Tracking

> **Fill ONLY if Constitution Check has violations that must be justified**

No violations - this is a straightforward verification pipeline with three standalone scripts. Complexity is minimal and justified by the verification requirements.
