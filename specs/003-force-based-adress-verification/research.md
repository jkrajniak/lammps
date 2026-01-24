# Research: Force-based AdResS Verification (Dimer System)

**Feature**: 003-force-based-adress-verification  
**Date**: 2025-01-24

## Technical Decisions

### Decision 1: Python Version and Library Constraints

**Decision**: Use Python 3.8+ with only standard libraries (numpy, matplotlib, subprocess)

**Rationale**: 
- Kaggle and Google Colab default to Python 3.8+ with numpy and matplotlib pre-installed
- No external dependencies needed - ensures maximum compatibility
- subprocess is part of standard library for executing LAMMPS binary

**Alternatives Considered**:
- Using pandas for data handling: Rejected - numpy is sufficient for array operations
- Using specialized LAMMPS Python wrappers: Rejected - adds dependencies, subprocess is sufficient
- Using pathlib instead of os.path: Considered but os.path is more universally compatible

### Decision 2: LAMMPS Data File Format

**Decision**: Generate standard LAMMPS data file format with explicit `Atoms`, `Bonds`, `Masses` sections

**Rationale**:
- Standard format ensures compatibility with all LAMMPS versions
- Explicit sections make file human-readable and debuggable
- Follows LAMMPS documentation conventions

**Alternatives Considered**:
- Using LAMMPS Python interface: Rejected - not available in cloud environments, requires compilation
- Binary format: Rejected - not human-readable, harder to debug

### Decision 3: Dimer System Geometry

**Decision**: Use elongated box (30 x 10 x 10) with lattice-based molecule placement

**Rationale**:
- Elongated box provides clear spatial separation for AdResS zones (Atomistic center, Hybrid flanks, CG edges)
- Lattice placement prevents initial overlaps and ensures proper density
- Random orientations after placement provide realistic initial conditions

**Alternatives Considered**:
- Cubic box: Rejected - insufficient spatial separation for zone definition
- Random placement: Rejected - risk of overlaps, harder to control initial density

### Decision 4: Thermodynamic Force Table Format

**Decision**: Use LAMMPS table format (N header lines, then x f(x) pairs)

**Rationale**:
- Standard LAMMPS table format is well-documented
- Easy to read/write from Python
- Compatible with LAMMPS `fix` commands that read tables

**Alternatives Considered**:
- JSON format: Rejected - LAMMPS doesn't natively read JSON tables
- Binary format: Rejected - harder to debug and modify manually

### Decision 5: Iterative Boltzmann Inversion Algorithm

**Decision**: Use simple update formula: `Force_new(x) = Force_old(x) - PreFactor * (Density(x) - Target)`

**Rationale**:
- Standard Iterative Boltzmann Inversion (IBI) approach for AdResS
- Simple and robust convergence behavior
- PreFactor can be tuned for stability (typically 0.1-1.0)

**Alternatives Considered**:
- More sophisticated optimization (e.g., gradient descent): Rejected - adds complexity, simple IBI is sufficient for verification
- Adaptive PreFactor: Considered but deferred - fixed PreFactor is simpler for initial implementation

### Decision 6: Density Profile Output Method

**Decision**: Use LAMMPS `fix ave/chunk` to compute spatial density profiles

**Rationale**:
- Native LAMMPS capability, no custom code needed
- Provides spatial resolution along specified axis (x-axis for elongated box)
- Output format is straightforward to parse from Python

**Alternatives Considered**:
- Custom compute: Rejected - adds complexity, `fix ave/chunk` is sufficient
- Post-processing dump files: Rejected - less efficient, requires parsing large dump files

### Decision 7: Convergence Criteria

**Decision**: Run fixed number of iterations (e.g., 20) and track convergence via plots

**Rationale**:
- Simple and predictable execution time
- Plots provide visual verification of convergence
- Can be extended later with automatic convergence detection

**Alternatives Considered**:
- Automatic convergence detection: Considered but deferred - fixed iterations simpler for initial implementation
- Statistical convergence tests: Considered but deferred - visual inspection sufficient for verification

## Unresolved Questions

None - all technical decisions resolved.

## References

- LAMMPS Documentation: https://docs.lammps.org/
- AdResS Methodology: Force-based adaptive resolution schemes
- Iterative Boltzmann Inversion: Standard technique for AdResS thermodynamic force optimization
