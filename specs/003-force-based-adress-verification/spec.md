# Feature Specification: Force-based AdResS Verification (Dimer System)

**Feature Branch**: `003-force-based-adress-verification`  
**Created**: 2025-01-24  
**Status**: Draft  
**Input**: User description: "I am implementing a verification pipeline for a custom C++ implementation of the Force-based Adaptive Resolution Scheme (AdResS) in LAMMPS. I need you to plan and generate the necessary Python scripts and LAMMPS input files to run this verification on a cloud environment (Kaggle/Colab)."

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Generate Dimer System Configuration (Priority: P1)

When a user needs to verify AdResS implementation, they must be able to generate a LAMMPS data file for a dimer fluid system with proper topology (2 atoms per molecule, bonds, masses).

**Why this priority**: This is the foundation for all verification. Without a properly configured dimer system, no verification can proceed. The system generator must create a valid LAMMPS data file that can be used for AdResS simulations.

**Independent Test**: Can be fully tested by running `gen_dimers.py` and verifying the output `system.data` file contains correct number of atoms, bonds, masses, and proper LAMMPS data format. Delivers a ready-to-use dimer system configuration.

**Acceptance Scenarios**:

1. **Given** a request to generate a dimer system, **When** `gen_dimers.py` is executed, **Then** it produces `system.data` with correct atom count, bond count, and LAMMPS data format
2. **Given** a generated dimer system, **When** the data file is loaded into LAMMPS, **Then** it parses without errors and reports correct number of molecules
3. **Given** a generated dimer system, **When** molecules are examined, **Then** each molecule has exactly 2 atoms connected by a bond with length approximately 0.97

---

### User Story 2 - Run AdResS Simulation with Thermodynamic Force (Priority: P1)

When a user needs to verify AdResS implementation, they must be able to run a LAMMPS simulation that applies a thermodynamic force from a table file to achieve flat density profile across zones.

**Why this priority**: This is the core verification mechanism. The LAMMPS input script must correctly set up AdResS zones (Atomistic, Hybrid, CG) and apply the thermodynamic force from the table file. Without this, the iterative verification loop cannot function.

**Independent Test**: Can be fully tested by running LAMMPS with `in.adress_dimer` and a `thermo_force.table` file, then verifying that the simulation completes and produces density profile output. Delivers a working AdResS simulation setup.

**Acceptance Scenarios**:

1. **Given** a `system.data` file and `thermo_force.table`, **When** LAMMPS is run with `in.adress_dimer`, **Then** the simulation completes without errors
2. **Given** a running AdResS simulation, **When** the simulation outputs density profiles, **Then** the output contains spatial density data for all zones
3. **Given** a LAMMPS input script, **When** it is executed, **Then** it correctly reads the thermodynamic force table and applies it to the appropriate zones

---

### User Story 3 - Iterative Density Profile Optimization (Priority: P1)

When a user needs to verify AdResS implementation, they must be able to iteratively tune the thermodynamic force to achieve a flat density profile (target density 0.844) across all zones using Iterative Boltzmann Inversion.

**Why this priority**: This is the verification goal - achieving a flat density profile proves the AdResS implementation is working correctly. The iterative loop must converge the density profile to the target value, demonstrating that the thermodynamic force is correctly applied.

**Independent Test**: Can be fully tested by running `run_verification.py` for multiple iterations and verifying that density profiles converge toward the target value (0.844) and that the thermodynamic force table is updated correctly. Delivers automated verification of AdResS correctness.

**Acceptance Scenarios**:

1. **Given** an initial zero thermodynamic force table, **When** `run_verification.py` runs for N iterations, **Then** the density profile converges toward target density (0.844)
2. **Given** a running verification loop, **When** each iteration completes, **Then** the thermodynamic force table is updated using the formula: `Force_new(x) = Force_old(x) - PreFactor * (Density(x) - Target)`
3. **Given** a verification run, **When** density profiles are plotted, **Then** the plots show convergence progress over iterations

---

### Edge Cases

- What happens when LAMMPS binary is not found? (System should provide clear error message)
- How does system handle malformed table files? (System should validate table format before running LAMMPS)
- What happens when density profile output is missing or corrupted? (System should detect and report errors)
- How does system handle convergence failures? (System should detect non-convergence and report)
- What happens when PreFactor is too large causing instability? (System should handle numerical stability)
- How does system handle different box geometries? (System should work with elongated boxes as specified)
- What happens when molecules cross zone boundaries during simulation? (LAMMPS should handle correctly)

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: System MUST generate LAMMPS data file (`system.data`) for dimer fluid with 2 atoms per molecule, bond type 1, mass 1.0 per atom, bond length 0.97
- **FR-002**: System MUST place molecule centers on a lattice to avoid overlaps, then assign random orientations
- **FR-003**: System MUST create elongated box geometry (e.g., 30 x 10 x 10) for AdResS zone separation
- **FR-004**: System MUST output standard LAMMPS data format with `Atoms`, `Bonds`, `Masses` sections
- **FR-005**: LAMMPS input script MUST use `units lj`, `atom_style molecular`, `boundary p p p`
- **FR-006**: LAMMPS input script MUST define regions for Atomistic (center), Hybrid (flanks), and CG (edges) zones
- **FR-007**: LAMMPS input script MUST use `bond_style harmonic` with k=500, r0=0.97
- **FR-008**: LAMMPS input script MUST use `pair_style lj/cut` with epsilon=1.0, sigma=1.0, cutoff=2.5
- **FR-009**: LAMMPS input script MUST exclude intra-molecular interactions (1-2 exclusion)
- **FR-010**: LAMMPS input script MUST read thermodynamic force from `thermo_force.table` file
- **FR-011**: LAMMPS input script MUST output dump file every 1000 steps
- **FR-012**: LAMMPS input script MUST output spatial average (density profile) using `fix ave/chunk`
- **FR-013**: Verification script MUST create initial `thermo_force.table` with all zeros
- **FR-014**: Verification script MUST execute LAMMPS using subprocess
- **FR-015**: Verification script MUST read density profile output from simulation
- **FR-016**: Verification script MUST compare measured density vs. target density (0.844)
- **FR-017**: Verification script MUST update thermodynamic force using formula: `Force_new(x) = Force_old(x) - PreFactor * (Density(x) - Target)`
- **FR-018**: Verification script MUST save updated `thermo_force.table` for next iteration
- **FR-019**: Verification script MUST generate density profile plots (optional .png) using Matplotlib to track progress
- **FR-020**: System MUST work in Linux-based cloud environment (Kaggle) with standard Python libraries only (numpy, matplotlib, subprocess)
- **FR-021**: System MUST use LAMMPS binary compiled as `lmp_serial`

### Key Entities *(include if feature involves data)*

- **Dimer System**: Represents a fluid of diatomic molecules with 2 atoms per molecule, connected by harmonic bonds
- **Thermodynamic Force Table**: Spatial force field read by LAMMPS to achieve flat density profile, updated iteratively
- **Density Profile**: Spatial distribution of density across simulation box, measured and compared to target (0.844)
- **AdResS Zones**: Three spatial regions - Atomistic (center), Hybrid (flanks), CG (edges) - where different resolution applies

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: `gen_dimers.py` generates valid LAMMPS data file that LAMMPS can parse without errors
- **SC-002**: `in.adress_dimer` runs successfully with provided `thermo_force.table` and produces density profile output
- **SC-003**: `run_verification.py` converges density profile to within 5% of target (0.844) within 20 iterations
- **SC-004**: All scripts run successfully in Kaggle/Colab environment using only standard Python libraries
- **SC-005**: Density profile plots show clear convergence trend over iterations
