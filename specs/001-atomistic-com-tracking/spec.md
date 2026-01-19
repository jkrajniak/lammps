# Feature Specification: Phase 3: Atomistic Region COM Tracking for AdResS

**Feature Branch**: `001-atomistic-com-tracking`  
**Created**: 2025-01-27  
**Status**: Draft  
**Input**: User description: "Implement Phase 3: Atomistic Region COM Tracking for AdResS"

## User Scenarios & Testing *(mandatory)*

### User Story 1 - CG Particles Track COM in Atomistic Region (Priority: P1)

When a molecule is in the atomistic region (lambda > threshold), the CG particle must accurately track the center-of-mass (COM) of the molecule's atomistic atoms. The CG particle position and velocity must be updated from the calculated COM after each timestep, ensuring the CG particle represents the molecule's collective motion.

**Why this priority**: This is the core functionality of Phase 3. Without accurate COM tracking, the CG particle cannot properly represent the molecule in the atomistic region, breaking the AdResS methodology. This enables proper force interpolation and smooth transitions between regions.

**Independent Test**: Can be fully tested by running a simulation with molecules in the atomistic region and verifying that CG particle positions match manually calculated COM positions at each timestep. Delivers accurate molecular representation in atomistic region.

**Acceptance Scenarios**:

1. **Given** a dimer molecule (2 atoms + 1 CG particle) in the atomistic region, **When** the simulation runs for multiple timesteps, **Then** the CG particle position matches the COM calculated from atom positions within numerical precision (1e-6)
2. **Given** a water molecule (3 atoms + 1 CG particle) in the atomistic region, **When** atoms move due to forces, **Then** the CG particle velocity matches the COM velocity calculated from atom velocities
3. **Given** molecules distributed across multiple MPI processors, **When** COM is calculated, **Then** the result is identical to a single-processor calculation

---

### User Story 2 - Periodic Boundary Handling (Priority: P2)

When molecules cross periodic boundaries, COM calculation must correctly unwrap coordinates and handle image flags to maintain accuracy across boundaries.

**Why this priority**: Periodic boundaries are fundamental to molecular dynamics simulations. Incorrect handling leads to incorrect COM positions and breaks the simulation. This ensures correctness in typical bulk system simulations.

**Independent Test**: Can be fully tested by placing a molecule such that atoms span a periodic boundary and verifying COM calculation matches expected unwrapped coordinates. Delivers correct behavior in periodic systems.

**Acceptance Scenarios**:

1. **Given** a molecule with atoms on opposite sides of a periodic boundary, **When** COM is calculated, **Then** coordinates are correctly unwrapped and COM is calculated accurately
2. **Given** a molecule crossing periodic boundaries during simulation, **When** CG particle is updated from COM, **Then** image flags are correctly handled and CG particle position is continuous

---

### User Story 3 - Integration with Existing CG Region Constraints (Priority: P3)

The atomistic region COM tracking must work seamlessly with existing CG region constraints (Phase 2), allowing molecules to transition between regions without errors.

**Why this priority**: Molecules move between regions during simulation. The system must handle both CG constraints (Phase 2) and atomistic COM tracking (Phase 3) correctly. This ensures smooth AdResS operation.

**Independent Test**: Can be fully tested by running a simulation where molecules move from CG region to atomistic region and verifying both constraint mechanisms work correctly. Delivers seamless region transitions.

**Acceptance Scenarios**:

1. **Given** a molecule transitioning from CG region to atomistic region, **When** lambda crosses the threshold, **Then** the system switches from CG constraints to atomistic COM tracking without errors
2. **Given** a simulation with molecules in both CG and atomistic regions, **When** the simulation runs, **Then** both constraint mechanisms operate correctly simultaneously

---

### Edge Cases

- What happens when a molecule has no CG particle? (System should handle gracefully with warning)
- How does system handle molecules split across many processors? (MPI communication must aggregate correctly)
- What happens when molecule mass is zero? (System should handle with error or warning)
- How does system handle molecules with only one atom? (COM should equal that atom's position)
- What happens when atoms in a molecule are very far apart (wrapping issues)? (Minimum image convention must handle correctly)

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: System MUST calculate center-of-mass (COM) position from atom positions for molecules in atomistic region using formula: `x_com = Σ(m_i * x_i) / Σ(m_i)`
- **FR-002**: System MUST calculate COM velocity from atom velocities for molecules in atomistic region using formula: `v_com = Σ(m_i * v_i) / Σ(m_i)`
- **FR-003**: System MUST aggregate COM calculations across all MPI processors for distributed molecules using MPI_Allreduce
- **FR-004**: System MUST update CG particle position to match calculated COM: `x_cg = x_com`
- **FR-005**: System MUST update CG particle velocity to match calculated COM velocity: `v_cg = v_com`
- **FR-006**: System MUST handle periodic boundaries correctly by unwrapping coordinates relative to a reference atom before COM calculation
- **FR-007**: System MUST handle image flags correctly when updating CG particle position across periodic boundaries
- **FR-008**: System MUST calculate COM in `post_force()` after forces have been computed and atoms have moved
- **FR-009**: System MUST only apply atomistic COM tracking to molecules with lambda > atomistic threshold (lambda_at_threshold, default 0.9)
- **FR-010**: System MUST correctly handle molecules with variable atom masses (both per-type and per-atom masses)
- **FR-011**: System MUST maintain numerical accuracy: COM calculations must match manual calculations within 1e-6 precision
- **FR-012**: System MUST handle edge case where molecule has no atoms (return zero COM with warning)

### Key Entities *(include if feature involves data)*

- **Molecule**: Represents a collection of atoms with a corresponding CG particle. Has molecule ID, associated atoms, CG particle index, and region classification (CG/transition/atomistic)
- **Center of Mass (COM)**: Calculated position and velocity representing the collective motion of a molecule's atoms. Has position (x_com[3]), velocity (v_com[3]), and total mass (mass_com)
- **CG Particle**: Dummy particle representing the molecule at coarse-grained resolution. Has position, velocity, and mass that must track COM in atomistic region

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: CG particle positions match manually calculated COM positions within 1e-6 precision for all molecules in atomistic region
- **SC-002**: CG particle velocities match manually calculated COM velocities within 1e-6 precision for all molecules in atomistic region
- **SC-003**: COM calculations produce identical results in single-processor and multi-processor (MPI) simulations
- **SC-004**: System correctly handles periodic boundaries: COM calculations remain accurate when molecules span boundaries
- **SC-005**: Molecules can transition between CG and atomistic regions without errors or discontinuities
- **SC-006**: Example test cases run successfully and pass automated verification scripts
- **SC-007**: Build succeeds without errors: `cmake --build . -j 4` completes successfully
- **SC-008**: Performance impact is acceptable: COM calculation overhead is less than 5% of total simulation time for systems with 1000+ molecules

## Assumptions

- Molecules are already mapped and CG particles are identified (from Phase 1 and Phase 2)
- Lambda values are correctly calculated by `fix adress/region` (from Phase 1)
- Atom masses are available via LAMMPS atom structure (mass[] or rmass[])
- MPI communication infrastructure is available and working
- Periodic boundary conditions are properly configured in LAMMPS domain
- Image flags are managed by LAMMPS domain decomposition

## Dependencies

- **Phase 1**: Lambda assignment and region management must be complete
- **Phase 2**: CG region constraints must be implemented (for testing region transitions)
- **LAMMPS Infrastructure**: Requires working MPI, domain decomposition, and atom data structures
- **Existing Code**: `calculate_com_from_atoms()` method exists but may need verification/completion

## Out of Scope

- Transition region interpolation (Phase 4)
- Bonded interactions support
- Time-dependent lambda changes
- Performance optimizations beyond basic MPI communication
- Energy conservation fixes (handled in Phase 2)
- New example creation (can reuse Phase 2 examples with atomistic region setup)
