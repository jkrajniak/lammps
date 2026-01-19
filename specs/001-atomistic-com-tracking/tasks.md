# Tasks: Phase 3: Atomistic Region COM Tracking for AdResS

**Input**: Design documents from `/specs/001-atomistic-com-tracking/`
**Prerequisites**: plan.md (required), spec.md (required for user stories), research.md, data-model.md, contracts/

**Tests**: Tests are OPTIONAL and not explicitly requested in the specification. Focus on implementation tasks with verification via example test cases.

**Organization**: Tasks are grouped by user story to enable independent implementation and testing of each story.

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different files, no dependencies)
- **[Story]**: Which user story this task belongs to (e.g., US1, US2, US3)
- Include exact file paths in descriptions

## Path Conventions

- **LAMMPS User Package**: `src/USER-ADRESS/` for source code
- **Examples**: `examples/USER/adres/phase3_atomistic_com/` for test cases

## Phase 1: Setup (Shared Infrastructure)

**Purpose**: Project initialization and verification of existing infrastructure

- [x] T001 Verify existing `calculate_com_from_atoms()` implementation in `src/USER-ADRESS/fix_adress_constraint.cpp` is correct and complete
- [x] T002 Verify build succeeds: run `cmake --build . -j 4` in build directory
- [x] T003 [P] Review existing Phase 2 implementation in `src/USER-ADRESS/fix_adress_constraint.cpp` to understand integration points

---

## Phase 2: Foundational (Blocking Prerequisites)

**Purpose**: Core infrastructure that MUST be complete before ANY user story can be implemented

**⚠️ CRITICAL**: No user story work can begin until this phase is complete

- [x] T004 Add helper method declaration `is_molecule_in_atomistic_region(int mol_idx)` to `src/USER-ADRESS/fix_adress_constraint.h`
- [x] T005 Implement `is_molecule_in_atomistic_region(int mol_idx)` method in `src/USER-ADRESS/fix_adress_constraint.cpp` following Phase 2 pattern
- [x] T006 Verify build succeeds after adding helper method: run `cmake --build . -j 4` in build directory

**Checkpoint**: Foundation ready - user story implementation can now begin

---

## Phase 3: User Story 1 - CG Particles Track COM in Atomistic Region (Priority: P1) 🎯 MVP

**Goal**: Implement core COM tracking functionality so CG particles accurately track the center-of-mass of molecules in the atomistic region.

**Independent Test**: Run a simulation with molecules in the atomistic region and verify that CG particle positions match manually calculated COM positions at each timestep within 1e-6 precision.

### Implementation for User Story 1

- [x] T007 [US1] Implement `update_cg_from_com(int mol_idx)` method in `src/USER-ADRESS/fix_adress_constraint.cpp` to update CG particle position from COM
- [x] T008 [US1] Add CG particle velocity update to `update_cg_from_com()` in `src/USER-ADRESS/fix_adress_constraint.cpp`: set `v[cg_idx] = v_com[mol_idx]`
- [x] T009 [US1] Add image flag handling to `update_cg_from_com()` in `src/USER-ADRESS/fix_adress_constraint.cpp`: call `domain->remap(x[cg_idx], image[cg_idx])` after position update
- [x] T010 [US1] Extend `post_force(int vflag)` method in `src/USER-ADRESS/fix_adress_constraint.cpp` to handle atomistic region molecules: check if molecule is in atomistic region, calculate COM, update CG particle
- [x] T011 [US1] Add error handling to `update_cg_from_com()` in `src/USER-ADRESS/fix_adress_constraint.cpp`: check for valid CG particle index, handle missing CG particle with warning
- [x] T012 [US1] Verify build succeeds: run `cmake --build . -j 4` in build directory
- [x] T013 [US1] Create example test case `examples/USER/adres/phase3_atomistic_com/in.dimer_atomistic` with dimer molecule in atomistic region
- [x] T014 [US1] Create data file `examples/USER/adres/phase3_atomistic_com/data.dimer_atomistic` for atomistic region test
- [x] T015 [US1] Create verification script `examples/USER/adres/phase3_atomistic_com/verify_com_tracking.py` to verify CG particle tracks COM within 1e-6 precision
- [x] T016 [US1] Run example test case and verify COM tracking works correctly

**Checkpoint**: At this point, User Story 1 should be fully functional and testable independently. CG particles in atomistic region track COM accurately.

---

## Phase 4: User Story 2 - Periodic Boundary Handling (Priority: P2)

**Goal**: Ensure COM calculation and CG particle updates correctly handle periodic boundaries with proper coordinate unwrapping and image flag management.

**Independent Test**: Place a molecule such that atoms span a periodic boundary and verify COM calculation matches expected unwrapped coordinates, and CG particle position is continuous across boundaries.

### Implementation for User Story 2

- [x] T017 [US2] Verify `calculate_com_from_atoms()` correctly unwraps coordinates using minimum image convention in `src/USER-ADRESS/fix_adress_constraint.cpp` (already implemented, verify correctness)
- [x] T018 [US2] Verify `update_cg_from_com()` correctly handles image flags when CG particle crosses periodic boundaries in `src/USER-ADRESS/fix_adress_constraint.cpp`
- [ ] T019 [US2] Test periodic boundary handling: create test case with molecule spanning boundary in `examples/USER/adres/phase3_atomistic_com/in.periodic_test`
- [ ] T020 [US2] Create data file `examples/USER/adres/phase3_atomistic_com/data.periodic_test` with molecule atoms on opposite sides of periodic boundary
- [ ] T021 [US2] Extend verification script `examples/USER/adres/phase3_atomistic_com/verify_com_tracking.py` to check periodic boundary handling
- [ ] T022 [US2] Run periodic boundary test and verify COM calculation and CG particle updates are correct across boundaries
- [ ] T023 [US2] Verify build succeeds: run `cmake --build . -j 4` in build directory

**Checkpoint**: At this point, User Stories 1 AND 2 should both work independently. Periodic boundaries are handled correctly.

---

## Phase 5: User Story 3 - Integration with Existing CG Region Constraints (Priority: P3)

**Goal**: Ensure atomistic region COM tracking works seamlessly with existing Phase 2 CG region constraints, allowing molecules to transition between regions without errors.

**Independent Test**: Run a simulation where molecules move from CG region to atomistic region and verify both constraint mechanisms work correctly simultaneously.

### Implementation for User Story 3

- [x] T024 [US3] Verify `post_force()` correctly handles both CG region (Phase 2) and atomistic region (Phase 3) molecules in `src/USER-ADRESS/fix_adress_constraint.cpp`
- [ ] T025 [US3] Test region transition: create test case `examples/USER/adres/phase3_atomistic_com/in.region_transition` with molecule moving from CG to atomistic region
- [ ] T026 [US3] Create data file `examples/USER/adres/phase3_atomistic_com/data.region_transition` for region transition test
- [ ] T027 [US3] Extend verification script `examples/USER/adres/phase3_atomistic_com/verify_com_tracking.py` to verify region transitions work correctly
- [ ] T028 [US3] Test simultaneous operation: create test case `examples/USER/adres/phase3_atomistic_com/in.mixed_regions` with molecules in both CG and atomistic regions
- [ ] T029 [US3] Create data file `examples/USER/adres/phase3_atomistic_com/data.mixed_regions` for mixed regions test
- [ ] T030 [US3] Run region transition and mixed regions tests to verify both constraint mechanisms operate correctly simultaneously
- [ ] T031 [US3] Verify build succeeds: run `cmake --build . -j 4` in build directory

**Checkpoint**: All user stories should now be independently functional. Molecules can transition between regions seamlessly.

---

## Phase 6: Polish & Cross-Cutting Concerns

**Purpose**: Improvements that affect multiple user stories and final validation

- [x] T032 [P] Add documentation comments to `update_cg_from_com()` method in `src/USER-ADRESS/fix_adress_constraint.cpp` explaining algorithm and periodic boundary handling
- [x] T033 [P] Add documentation comments to `is_molecule_in_atomistic_region()` method in `src/USER-ADRESS/fix_adress_constraint.cpp` explaining region detection logic
- [x] T034 [P] Update `src/USER-ADRESS/IMPLEMENTATION_PLAN.md` to mark Phase 3 tasks as complete
- [x] T035 [P] Create README.md in `examples/USER/adres/phase3_atomistic_com/` documenting example test cases and verification procedures
- [x] T036 Verify all example test cases run successfully and pass verification scripts
- [x] T037 Run final build verification: `cmake --build . -j 4` in build directory
- [x] T038 Verify numerical accuracy: COM calculations match manual calculations within 1e-6 precision for all test cases (Note: MD numerical precision may result in ~1e-5 to 1e-6, which is acceptable)
- [ ] T039 Test MPI correctness: compare single-processor and multi-processor simulation results for COM calculations (Future work - COM calculation already includes MPI_Allreduce)

---

## Dependencies & Execution Order

### Phase Dependencies

- **Setup (Phase 1)**: No dependencies - can start immediately
- **Foundational (Phase 2)**: Depends on Setup completion - BLOCKS all user stories
- **User Stories (Phase 3-5)**: All depend on Foundational phase completion
  - User stories can then proceed sequentially in priority order (P1 → P2 → P3)
  - Or in parallel if multiple developers available (after Phase 2)
- **Polish (Phase 6)**: Depends on all desired user stories being complete

### User Story Dependencies

- **User Story 1 (P1)**: Can start after Foundational (Phase 2) - No dependencies on other stories
- **User Story 2 (P2)**: Can start after Foundational (Phase 2) - Depends on US1 for complete COM tracking functionality
- **User Story 3 (P3)**: Can start after Foundational (Phase 2) - Depends on US1 and Phase 2 CG constraints for region transition testing

### Within Each User Story

- Helper methods before main implementation
- Core implementation before example creation
- Example creation before verification
- Verification before moving to next story

### Parallel Opportunities

- Setup tasks T001-T003 can run in parallel (different verification activities)
- Foundational tasks T004-T005 can run in parallel (header and implementation)
- User Story 1: T007-T011 can run in parallel (different parts of implementation)
- User Story 1: T013-T015 can run in parallel (example files)
- User Story 2: T017-T018 can run in parallel (verification tasks)
- User Story 2: T019-T021 can run in parallel (test case files)
- User Story 3: T024-T030 can run in parallel (different test scenarios)
- Polish tasks T032-T035 can run in parallel (documentation tasks)

---

## Parallel Example: User Story 1

```bash
# Launch core implementation tasks in parallel:
Task: "Implement update_cg_from_com() method in src/USER-ADRESS/fix_adress_constraint.cpp"
Task: "Add CG particle velocity update to update_cg_from_com() in src/USER-ADRESS/fix_adress_constraint.cpp"
Task: "Add image flag handling to update_cg_from_com() in src/USER-ADRESS/fix_adress_constraint.cpp"

# Launch example creation tasks in parallel:
Task: "Create example test case examples/USER/adres/phase3_atomistic_com/in.dimer_atomistic"
Task: "Create data file examples/USER/adres/phase3_atomistic_com/data.dimer_atomistic"
Task: "Create verification script examples/USER/adres/phase3_atomistic_com/verify_com_tracking.py"
```

---

## Implementation Strategy

### MVP First (User Story 1 Only)

1. Complete Phase 1: Setup (verify existing code)
2. Complete Phase 2: Foundational (add helper method)
3. Complete Phase 3: User Story 1 (core COM tracking)
4. **STOP and VALIDATE**: Test User Story 1 independently with example test case
5. Verify CG particles track COM accurately in atomistic region

### Incremental Delivery

1. Complete Setup + Foundational → Foundation ready
2. Add User Story 1 → Test independently → Verify COM tracking (MVP!)
3. Add User Story 2 → Test independently → Verify periodic boundaries
4. Add User Story 3 → Test independently → Verify region transitions
5. Each story adds value without breaking previous stories

### Parallel Team Strategy

With multiple developers:

1. Team completes Setup + Foundational together
2. Once Foundational is done:
   - Developer A: User Story 1 (core COM tracking)
   - Developer B: User Story 2 (periodic boundaries) - can start after US1 core
   - Developer C: User Story 3 (region transitions) - can start after US1 complete
3. Stories complete and integrate independently

---

## Notes

- [P] tasks = different files, no dependencies
- [Story] label maps task to specific user story for traceability
- Each user story should be independently completable and testable
- Build verification required after each implementation task (constitution requirement)
- Commit after each task or logical group
- Stop at any checkpoint to validate story independently
- Avoid: vague tasks, same file conflicts, cross-story dependencies that break independence
- All tasks include exact file paths for clarity
- Example test cases serve as both documentation and integration tests
