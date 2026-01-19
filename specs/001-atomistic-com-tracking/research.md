# Research: Phase 3 Atomistic Region COM Tracking

**Date**: 2025-01-27  
**Feature**: Phase 3: Atomistic Region COM Tracking for AdResS

## Research Questions

### 1. COM Calculation Algorithm Verification

**Question**: Is the existing `calculate_com_from_atoms()` implementation correct for atomistic region COM tracking?

**Research**: Reviewed existing implementation in `fix_adress_constraint.cpp` lines 547-646.

**Findings**:
- ✅ Algorithm correctly implements: `x_com = Σ(m_i * x_i) / Σ(m_i)`
- ✅ Algorithm correctly implements: `v_com = Σ(m_i * v_i) / Σ(m_i)`
- ✅ MPI communication uses `MPI_Allreduce` correctly for distributed molecules
- ✅ Periodic boundary handling uses minimum image convention with reference atom
- ✅ Handles both per-type (`mass[]`) and per-atom (`rmass[]`) masses

**Decision**: Existing `calculate_com_from_atoms()` implementation is correct and complete. No changes needed to the algorithm itself.

**Rationale**: The implementation follows standard COM calculation formulas and correctly handles MPI and periodic boundaries. The code is production-ready.

**Alternatives Considered**: 
- Rewriting the method: Rejected - existing implementation is correct
- Using different reference atom selection: Current approach (CG particle or first atom) is appropriate

---

### 2. CG Particle Update from COM

**Question**: How should CG particle position and velocity be updated from COM in the atomistic region?

**Research**: Reviewed LAMMPS fix patterns and Phase 2 implementation.

**Findings**:
- CG particles are regular atoms in LAMMPS (not special dummy particles)
- Position update: `x_cg = x_com` (direct assignment)
- Velocity update: `v_cg = v_com` (direct assignment)
- Image flags must be handled using `domain->remap()` for periodic boundaries
- Update should happen in `post_force()` after atoms have moved

**Decision**: Implement `update_cg_from_com()` to:
1. Set CG particle position: `x[cg_idx] = x_com`
2. Set CG particle velocity: `v[cg_idx] = v_com`
3. Handle image flags: `domain->remap(x[cg_idx], image[cg_idx])`
4. Only update if CG particle exists and molecule is in atomistic region

**Rationale**: Direct assignment is correct since CG particles are regular atoms. Image flag handling ensures continuity across periodic boundaries.

**Alternatives Considered**:
- Using LAMMPS atom->set_position(): Rejected - direct assignment is simpler and correct
- Updating in different fix method: Rejected - `post_force()` is correct timing (after forces, before final_integrate)

---

### 3. Atomistic Region Detection

**Question**: How to determine if a molecule is in the atomistic region?

**Research**: Reviewed Phase 2 implementation and `fix adress/region` lambda calculation.

**Findings**:
- Lambda values are stored in `fix adress/region` and accessible via `fix_region->get_lambda(i)`
- Atomistic region threshold: `lambda_at_threshold` (default 0.9)
- Molecule is in atomistic region if: `lambda > lambda_at_threshold`
- Can check lambda of representative atom (CG particle or first atom)

**Decision**: Add helper method `is_molecule_in_atomistic_region(int mol_idx)` similar to existing `is_molecule_in_cg_region()`. Check lambda of representative atom and compare to threshold.

**Rationale**: Consistent with Phase 2 pattern. Lambda-based region detection is the standard AdResS approach.

**Alternatives Considered**:
- Position-based detection: Rejected - lambda is more accurate and already calculated
- Storing region flags: Rejected - lambda check is simple and accurate

---

### 4. Integration with post_force()

**Question**: How should atomistic COM tracking integrate with existing `post_force()` that handles CG region?

**Research**: Reviewed current `post_force()` implementation (lines 466-506).

**Findings**:
- Current `post_force()` only handles CG region molecules
- Need to add atomistic region handling
- Both can coexist: check region and apply appropriate logic
- COM calculation happens after forces, which is correct timing

**Decision**: Extend `post_force()` to:
1. Loop through all molecules
2. For CG region: apply existing CG constraints (already implemented)
3. For atomistic region: calculate COM from atoms, update CG particle from COM
4. For transition region: Skip (handled in Phase 4)

**Rationale**: Single method handles both regions efficiently. Clear separation of logic based on lambda threshold.

**Alternatives Considered**:
- Separate method for atomistic region: Rejected - `post_force()` is the correct integration point
- Different timing: Rejected - after forces is correct (atoms have moved, COM can be calculated)

---

### 5. Periodic Boundary Image Flag Handling

**Question**: How to correctly handle image flags when updating CG particle position across periodic boundaries?

**Research**: Reviewed LAMMPS domain documentation and Phase 2 implementation.

**Findings**:
- LAMMPS uses `image` flags to track periodic boundary crossings
- `domain->remap(x, image)` updates image flags based on position
- Must call remap after updating position
- Image flags are integers encoding periodic crossings in each dimension

**Decision**: After updating CG particle position from COM, call `domain->remap(atom->x[cg_idx], atom->image[cg_idx])` to ensure image flags are correct.

**Rationale**: Standard LAMMPS pattern for handling periodic boundaries. Ensures continuity and correct neighbor list updates.

**Alternatives Considered**:
- Manual image flag calculation: Rejected - `domain->remap()` is the standard LAMMPS method
- Not updating image flags: Rejected - would break periodic boundary handling

---

## Summary of Decisions

1. **COM Calculation**: Existing `calculate_com_from_atoms()` is correct - no changes needed
2. **CG Particle Update**: Implement `update_cg_from_com()` with direct assignment and image flag handling
3. **Region Detection**: Add `is_molecule_in_atomistic_region()` helper method
4. **Integration**: Extend `post_force()` to handle atomistic region molecules
5. **Periodic Boundaries**: Use `domain->remap()` after updating CG particle position

All research questions resolved. Ready for Phase 1 design.
