# API Contract: fix adress/constraint (Atomistic Region Extension)

**Date**: 2025-01-27  
**Feature**: Phase 3: Atomistic Region COM Tracking

## Overview

This document describes the internal API contract for the atomistic region COM tracking functionality added to `fix adress/constraint`.

## Methods

### `calculate_com_from_atoms(int mol_idx)`

**Purpose**: Calculate center-of-mass position and velocity from atom positions/velocities.

**Preconditions**:
- `mol_idx` is valid (0 <= mol_idx < nmolecules)
- Molecule has at least one atom
- Atom positions and velocities are current

**Postconditions**:
- `x_com[mol_idx]` contains calculated COM position
- `v_com[mol_idx]` contains calculated COM velocity
- `mass_com[mol_idx]` contains total molecule mass
- COM calculation is accurate within 1e-6 precision

**Algorithm**:
1. Find reference atom (CG particle or first atom)
2. For each atom in molecule:
   - Unwrap coordinates relative to reference using minimum image convention
   - Accumulate `m_i * x_i` and `m_i * v_i`
3. Aggregate across MPI processors using `MPI_Allreduce`
4. Divide by total mass to get COM

**Returns**: void

**Side Effects**: Updates `x_com`, `v_com`, `mass_com` arrays

---

### `update_cg_from_com(int mol_idx)`

**Purpose**: Update CG particle position and velocity from calculated COM.

**Preconditions**:
- `mol_idx` is valid (0 <= mol_idx < nmolecules)
- `x_com[mol_idx]` and `v_com[mol_idx]` are calculated
- CG particle exists for molecule (`cg_atom_index[mol_idx] >= 0`)
- Molecule is in atomistic region

**Postconditions**:
- CG particle position equals COM position: `x[cg_idx] = x_com[mol_idx]`
- CG particle velocity equals COM velocity: `v[cg_idx] = v_com[mol_idx]`
- Image flags are updated for periodic boundaries

**Algorithm**:
1. Get CG particle index: `cg_idx = cg_atom_index[mol_idx]`
2. Set position: `x[cg_idx] = x_com[mol_idx]`
3. Set velocity: `v[cg_idx] = v_com[mol_idx]`
4. Update image flags: `domain->remap(x[cg_idx], image[cg_idx])`

**Returns**: void

**Side Effects**: Updates CG particle position, velocity, and image flags

---

### `is_molecule_in_atomistic_region(int mol_idx)`

**Purpose**: Check if molecule is in atomistic region based on lambda value.

**Preconditions**:
- `mol_idx` is valid (0 <= mol_idx < nmolecules)
- `fix_region` is initialized and lambda values are available

**Postconditions**:
- Returns true if molecule is in atomistic region (lambda > lambda_at_threshold)
- Returns false otherwise

**Algorithm**:
1. Get representative atom (CG particle or first atom in molecule)
2. Get lambda value from `fix_region->get_lambda(atom_idx)`
3. Compare with threshold: `lambda > lambda_at_threshold`

**Returns**: bool (true if in atomistic region, false otherwise)

**Side Effects**: None

---

### `post_force(int vflag)` (Extended)

**Purpose**: Calculate COM and update CG particles for atomistic region molecules.

**Preconditions**:
- Forces have been computed
- Atoms have moved from forces
- `fix_region` is initialized

**Postconditions**:
- For atomistic region molecules: COM calculated and CG particle updated
- For CG region molecules: Existing Phase 2 constraints applied

**Algorithm**:
1. Loop through all molecules
2. For each molecule:
   - If in CG region: Apply Phase 2 constraints (existing logic)
   - If in atomistic region:
     - Calculate COM: `calculate_com_from_atoms(mol_idx)`
     - Update CG particle: `update_cg_from_com(mol_idx)`
   - If in transition region: Skip (Phase 4)

**Returns**: void

**Side Effects**: Updates CG particle positions/velocities for atomistic region molecules

---

## Data Structures

### COM Arrays

- `x_com[nmolecules][3]`: COM positions
- `v_com[nmolecules][3]`: COM velocities
- `mass_com[nmolecules]`: Total molecule masses

### Molecule Mapping

- `molecule_map[natoms]`: Maps atom index to molecule index
- `cg_atom_index[nmolecules]`: Maps molecule index to CG particle atom index

---

## Error Handling

### Invalid Molecule Index

- If `mol_idx < 0` or `mol_idx >= nmolecules`: Return early, no error (defensive)

### Missing CG Particle

- If `cg_atom_index[mol_idx] < 0`: Skip update, log warning

### Zero Mass

- If `mass_com[mol_idx] == 0`: Set COM to zero, log warning

### Missing Reference Atom

- If no reference atom found: Set COM to zero, return early

---

## Performance Contract

- COM calculation: O(N) where N is number of atoms in molecule
- MPI communication: O(1) per molecule (fixed-size Allreduce)
- Total overhead: < 5% of simulation time for 1000+ molecules

---

## Testing Contract

- Unit tests: Test each method in isolation
- Integration tests: Test with full LAMMPS simulation
- Verification: Compare with manual COM calculations (1e-6 precision)
