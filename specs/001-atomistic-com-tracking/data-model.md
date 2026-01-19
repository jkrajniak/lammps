# Data Model: Phase 3 Atomistic Region COM Tracking

**Date**: 2025-01-27  
**Feature**: Phase 3: Atomistic Region COM Tracking for AdResS

## Entities

### Molecule

**Purpose**: Represents a collection of atoms with a corresponding CG particle in AdResS simulation.

**Attributes**:
- `mol_idx` (int): Internal molecule index (0 to nmolecules-1)
- `mol_id` (tagint): LAMMPS molecule ID (from atom->molecule)
- `cg_atom_index` (int): Local atom index of CG particle (-1 if none)
- `region` (enum): Current region classification (CG, Transition, Atomistic) - derived from lambda
- `x_com[3]` (double): Center-of-mass position (calculated)
- `v_com[3]` (double): Center-of-mass velocity (calculated)
- `mass_com` (double): Total mass of molecule (sum of atom masses)

**State Transitions**:
- Molecule can transition between regions based on lambda value:
  - CG region: `lambda < lambda_cg_threshold` (default 0.1)
  - Transition region: `lambda_cg_threshold <= lambda <= lambda_at_threshold`
  - Atomistic region: `lambda > lambda_at_threshold` (default 0.9)

**Validation Rules**:
- Molecule must have at least one atom
- CG particle index must be valid if molecule has CG particle
- Total mass must be positive (sum of atom masses > 0)
- COM position/velocity must be finite (not NaN or Inf)

**Relationships**:
- Has many: Atoms (via molecule_map)
- Has one: CG particle (via cg_atom_index, optional)

---

### Center of Mass (COM)

**Purpose**: Calculated position and velocity representing the collective motion of a molecule's atoms.

**Attributes**:
- `x_com[3]` (double): COM position in simulation box coordinates
- `v_com[3]` (double): COM velocity
- `mass_com` (double): Total mass used in calculation

**Calculation**:
- Position: `x_com = Σ(m_i * x_i) / Σ(m_i)` where sum is over all atoms in molecule
- Velocity: `v_com = Σ(m_i * v_i) / Σ(m_i)` where sum is over all atoms in molecule
- Mass: `mass_com = Σ(m_i)` where sum is over all atoms in molecule

**Validation Rules**:
- COM position must be within simulation box (or correctly unwrapped for periodic boundaries)
- COM velocity must be finite
- Mass must be positive
- Calculation must match manual calculation within 1e-6 precision

**Relationships**:
- Calculated from: Molecule's atoms
- Used to update: CG particle position and velocity

---

### CG Particle

**Purpose**: Dummy particle representing the molecule at coarse-grained resolution.

**Attributes**:
- `cg_idx` (int): Local atom index in LAMMPS atom arrays
- `x[3]` (double): Position (must track COM in atomistic region)
- `v[3]` (double): Velocity (must track COM velocity in atomistic region)
- `image` (imageint): Image flags for periodic boundaries
- `type` (int): Atom type (CG type, e.g., 4)
- `mass` (double): Mass (typically sum of atom masses)

**State in Atomistic Region**:
- Position: `x_cg = x_com` (updated from COM each timestep)
- Velocity: `v_cg = v_com` (updated from COM each timestep)
- Image flags: Updated via `domain->remap()` after position update

**Validation Rules**:
- CG particle must exist for molecule (cg_atom_index >= 0)
- Position must match COM within 1e-6 precision after update
- Velocity must match COM velocity within 1e-6 precision after update
- Image flags must be correct for periodic boundaries

**Relationships**:
- Belongs to: One molecule (via cg_atom_index)
- Updated from: COM (in atomistic region)

---

## Data Flow

### Atomistic Region COM Tracking Flow

1. **Input**: Molecule with atoms in atomistic region (lambda > threshold)
2. **Calculate COM**: 
   - Sum `m_i * x_i` and `m_i * v_i` for all atoms in molecule
   - Aggregate across MPI processors using `MPI_Allreduce`
   - Divide by total mass to get COM position and velocity
3. **Update CG Particle**:
   - Set `x[cg_idx] = x_com`
   - Set `v[cg_idx] = v_com`
   - Update image flags: `domain->remap(x[cg_idx], image[cg_idx])`
4. **Output**: CG particle position and velocity match COM

### Periodic Boundary Handling

1. **Reference Selection**: Choose reference atom (CG particle or first atom)
2. **Unwrap Coordinates**: For each atom, calculate relative position to reference using minimum image convention
3. **Calculate COM**: Use unwrapped coordinates in COM calculation
4. **Update CG Particle**: Set position and remap image flags

---

## Storage

### In-Memory Arrays (FixAdResSConstraint class)

- `x_com[nmolecules][3]`: COM positions for all molecules
- `v_com[nmolecules][3]`: COM velocities for all molecules
- `mass_com[nmolecules]`: Total masses for all molecules
- `cg_atom_index[nmolecules]`: CG particle indices for all molecules
- `molecule_map[natoms]`: Maps atom index to molecule index

### LAMMPS Atom Arrays (accessed via atom->)

- `atom->x[natoms][3]`: Atom positions (including CG particles)
- `atom->v[natoms][3]`: Atom velocities (including CG particles)
- `atom->image[natoms]`: Image flags for periodic boundaries
- `atom->molecule[natoms]`: Molecule IDs
- `atom->mass[ntypes]`: Per-type masses
- `atom->rmass[natoms]`: Per-atom masses (if used)

---

## Validation

### COM Calculation Validation

- Compare calculated COM with manual calculation: `|x_com_calculated - x_com_manual| < 1e-6`
- Verify MPI correctness: Single-processor and multi-processor results must match
- Check periodic boundaries: COM must be correct when molecules span boundaries

### CG Particle Update Validation

- Verify position match: `|x_cg - x_com| < 1e-6`
- Verify velocity match: `|v_cg - v_com| < 1e-6`
- Check image flags: CG particle image flags must be correct after remap
