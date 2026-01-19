# Quickstart: Phase 3 Atomistic Region COM Tracking

**Date**: 2025-01-27  
**Feature**: Phase 3: Atomistic Region COM Tracking for AdResS

## Overview

Phase 3 implements center-of-mass (COM) tracking for molecules in the atomistic region. When molecules are in the atomistic region (lambda > 0.9), the CG particle position and velocity are updated from the calculated COM of the molecule's atomistic atoms.

## Key Concepts

### Atomistic Region

Molecules with `lambda > lambda_at_threshold` (default 0.9) are in the atomistic region. In this region:
- Atoms move freely under atomistic forces
- CG particle tracks the COM of the atoms
- COM is calculated from atom positions and velocities

### COM Calculation

Center-of-mass is calculated using standard formulas:
- Position: `x_com = Σ(m_i * x_i) / Σ(m_i)`
- Velocity: `v_com = Σ(m_i * v_i) / Σ(m_i)`

Where the sum is over all atoms in the molecule.

### CG Particle Update

In the atomistic region, the CG particle is updated from COM:
- `x_cg = x_com`
- `v_cg = v_com`

This happens in `post_force()` after atoms have moved.

## Usage

### Basic Setup

```lammps
# Define regions
region atomistic block 10.0 20.0 0.0 20.0 0.0 20.0
region transition block 5.0 10.0 0.0 20.0 0.0 20.0
region cg block 0.0 5.0 0.0 20.0 0.0 20.0

# Setup AdResS region management
fix 1 all adress/region atomistic transition cg x

# Setup COM tracking (includes both CG constraints and atomistic COM tracking)
fix 2 all adress/constraint 1

# Normal integration
fix 3 all nve
```

### How It Works

1. **Region Detection**: `fix adress/region` calculates lambda for each atom based on position
2. **Molecule Classification**: `fix adress/constraint` determines if molecule is in atomistic region (lambda > 0.9)
3. **COM Calculation**: For atomistic region molecules, COM is calculated from atom positions/velocities
4. **CG Particle Update**: CG particle position and velocity are set to match COM
5. **Periodic Boundaries**: Image flags are updated to handle periodic boundary crossings

### Verification

To verify COM tracking is working:

```python
# In verification script
for timestep in dump_data:
    for molecule in molecules:
        if molecule.region == "atomistic":
            com_manual = calculate_com_manually(molecule.atoms)
            com_cg = molecule.cg_particle.position
            assert abs(com_manual - com_cg) < 1e-6
```

## Example Test Case

See `examples/USER/adres/phase3_atomistic_com/in.dimer_atomistic` for a complete example with:
- Dimer molecule in atomistic region
- COM tracking verification
- Periodic boundary handling

## Troubleshooting

### CG Particle Not Tracking COM

- Check that molecule is in atomistic region: `lambda > 0.9`
- Verify `fix adress/constraint` is called after `fix adress/region`
- Check that CG particle exists for molecule

### Incorrect COM Calculation

- Verify MPI communication: Compare single-processor vs multi-processor results
- Check periodic boundaries: Ensure coordinates are unwrapped correctly
- Verify masses: Check that atom masses are correct

### Performance Issues

- Profile COM calculation overhead
- Check MPI communication efficiency
- Consider caching lambda values if needed

## Integration with Phase 2

Phase 3 works seamlessly with Phase 2 (CG region constraints):
- Molecules in CG region: Use Phase 2 constraints (atoms follow CG particle)
- Molecules in atomistic region: Use Phase 3 COM tracking (CG particle follows COM)
- Molecules in transition region: Handled in Phase 4 (interpolation)

Both mechanisms operate simultaneously for different molecules in the same simulation.
