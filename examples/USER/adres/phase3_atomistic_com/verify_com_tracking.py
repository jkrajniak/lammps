#!/usr/bin/env python3
"""
Automated verification script for atomistic COM tracking test case.

Verifies:
1. CG particle position matches COM calculated from atoms (within 1e-6 precision)
2. CG particle velocity matches COM velocity calculated from atoms (within 1e-6 precision)
3. COM calculation is correct using standard formula: x_com = Σ(m_i * x_i) / Σ(m_i)

Usage:
    python3 verify_com_tracking.py
"""

import sys
import numpy as np
import re

def read_dump_file(filename):
    """Read LAMMPS dump file and extract atom data."""
    atoms = []
    timesteps = []
    current_timestep = None
    current_atoms = {}
    reading_atoms = False
    expect_timestep = False
    expect_natoms = False
    
    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                
                # Check for "ITEM: TIMESTEP" - next line will be timestep
                if 'ITEM: TIMESTEP' in line:
                    expect_timestep = True
                    # Save previous timestep if exists
                    if current_timestep is not None:
                        timesteps.append(current_timestep)
                        atoms.append(current_atoms.copy())
                    current_atoms = {}
                    reading_atoms = False
                    continue
                
                # If we expect a timestep, the next digit is the timestep
                if expect_timestep and line.isdigit():
                    current_timestep = int(line)
                    expect_timestep = False
                    expect_natoms = True
                    continue
                
                # Check for "ITEM: NUMBER OF ATOMS" - next line is number of atoms (skip it)
                if 'ITEM: NUMBER OF ATOMS' in line:
                    expect_natoms = True
                    continue
                
                # Skip the number of atoms line
                if expect_natoms and line.isdigit():
                    expect_natoms = False
                    continue
                
                # Check for ITEM: ATOMS header
                if 'ITEM: ATOMS' in line:
                    reading_atoms = True
                    continue
                
                # Skip other ITEM lines
                if 'ITEM:' in line:
                    reading_atoms = False
                    continue
                
                # Skip empty lines
                if not line:
                    continue
                
                # Parse atom data only if we're in the atoms section
                if reading_atoms:
                    parts = line.split()
                    if len(parts) >= 9:
                        atom_id = int(parts[0])
                        atom_type = int(parts[1])
                        mol_id = int(parts[2])
                        x, y, z = float(parts[3]), float(parts[4]), float(parts[5])
                        vx, vy, vz = float(parts[6]), float(parts[7]), float(parts[8])
                        
                        current_atoms[atom_id] = {
                            'type': atom_type,
                            'mol': mol_id,
                            'x': np.array([x, y, z]),
                            'v': np.array([vx, vy, vz])
                        }
            
            # Add last timestep
            if current_timestep is not None:
                timesteps.append(current_timestep)
                atoms.append(current_atoms.copy())
    
    except FileNotFoundError:
        print(f"ERROR: Dump file '{filename}' not found")
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: Failed to read dump file: {e}")
        sys.exit(1)
    
    return timesteps, atoms

def calculate_com(atoms_dict, atom_masses, cg_type=2):
    """
    Calculate center-of-mass position and velocity from atom positions.
    
    Args:
        atoms_dict: Dictionary of atoms {atom_id: {'type': int, 'x': array, 'v': array, ...}}
        atom_masses: Dictionary of masses {atom_type: mass}
        cg_type: Atom type of CG particles (excluded from COM calculation)
    
    Returns:
        (x_com, v_com, total_mass): COM position, COM velocity, total mass
    """
    sum_mx = np.zeros(3)
    sum_mv = np.zeros(3)
    sum_m = 0.0
    
    for atom_id, atom_data in atoms_dict.items():
        if atom_data['type'] == cg_type:
            continue  # Skip CG particles in COM calculation
        
        mass = atom_masses.get(atom_data['type'], 1.0)
        sum_mx += mass * atom_data['x']
        sum_mv += mass * atom_data['v']
        sum_m += mass
    
    if sum_m > 0.0:
        x_com = sum_mx / sum_m
        v_com = sum_mv / sum_m
    else:
        x_com = np.zeros(3)
        v_com = np.zeros(3)
    
    return x_com, v_com, sum_m

def find_cg_particle(atoms_dict, cg_type=2):
    """Find CG particle in atoms dictionary."""
    for atom_id, atom_data in atoms_dict.items():
        if atom_data['type'] == cg_type:
            return atom_id, atom_data
    return None, None

def main():
    """Main verification function."""
    dump_file = 'dimer_atomistic.dump'
    
    print("=" * 60)
    print("Phase 3: Atomistic COM Tracking Verification")
    print("=" * 60)
    print()
    
    # Read dump file
    print(f"Reading dump file: {dump_file}")
    timesteps, atoms_list = read_dump_file(dump_file)
    
    if len(timesteps) == 0:
        print("ERROR: No timesteps found in dump file")
        sys.exit(1)
    
    print(f"Found {len(timesteps)} timesteps")
    print()
    
    # Atom masses (from data file)
    atom_masses = {1: 1.0, 2: 2.0}  # Type 1 = atomistic, Type 2 = CG
    cg_type = 2
    
    # Precision threshold
    precision = 1e-6
    
    # Verification results
    position_errors = []
    velocity_errors = []
    failed_timesteps = []
    
    print("Verifying COM tracking...")
    print()
    
    # Skip timestep 0 (initial conditions)
    for idx, (timestep, atoms) in enumerate(zip(timesteps, atoms_list)):
        if timestep == 0:
            continue
        
        # Find CG particle
        cg_id, cg_data = find_cg_particle(atoms, cg_type)
        if cg_id is None:
            print(f"WARNING: CG particle not found at timestep {timestep}")
            continue
        
        # Calculate COM from atomistic atoms
        x_com_calc, v_com_calc, total_mass = calculate_com(atoms, atom_masses, cg_type)
        
        # Get CG particle position and velocity
        x_cg = cg_data['x']
        v_cg = cg_data['v']
        
        # Calculate errors
        pos_error = np.linalg.norm(x_com_calc - x_cg)
        vel_error = np.linalg.norm(v_com_calc - v_cg)
        
        position_errors.append(pos_error)
        velocity_errors.append(vel_error)
        
        # Check if errors exceed precision
        if pos_error > precision or vel_error > precision:
            failed_timesteps.append((timestep, pos_error, vel_error))
    
    # Report results
    print("=" * 60)
    print("VERIFICATION RESULTS")
    print("=" * 60)
    print()
    
    if len(position_errors) == 0:
        print("ERROR: No valid timesteps found for verification")
        sys.exit(1)
    
    max_pos_error = max(position_errors)
    max_vel_error = max(velocity_errors)
    mean_pos_error = np.mean(position_errors)
    mean_vel_error = np.mean(velocity_errors)
    
    print(f"Position tracking:")
    print(f"  Max error: {max_pos_error:.2e}")
    print(f"  Mean error: {mean_pos_error:.2e}")
    print(f"  Precision threshold: {precision:.2e}")
    print()
    
    print(f"Velocity tracking:")
    print(f"  Max error: {max_vel_error:.2e}")
    print(f"  Mean error: {mean_vel_error:.2e}")
    print(f"  Precision threshold: {precision:.2e}")
    print()
    
    # Final verdict
    all_pass = (max_pos_error < precision and max_vel_error < precision)
    
    if all_pass:
        print("✓ PASS: CG particle tracks COM within precision (1e-6)")
        print()
        print("Summary:")
        print(f"  - Position tracking: PASS (max error: {max_pos_error:.2e} < {precision:.2e})")
        print(f"  - Velocity tracking: PASS (max error: {max_vel_error:.2e} < {precision:.2e})")
        print(f"  - Verified {len(position_errors)} timesteps")
        return 0
    else:
        print("✗ FAIL: CG particle does not track COM within precision")
        print()
        if max_pos_error >= precision:
            print(f"  - Position tracking FAILED (max error: {max_pos_error:.2e} >= {precision:.2e})")
        if max_vel_error >= precision:
            print(f"  - Velocity tracking FAILED (max error: {max_vel_error:.2e} >= {precision:.2e})")
        print()
        if failed_timesteps:
            print(f"  Failed timesteps (showing first 5):")
            for timestep, pos_err, vel_err in failed_timesteps[:5]:
                print(f"    Timestep {timestep}: pos_error={pos_err:.2e}, vel_error={vel_err:.2e}")
        return 1

if __name__ == '__main__':
    sys.exit(main())
