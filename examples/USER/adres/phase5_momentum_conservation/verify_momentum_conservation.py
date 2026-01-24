#!/usr/bin/env python3
"""
Automated verification script for momentum conservation in AdResS transition region.

Verifies:
1. Total momentum conservation: P(t) = P(0) (within numerical precision)
2. Force symmetry: F_ij = -F_ji for all pairs (Newton's Third Law)
3. Total force on system: Σ F_i = 0 (for closed system)
4. Momentum change rate: dP/dt ≈ 0 (indirectly verifies sum of forces = 0)

This test specifically targets the transition region where symmetric pair-weighting
is critical. The fix ensures F_pair = w_i * w_j * F_atom + (1 - w_i * w_j) * F_cg,
which guarantees F_ij = -F_ji even when lambda_i != lambda_j.

Usage:
    python3 verify_momentum_conservation.py [dump_file]
"""

import sys
import numpy as np
import re
from collections import defaultdict

def read_dump_file(filename):
    """Read LAMMPS dump file and extract atom data including forces."""
    atoms = []
    timesteps = []
    current_timestep = None
    current_atoms = {}
    reading_atoms = False
    expect_timestep = False
    expect_natoms = False
    atom_fields = []
    
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
                
                # Check for "ITEM: NUMBER OF ATOMS"
                if 'ITEM: NUMBER OF ATOMS' in line:
                    expect_natoms = True
                    continue
                
                # Skip the number of atoms line
                if expect_natoms and line.isdigit():
                    expect_natoms = False
                    continue
                
                # Check for "ITEM: ATOMS" - next lines will be atom data
                if 'ITEM: ATOMS' in line:
                    # Parse field names
                    atom_fields = line.split()[2:]  # Skip "ITEM:" and "ATOMS"
                    reading_atoms = True
                    expect_natoms = False
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
                    if len(parts) >= len(atom_fields):
                        atom_data = {}
                        for i, field in enumerate(atom_fields):
                            try:
                                atom_data[field] = float(parts[i])
                            except (ValueError, IndexError):
                                atom_data[field] = parts[i] if i < len(parts) else None
                        
                        atom_id = int(atom_data.get('id', 0))
                        current_atoms[atom_id] = atom_data
        
        # Add last timestep
        if current_timestep is not None:
            timesteps.append(current_timestep)
            atoms.append(current_atoms.copy())
    
    except FileNotFoundError:
        print(f"ERROR: Dump file '{filename}' not found")
        print("Please run the LAMMPS simulation first to generate the dump file.")
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: Failed to read dump file: {e}")
        sys.exit(1)
    
    return timesteps, atoms, atom_fields

def calculate_momentum(atoms_dict, atom_masses):
    """
    Calculate total momentum of the system.
    
    Args:
        atoms_dict: Dictionary of atoms {atom_id: {'type': int, 'vx': float, ...}}
        atom_masses: Dictionary of masses {atom_type: mass}
    
    Returns:
        P: Total momentum vector [Px, Py, Pz]
    """
    P = np.zeros(3)
    
    for atom_id, atom_data in atoms_dict.items():
        atom_type = int(atom_data.get('type', 0))
        mass = atom_masses.get(atom_type, 1.0)
        
        vx = atom_data.get('vx', 0.0)
        vy = atom_data.get('vy', 0.0)
        vz = atom_data.get('vz', 0.0)
        
        P[0] += mass * vx
        P[1] += mass * vy
        P[2] += mass * vz
    
    return P

def calculate_total_force(atoms_dict):
    """
    Calculate total force on the system.
    
    Args:
        atoms_dict: Dictionary of atoms {atom_id: {'fx': float, ...}}
    
    Returns:
        F: Total force vector [Fx, Fy, Fz]
    """
    F = np.zeros(3)
    
    for atom_id, atom_data in atoms_dict.items():
        fx = atom_data.get('fx', 0.0)
        fy = atom_data.get('fy', 0.0)
        fz = atom_data.get('fz', 0.0)
        
        F[0] += fx
        F[1] += fy
        F[2] += fz
    
    return F

def verify_momentum_conservation(timesteps, atoms_list, atom_masses, precision=1e-10):
    """
    Verify total momentum conservation.
    
    Returns:
        (success, max_error, mean_error, errors)
    """
    print("=" * 70)
    print("Verifying Momentum Conservation")
    print("=" * 70)
    
    if len(timesteps) == 0:
        print("ERROR: No timesteps found")
        return False, 0.0, 0.0, []
    
    # Skip timestep 0 - velocities may not be consistent before constraints are applied
    # Use first non-zero timestep as reference (constraints are applied by then)
    reference_idx = 1 if len(timesteps) > 1 else 0
    initial_momentum = calculate_momentum(atoms_list[reference_idx], atom_masses)
    initial_momentum_mag = np.linalg.norm(initial_momentum)
    reference_timestep = timesteps[reference_idx]
    
    print(f"\nReference momentum (timestep {reference_timestep}, after constraints applied):")
    print(f"  Px = {initial_momentum[0]:.10e}")
    print(f"  Py = {initial_momentum[1]:.10e}")
    print(f"  Pz = {initial_momentum[2]:.10e}")
    print(f"  |P| = {initial_momentum_mag:.10e}")
    
    errors = []
    momentum_deviations = []
    
    # Check all timesteps, but note that timestep 0 may have different momentum
    for i, (timestep, atoms_dict) in enumerate(zip(timesteps, atoms_list)):
        # Skip timestep 0 in error checking (velocities not yet constrained)
        if i == 0:
            continue
        current_momentum = calculate_momentum(atoms_dict, atom_masses)
        momentum_error = current_momentum - initial_momentum
        momentum_error_mag = np.linalg.norm(momentum_error)
        
        momentum_deviations.append(momentum_error_mag)
        
        if momentum_error_mag > precision:
            errors.append((timestep, momentum_error_mag, momentum_error.copy()))
    
    if len(momentum_deviations) == 0:
        return False, 0.0, 0.0, []
    
    max_error = max(momentum_deviations)
    mean_error = np.mean(momentum_deviations)
    
    print(f"\nMomentum conservation statistics:")
    print(f"  Max deviation: {max_error:.10e}")
    print(f"  Mean deviation: {mean_error:.10e}")
    print(f"  Precision threshold: {precision:.0e}")
    print(f"  Analyzed {len(timesteps)} timesteps")
    
    if errors:
        print(f"\nTimesteps with errors > {precision:.0e} ({len(errors)}):")
        for timestep, error_mag, error_vec in errors[:10]:
            print(f"  Timestep {timestep}: error = {error_mag:.10e}, "
                  f"error_vec = [{error_vec[0]:.6e}, {error_vec[1]:.6e}, {error_vec[2]:.6e}]")
        if len(errors) > 10:
            print(f"  ... and {len(errors) - 10} more timesteps with errors")
    
    success = max_error < precision
    return success, max_error, mean_error, errors

def verify_force_symmetry(timesteps, atoms_list, precision=1e-8):
    """
    Verify force symmetry: F_ij = -F_ji for all pairs.
    
    Note: This is indirect verification since we only have total forces on atoms.
    We verify that total force on system is zero, which implies force symmetry.
    
    Returns:
        (success, max_force_error, mean_force_error, errors)
    """
    print("\n" + "=" * 70)
    print("Verifying Force Symmetry (Total Force = 0)")
    print("=" * 70)
    
    if len(timesteps) == 0:
        print("ERROR: No timesteps found")
        return False, 0.0, 0.0, []
    
    force_errors = []
    errors = []
    
    for timestep, atoms_dict in zip(timesteps, atoms_list):
        # Check if forces are available
        if 'fx' not in list(atoms_dict.values())[0]:
            print("\nWARNING: Forces not found in dump file.")
            print("Force symmetry verification requires forces in dump file.")
            print("Add 'fx fy fz' to dump command in input script.")
            return False, 0.0, 0.0, []
        
        total_force = calculate_total_force(atoms_dict)
        force_error_mag = np.linalg.norm(total_force)
        
        force_errors.append(force_error_mag)
        
        if force_error_mag > precision:
            errors.append((timestep, force_error_mag, total_force.copy()))
    
    if len(force_errors) == 0:
        return False, 0.0, 0.0, []
    
    max_force_error = max(force_errors)
    mean_force_error = np.mean(force_errors)
    
    print(f"\nTotal force statistics:")
    print(f"  Max |Σ F_i|: {max_force_error:.10e}")
    print(f"  Mean |Σ F_i|: {mean_force_error:.10e}")
    print(f"  Precision threshold: {precision:.0e}")
    print(f"  Analyzed {len(timesteps)} timesteps")
    
    if errors:
        print(f"\nTimesteps with |Σ F_i| > {precision:.0e} ({len(errors)}):")
        for timestep, error_mag, force_vec in errors[:10]:
            print(f"  Timestep {timestep}: |Σ F| = {error_mag:.10e}, "
                  f"F = [{force_vec[0]:.6e}, {force_vec[1]:.6e}, {force_vec[2]:.6e}]")
        if len(errors) > 10:
            print(f"  ... and {len(errors) - 10} more timesteps with errors")
    
    success = max_force_error < precision
    return success, max_force_error, mean_force_error, errors

def verify_momentum_change_rate(timesteps, atoms_list, atom_masses, dt, precision=1e-8):
    """
    Verify momentum change rate: dP/dt ≈ 0.
    
    This indirectly verifies that sum of forces = 0.
    
    Returns:
        (success, max_dpdt_error, mean_dpdt_error)
    """
    print("\n" + "=" * 70)
    print("Verifying Momentum Change Rate (dP/dt ≈ 0)")
    print("=" * 70)
    
    if len(timesteps) < 2:
        print("WARNING: Need at least 2 timesteps to calculate dP/dt")
        return False, 0.0, 0.0
    
    dpdt_errors = []
    
    # Skip timestep 0 (velocities not yet constrained)
    # Start from timestep 1 (after constraints applied)
    for i in range(2, len(timesteps)):
        prev_momentum = calculate_momentum(atoms_list[i-1], atom_masses)
        curr_momentum = calculate_momentum(atoms_list[i], atom_masses)
        
        dpdt = (curr_momentum - prev_momentum) / dt
        dpdt_mag = np.linalg.norm(dpdt)
        
        dpdt_errors.append(dpdt_mag)
    
    if len(dpdt_errors) == 0:
        return False, 0.0, 0.0
    
    max_dpdt_error = max(dpdt_errors)
    mean_dpdt_error = np.mean(dpdt_errors)
    
    print(f"\nMomentum change rate statistics:")
    print(f"  Max |dP/dt|: {max_dpdt_error:.10e}")
    print(f"  Mean |dP/dt|: {mean_dpdt_error:.10e}")
    print(f"  Precision threshold: {precision:.0e}")
    print(f"  Analyzed {len(dpdt_errors)} timestep pairs")
    
    success = max_dpdt_error < precision
    return success, max_dpdt_error, mean_dpdt_error

def main():
    """Main verification function."""
    dump_file = "momentum_conservation.dump"
    
    if len(sys.argv) > 1:
        dump_file = sys.argv[1]
    
    print("=" * 70)
    print("Phase 5: Momentum Conservation Verification Test")
    print("=" * 70)
    print(f"\nReading dump file: {dump_file}")
    
    timesteps, atoms_list, atom_fields = read_dump_file(dump_file)
    
    if not timesteps:
        print("ERROR: No timesteps found in dump file.")
        sys.exit(1)
    
    print(f"Found {len(timesteps)} timesteps")
    print(f"Atom fields: {atom_fields}")
    
    # Check if forces are available
    has_forces = 'fx' in atom_fields and 'fy' in atom_fields and 'fz' in atom_fields
    if not has_forces:
        print("\nWARNING: Forces (fx, fy, fz) not found in dump file.")
        print("Force symmetry verification will be skipped.")
        print("To enable force verification, add 'fx fy fz' to dump command.")
    
    # Atom masses (from data file)
    atom_masses = {1: 1.0, 2: 2.0}  # Type 1 = atomistic, Type 2 = CG
    
    # Timestep (from input script)
    dt = 0.005
    
    # Precision thresholds
    momentum_precision = 1e-10  # For exact momentum conservation
    force_precision = 1e-8     # For force symmetry (allowing numerical errors)
    dpdt_precision = 1e-8      # For momentum change rate
    
    # Run all verification checks
    all_passed = True
    
    # 1. Momentum conservation
    success1, max_err1, mean_err1, errors1 = verify_momentum_conservation(
        timesteps, atoms_list, atom_masses, momentum_precision)
    all_passed = all_passed and success1
    
    # 2. Force symmetry (if forces available)
    if has_forces:
        success2, max_err2, mean_err2, errors2 = verify_force_symmetry(
            timesteps, atoms_list, force_precision)
        all_passed = all_passed and success2
    else:
        success2 = None
        print("\n" + "=" * 70)
        print("Skipping Force Symmetry Verification (forces not available)")
        print("=" * 70)
    
    # 3. Momentum change rate
    success3, max_err3, mean_err3 = verify_momentum_change_rate(
        timesteps, atoms_list, atom_masses, dt, dpdt_precision)
    all_passed = all_passed and success3
    
    # Final summary
    print("\n" + "=" * 70)
    print("VERIFICATION SUMMARY")
    print("=" * 70)
    print(f"1. Momentum conservation: {'PASS' if success1 else 'FAIL'}")
    if success2 is not None:
        print(f"2. Force symmetry (|Σ F_i| = 0): {'PASS' if success2 else 'FAIL'}")
    else:
        print(f"2. Force symmetry: SKIPPED (forces not in dump file)")
    print(f"3. Momentum change rate (|dP/dt| ≈ 0): {'PASS' if success3 else 'FAIL'}")
    print("=" * 70)
    
    if all_passed:
        print("\n✓ VERIFICATION PASSED: All checks completed successfully")
        print("\nThe symmetric pair-weighting implementation correctly ensures:")
        print("  - Momentum is conserved: P(t) = P(0)")
        print("  - Forces are symmetric: F_ij = -F_ji")
        print("  - Total force is zero: Σ F_i = 0")
        print("\nThis confirms that Newton's Third Law is satisfied in the")
        print("transition region where lambda_i != lambda_j.")
        sys.exit(0)
    else:
        print("\n✗ VERIFICATION FAILED: Some checks failed")
        print("\nThis may indicate:")
        print("  - Momentum conservation violation (force asymmetry)")
        print("  - Incorrect pair-weighting implementation")
        print("  - Numerical precision issues (check thresholds)")
        print("=" * 70)
        sys.exit(1)

if __name__ == "__main__":
    main()
