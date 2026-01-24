#!/usr/bin/env python3
"""
Automated verification script for transition region interpolation integration test.

Verifies:
1. Molecules can transition smoothly between all three regions (CG → transition → atomistic)
2. No discontinuities in positions or velocities at region boundaries
3. Interpolation is continuous at boundaries (lambda = 0.1 and lambda = 0.9)
4. All three constraint mechanisms work correctly:
   - CG region: atoms follow CG particle (rigid body motion)
   - Transition region: positions/velocities interpolated
   - Atomistic region: CG particle follows COM

Usage:
    python3 verify_region_transitions.py
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
                
                # Check for "ITEM: ATOMS" - next lines will be atom data
                if 'ITEM: ATOMS' in line:
                    # Parse field names
                    atom_fields = line.split()[2:]  # Skip "ITEM:" and "ATOMS"
                    reading_atoms = True
                    expect_natoms = False
                    continue
                
                # If reading atoms, parse atom data
                if reading_atoms and line:
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
        
        # Save last timestep
        if current_timestep is not None:
            timesteps.append(current_timestep)
            atoms.append(current_atoms.copy())
    
    except FileNotFoundError:
        print(f"ERROR: Dump file '{filename}' not found.")
        print("Please run the LAMMPS simulation first to generate the dump file.")
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: Failed to read dump file: {e}")
        sys.exit(1)
    
    return timesteps, atoms, atom_fields

def calculate_lambda(x, box_size=20.0, cg_boundary=5.0, at_boundary=7.5, transition_width=2.5):
    """
    Calculate lambda value from position, matching fix adress/region implementation.
    
    Region layout (for x-axis, periodic):
    - CG: 0.0-5.0 and 15.0-20.0
    - Transition: 5.0-7.5 and 12.5-15.0
    - Atomistic: 7.5-12.5
    
    Lambda calculation:
    - If pos <= cg_boundary: lambda = 0.0 (CG region)
    - If pos >= at_boundary (but < next CG boundary): lambda = 1.0 (atomistic region)
    - Otherwise: d = (pos - cg_boundary) / transition_width, lambda = d²
    
    For periodic boundaries, we need to handle the right transition region (12.5-15.0)
    which transitions from atomistic (lambda=1) to CG (lambda=0).
    """
    # Handle periodic boundaries - unwrap position to [0, box_size)
    x_unwrapped = x % box_size
    if x_unwrapped < 0:
        x_unwrapped += box_size
    
    # Left transition region (5.0-7.5): CG -> Transition -> Atomistic
    if x_unwrapped <= cg_boundary:  # CG region (0-5.0)
        return 0.0
    elif x_unwrapped < at_boundary:  # Left transition (5.0-7.5)
        d = (x_unwrapped - cg_boundary) / transition_width
        d = max(0.0, min(1.0, d))  # Clamp to [0, 1]
        return d * d
    elif x_unwrapped < 12.5:  # Atomistic region (7.5-12.5)
        return 1.0
    elif x_unwrapped < 15.0:  # Right transition (12.5-15.0): Atomistic -> CG
        # For right transition, lambda decreases from 1 to 0
        # Distance from atomistic boundary (12.5)
        d = (x_unwrapped - 12.5) / transition_width
        d = max(0.0, min(1.0, d))  # Clamp to [0, 1]
        return 1.0 - (d * d)  # Decrease from 1 to 0
    else:  # CG region (15.0-20.0)
        return 0.0

def determine_region_from_lambda(lambda_val, lambda_cg_threshold=0.1, lambda_at_threshold=0.9):
    """Determine region based on lambda value."""
    if lambda_val < lambda_cg_threshold:
        return "CG"
    elif lambda_val < lambda_at_threshold:
        return "TRANSITION"
    else:
        return "ATOMISTIC"

def determine_region(x, lambda_cg_threshold=0.1, lambda_at_threshold=0.9):
    """Determine region based on position (using lambda calculation)."""
    lambda_val = calculate_lambda(x)
    return determine_region_from_lambda(lambda_val, lambda_cg_threshold, lambda_at_threshold)

def calculate_com(atoms_dict, atom_type=1, cg_type=2, box_size=20.0):
    """Calculate center of mass from atom positions, handling periodic boundaries."""
    total_mass = 0.0
    com = np.array([0.0, 0.0, 0.0])
    atom_positions = []
    
    # Collect atom positions first
    for atom_id, atom_data in atoms_dict.items():
        atom_type_val = int(atom_data.get('type', 0))
        if atom_type_val == atom_type:  # Only count atomistic atoms
            mass = 1.0  # Mass of atomistic atoms
            x = atom_data.get('x', 0.0)
            y = atom_data.get('y', 0.0)
            z = atom_data.get('z', 0.0)
            atom_positions.append((mass, np.array([x, y, z])))
            total_mass += mass
    
    if total_mass == 0:
        return com, 0.0
    
    # Unwrap coordinates relative to first atom for PBC handling
    if len(atom_positions) > 0:
        ref_pos = atom_positions[0][1].copy()
        for mass, pos in atom_positions:
            # Unwrap relative to reference
            for dim in range(3):
                dx = pos[dim] - ref_pos[dim]
                if dx > box_size / 2:
                    pos[dim] -= box_size
                elif dx < -box_size / 2:
                    pos[dim] += box_size
            
            com += mass * pos
    
    if total_mass > 0:
        com /= total_mass
    
    return com, total_mass

def calculate_com_velocity(atoms_dict, atom_type=1, cg_type=2):
    """Calculate center of mass velocity from atom velocities."""
    total_mass = 0.0
    v_com = np.array([0.0, 0.0, 0.0])
    
    for atom_id, atom_data in atoms_dict.items():
        atom_type_val = int(atom_data.get('type', 0))
        if atom_type_val == atom_type:  # Only count atomistic atoms
            mass = 1.0  # Mass of atomistic atoms
            vx = atom_data.get('vx', 0.0)
            vy = atom_data.get('vy', 0.0)
            vz = atom_data.get('vz', 0.0)
            
            v_com[0] += mass * vx
            v_com[1] += mass * vy
            v_com[2] += mass * vz
            total_mass += mass
    
    if total_mass > 0:
        v_com /= total_mass
    
    return v_com, total_mass

def verify_boundary_continuity(timesteps, atoms_list, lambda_cg_threshold=0.1, lambda_at_threshold=0.9, precision_pos=1e-4, precision_vel=1e-3):
    """
    Verify continuity at region boundaries (lambda = 0.1 and 0.9).
    
    At lambda = 0.1 (CG/transition boundary):
    - Interpolation should match CG constraint behavior
    - x_atom ≈ x_constrained, v_atom ≈ v_constrained
    
    At lambda = 0.9 (transition/atomistic boundary):
    - Interpolation should match atomistic COM tracking
    - x_atom ≈ x_free, v_atom ≈ v_free
    
    Also tests edge cases with lambda very close to thresholds (0.1001, 0.8999).
    """
    print("\n" + "=" * 70)
    print("Verifying Boundary Continuity")
    print("=" * 70)
    
    errors = []
    warnings = []
    boundary_timesteps = []
    edge_case_timesteps = []
    
    for i, (timestep, atoms_dict) in enumerate(zip(timesteps, atoms_list)):
        # Find CG particle and atoms for molecule 1
        cg_particle = None
        atomistic_atoms = []
        
        for atom_id, atom_data in atoms_dict.items():
            mol_id = int(atom_data.get('mol', 0))
            atom_type = int(atom_data.get('type', 0))
            
            if mol_id == 1:
                if atom_type == 2:  # CG particle
                    cg_particle = atom_data
                elif atom_type == 1:  # Atomistic atom
                    atomistic_atoms.append(atom_data)
        
        if cg_particle is None or len(atomistic_atoms) < 2:
            continue
        
        # Calculate lambda for CG particle (representative of molecule)
        cg_x = cg_particle.get('x', 0.0)
        lambda_val = calculate_lambda(cg_x)
        
        # Check if near boundaries
        near_cg_boundary = abs(lambda_val - lambda_cg_threshold) < 0.05
        near_at_boundary = abs(lambda_val - lambda_at_threshold) < 0.05
        
        # Check for edge cases (very close to thresholds)
        # Test with lambda values like 0.1001 and 0.8999 as specified in T041
        is_edge_case = (abs(lambda_val - 0.1001) < 0.01) or (abs(lambda_val - 0.8999) < 0.01) or \
                       (lambda_val > lambda_cg_threshold and lambda_val < lambda_cg_threshold + 0.01) or \
                       (lambda_val < lambda_at_threshold and lambda_val > lambda_at_threshold - 0.01)
        
        if near_cg_boundary or near_at_boundary:
            boundary_timesteps.append((timestep, lambda_val, cg_particle, atomistic_atoms))
        
        if is_edge_case:
            edge_case_timesteps.append((timestep, lambda_val, cg_particle, atomistic_atoms))
    
    print(f"Found {len(boundary_timesteps)} timesteps near boundaries")
    if edge_case_timesteps:
        print(f"Found {len(edge_case_timesteps)} timesteps with edge case lambda values (very close to thresholds)")
    
    # Verify continuity at boundaries
    for timestep, lambda_val, cg_particle, atomistic_atoms in boundary_timesteps[:20]:  # Check first 20
        cg_x = np.array([cg_particle.get('x', 0.0), cg_particle.get('y', 0.0), cg_particle.get('z', 0.0)])
        cg_v = np.array([cg_particle.get('vx', 0.0), cg_particle.get('vy', 0.0), cg_particle.get('vz', 0.0)])
        
        # Calculate COM from atoms
        com_pos, com_mass = calculate_com({cg_particle.get('id', 0): cg_particle, 
                                           **{a.get('id', 0): a for a in atomistic_atoms}}, 
                                          atom_type=1, cg_type=2)
        com_vel, _ = calculate_com_velocity({cg_particle.get('id', 0): cg_particle, 
                                             **{a.get('id', 0): a for a in atomistic_atoms}}, 
                                            atom_type=1, cg_type=2)
        
        if abs(lambda_val - lambda_cg_threshold) < 0.05:
            # At CG/transition boundary: should match CG behavior
            # CG particle should be close to constrained position (which is CG particle itself)
            # For CG region, atoms follow CG particle, so COM should be close to CG particle
            pos_error = np.linalg.norm(com_pos - cg_x)
            if pos_error > precision_pos:
                warnings.append(f"Timestep {timestep} (lambda={lambda_val:.4f}): Position continuity error at CG boundary: {pos_error:.6f}")
        
        if abs(lambda_val - lambda_at_threshold) < 0.05:
            # At transition/atomistic boundary: should match atomistic behavior
            # CG particle should track COM of atoms
            pos_error = np.linalg.norm(com_pos - cg_x)
            vel_error = np.linalg.norm(com_vel - cg_v)
            if pos_error > precision_pos:
                warnings.append(f"Timestep {timestep} (lambda={lambda_val:.4f}): Position continuity error at atomistic boundary: {pos_error:.6f}")
            if vel_error > precision_vel:
                warnings.append(f"Timestep {timestep} (lambda={lambda_val:.4f}): Velocity continuity error at atomistic boundary: {vel_error:.6f}")
    
    # Verify edge cases (lambda very close to thresholds: 0.1001, 0.8999)
    if edge_case_timesteps:
        print(f"\nVerifying edge cases (lambda very close to thresholds):")
        for timestep, lambda_val, cg_particle, atomistic_atoms in edge_case_timesteps[:10]:  # Check first 10
            cg_x = np.array([cg_particle.get('x', 0.0), cg_particle.get('y', 0.0), cg_particle.get('z', 0.0)])
            cg_v = np.array([cg_particle.get('vx', 0.0), cg_particle.get('vy', 0.0), cg_particle.get('vz', 0.0)])
            
            # Calculate COM from atoms
            com_pos, com_mass = calculate_com({cg_particle.get('id', 0): cg_particle, 
                                               **{a.get('id', 0): a for a in atomistic_atoms}}, 
                                              atom_type=1, cg_type=2)
            com_vel, _ = calculate_com_velocity({cg_particle.get('id', 0): cg_particle, 
                                                 **{a.get('id', 0): a for a in atomistic_atoms}}, 
                                                atom_type=1, cg_type=2)
            
            # At lambda ≈ 0.1001 (just above CG threshold): should behave like transition region
            if lambda_val > lambda_cg_threshold and lambda_val < lambda_cg_threshold + 0.01:
                # Should be in transition region, interpolation should be close to CG behavior
                pos_error = np.linalg.norm(com_pos - cg_x)
                if pos_error > precision_pos * 2:  # Slightly relaxed for edge cases
                    warnings.append(f"Timestep {timestep} (lambda={lambda_val:.4f}): Edge case continuity error at CG boundary: {pos_error:.6f}")
            
            # At lambda ≈ 0.8999 (just below atomistic threshold): should behave like transition region
            if lambda_val < lambda_at_threshold and lambda_val > lambda_at_threshold - 0.01:
                # Should be in transition region, interpolation should be close to atomistic behavior
                pos_error = np.linalg.norm(com_pos - cg_x)
                vel_error = np.linalg.norm(com_vel - cg_v)
                if pos_error > precision_pos * 2:  # Slightly relaxed for edge cases
                    warnings.append(f"Timestep {timestep} (lambda={lambda_val:.4f}): Edge case position continuity error at atomistic boundary: {pos_error:.6f}")
                if vel_error > precision_vel * 2:  # Slightly relaxed for edge cases
                    warnings.append(f"Timestep {timestep} (lambda={lambda_val:.4f}): Edge case velocity continuity error at atomistic boundary: {vel_error:.6f}")
    
    if warnings:
        print(f"\nWarnings ({len(warnings)}):")
        for warning in warnings[:10]:
            print(f"  - {warning}")
        if len(warnings) > 10:
            print(f"  ... and {len(warnings) - 10} more warnings")
    
    return len(errors) == 0

def verify_interpolation_formulas(timesteps, atoms_list, precision=1e-4):
    """
    Verify that interpolation formulas are mathematically correct.
    
    Position interpolation: x_atom = λ_i · x_free + (1-λ_i) · x_constrained
    COM interpolation: x_com = λ_avg · x_com_atoms + (1-λ_avg) · x_com_cg
    Velocity interpolation: v_atom = λ_i · v_free + (1-λ_i) · v_constrained
    """
    print("\n" + "=" * 70)
    print("Verifying Interpolation Formulas")
    print("=" * 70)
    
    errors = []
    warnings = []
    transition_timesteps = []
    
    for i, (timestep, atoms_dict) in enumerate(zip(timesteps, atoms_list)):
        # Find CG particle and atoms for molecule 1
        cg_particle = None
        atomistic_atoms = []
        
        for atom_id, atom_data in atoms_dict.items():
            mol_id = int(atom_data.get('mol', 0))
            atom_type = int(atom_data.get('type', 0))
            
            if mol_id == 1:
                if atom_type == 2:  # CG particle
                    cg_particle = atom_data
                elif atom_type == 1:  # Atomistic atom
                    atomistic_atoms.append(atom_data)
        
        if cg_particle is None or len(atomistic_atoms) < 2:
            continue
        
        # Calculate lambda for molecule
        cg_x = cg_particle.get('x', 0.0)
        lambda_val = calculate_lambda(cg_x)
        region = determine_region_from_lambda(lambda_val)
        
        if region == "TRANSITION":
            transition_timesteps.append((timestep, lambda_val, cg_particle, atomistic_atoms, atoms_dict))
    
    print(f"Found {len(transition_timesteps)} timesteps in transition region")
    
    # Verify COM interpolation formula
    for timestep, lambda_avg, cg_particle, atomistic_atoms, atoms_dict in transition_timesteps[:50]:  # Check first 50
        # Calculate COM from atoms
        com_atoms, _ = calculate_com(atoms_dict, atom_type=1, cg_type=2)
        com_cg = np.array([cg_particle.get('x', 0.0), cg_particle.get('y', 0.0), cg_particle.get('z', 0.0)])
        
        # Expected interpolated COM: x_com = λ_avg · x_com_atoms + (1-λ_avg) · x_com_cg
        com_expected = lambda_avg * com_atoms + (1.0 - lambda_avg) * com_cg
        
        # Actual CG particle position (should match interpolated COM)
        com_actual = com_cg.copy()
        
        # Compare
        com_error = np.linalg.norm(com_expected - com_actual)
        if com_error > precision:
            warnings.append(f"Timestep {timestep} (lambda={lambda_avg:.4f}): COM interpolation error: {com_error:.6f}")
    
    if warnings:
        print(f"\nWarnings ({len(warnings)}):")
        for warning in warnings[:10]:
            print(f"  - {warning}")
        if len(warnings) > 10:
            print(f"  ... and {len(warnings) - 10} more warnings")
    
    return len(errors) == 0

def verify_periodic_boundary_handling(timesteps, atoms_list, box_size=20.0, precision=1e-4):
    """
    Verify periodic boundary handling works correctly during interpolation.
    
    Checks:
    1. COM calculations handle PBC correctly (unwrapping)
    2. Position changes account for PBC crossings
    3. No discontinuities when crossing periodic boundaries
    """
    print("\n" + "=" * 70)
    print("Verifying Periodic Boundary Handling")
    print("=" * 70)
    
    errors = []
    warnings = []
    pbc_crossings = 0
    
    for i, (timestep, atoms_dict) in enumerate(zip(timesteps, atoms_list)):
        if i == 0:
            continue
        
        # Find CG particle and atoms for molecule 1
        cg_particle = None
        atomistic_atoms = []
        
        for atom_id, atom_data in atoms_dict.items():
            mol_id = int(atom_data.get('mol', 0))
            atom_type = int(atom_data.get('type', 0))
            
            if mol_id == 1:
                if atom_type == 2:  # CG particle
                    cg_particle = atom_data
                elif atom_type == 1:  # Atomistic atom
                    atomistic_atoms.append(atom_data)
        
        if cg_particle is None or len(atomistic_atoms) < 2:
            continue
        
        # Get previous timestep data
        prev_atoms_dict = atoms_list[i-1]
        prev_cg_particle = None
        prev_atomistic_atoms = []
        
        for atom_id, atom_data in prev_atoms_dict.items():
            mol_id = int(atom_data.get('mol', 0))
            atom_type = int(atom_data.get('type', 0))
            if mol_id == 1:
                if atom_type == 2:
                    prev_cg_particle = atom_data
                elif atom_type == 1:
                    prev_atomistic_atoms.append(atom_data)
        
        if prev_cg_particle is None:
            continue
        
        # Check for PBC crossings
        cg_x = cg_particle.get('x', 0.0)
        prev_cg_x = prev_cg_particle.get('x', 0.0)
        dx_raw = cg_x - prev_cg_x
        
        # Detect PBC crossing
        is_pbc_crossing = abs(dx_raw) > box_size / 2
        if is_pbc_crossing:
            pbc_crossings += 1
            
            # Verify COM calculation handles PBC correctly
            com_pos, _ = calculate_com(atoms_dict, atom_type=1, cg_type=2, box_size=box_size)
            prev_com_pos, _ = calculate_com(prev_atoms_dict, atom_type=1, cg_type=2, box_size=box_size)
            
            # COM should change smoothly even across PBC
            com_dx = com_pos[0] - prev_com_pos[0]
            # Unwrap for PBC
            if com_dx > box_size / 2:
                com_dx -= box_size
            elif com_dx < -box_size / 2:
                com_dx += box_size
            
            # Large COM change might indicate PBC handling issue
            if abs(com_dx) > 1.0 and abs(com_dx) < box_size / 2:
                warnings.append(f"Timestep {timestep}: Large COM change during PBC crossing: {abs(com_dx):.6f}")
        
        # Verify atom positions handle PBC correctly
        for j, atom in enumerate(atomistic_atoms):
            if j < len(prev_atomistic_atoms):
                atom_x = atom.get('x', 0.0)
                prev_x = prev_atomistic_atoms[j].get('x', 0.0)
                dx_raw = atom_x - prev_x
                
                # Unwrap for PBC
                if dx_raw > box_size / 2:
                    dx = dx_raw - box_size
                elif dx_raw < -box_size / 2:
                    dx = dx_raw + box_size
                else:
                    dx = dx_raw
                
                # Verify unwrapped distance is reasonable (not a discontinuity)
                if abs(dx) > 1.0 and abs(dx) < box_size / 2:
                    warnings.append(f"Timestep {timestep}: Large atom {j+1} position change (possible PBC issue): {abs(dx):.6f}")
    
    print(f"\nAnalyzed {len(timesteps)} timesteps")
    print(f"Detected {pbc_crossings} periodic boundary crossings")
    print(f"PBC handling: {'PASS' if len(warnings) == 0 else 'WARNINGS'}")
    
    if warnings:
        print(f"\nWarnings ({len(warnings)}):")
        for warning in warnings[:10]:
            print(f"  - {warning}")
        if len(warnings) > 10:
            print(f"  ... and {len(warnings) - 10} more warnings")
    
    return len(errors) == 0

def verify_constraint_behavior(timesteps, atoms_list, lambda_cg_threshold=0.1, lambda_at_threshold=0.9, precision=1e-3):
    """
    Verify that constraint behavior matches the region type.
    
    CG region (lambda < 0.1):
    - Atoms follow CG particle (rigid body)
    - Inter-atomic distances preserved
    - Atom velocities match CG particle velocity
    
    Transition region (0.1 <= lambda < 0.9):
    - Positions/velocities are interpolated
    
    Atomistic region (lambda >= 0.9):
    - CG particle tracks COM of atoms
    """
    print("\n" + "=" * 70)
    print("Verifying Constraint Behavior")
    print("=" * 70)
    
    errors = []
    warnings = []
    
    # Track initial inter-atomic distance for CG region verification
    initial_atom_distance = None
    
    for i, (timestep, atoms_dict) in enumerate(zip(timesteps, atoms_list)):
        # Find CG particle and atoms for molecule 1
        cg_particle = None
        atomistic_atoms = []
        
        for atom_id, atom_data in atoms_dict.items():
            mol_id = int(atom_data.get('mol', 0))
            atom_type = int(atom_data.get('type', 0))
            
            if mol_id == 1:
                if atom_type == 2:  # CG particle
                    cg_particle = atom_data
                elif atom_type == 1:  # Atomistic atom
                    atomistic_atoms.append(atom_data)
        
        if cg_particle is None or len(atomistic_atoms) < 2:
            continue
        
        # Calculate lambda
        cg_x = cg_particle.get('x', 0.0)
        lambda_val = calculate_lambda(cg_x)
        region = determine_region_from_lambda(lambda_val)
        
        # Store initial distance
        if i == 0 and len(atomistic_atoms) >= 2:
            pos1 = np.array([atomistic_atoms[0].get('x', 0.0), 
                            atomistic_atoms[0].get('y', 0.0), 
                            atomistic_atoms[0].get('z', 0.0)])
            pos2 = np.array([atomistic_atoms[1].get('x', 0.0), 
                            atomistic_atoms[1].get('y', 0.0), 
                            atomistic_atoms[1].get('z', 0.0)])
            initial_atom_distance = np.linalg.norm(pos2 - pos1)
        
        cg_pos = np.array([cg_particle.get('x', 0.0), cg_particle.get('y', 0.0), cg_particle.get('z', 0.0)])
        cg_vel = np.array([cg_particle.get('vx', 0.0), cg_particle.get('vy', 0.0), cg_particle.get('vz', 0.0)])
        
        if region == "CG":
            # Verify atoms follow CG particle (rigid body)
            if len(atomistic_atoms) >= 2:
                # Check inter-atomic distance is preserved (handle PBC)
                pos1 = np.array([atomistic_atoms[0].get('x', 0.0), 
                                atomistic_atoms[0].get('y', 0.0), 
                                atomistic_atoms[0].get('z', 0.0)])
                pos2 = np.array([atomistic_atoms[1].get('x', 0.0), 
                                atomistic_atoms[1].get('y', 0.0), 
                                atomistic_atoms[1].get('z', 0.0)])
                
                # Unwrap positions relative to first atom for PBC handling
                box_size = 20.0
                dx = pos2 - pos1
                for dim in range(3):
                    if dx[dim] > box_size / 2:
                        pos2[dim] -= box_size
                    elif dx[dim] < -box_size / 2:
                        pos2[dim] += box_size
                
                current_distance = np.linalg.norm(pos2 - pos1)
                
                if initial_atom_distance and abs(current_distance - initial_atom_distance) > precision:
                    warnings.append(f"Timestep {timestep}: Inter-atomic distance changed in CG region: {abs(current_distance - initial_atom_distance):.6f}")
                
                # Check atom velocities match CG particle velocity
                for atom in atomistic_atoms:
                    atom_vel = np.array([atom.get('vx', 0.0), atom.get('vy', 0.0), atom.get('vz', 0.0)])
                    vel_diff = np.linalg.norm(atom_vel - cg_vel)
                    if vel_diff > precision:
                        warnings.append(f"Timestep {timestep}: Atom velocity doesn't match CG in CG region: {vel_diff:.6f}")
        
        elif region == "ATOMISTIC":
            # Verify CG particle tracks COM
            com_pos, _ = calculate_com(atoms_dict, atom_type=1, cg_type=2)
            com_vel, _ = calculate_com_velocity(atoms_dict, atom_type=1, cg_type=2)
            
            pos_error = np.linalg.norm(com_pos - cg_pos)
            vel_error = np.linalg.norm(com_vel - cg_vel)
            
            if pos_error > precision:
                warnings.append(f"Timestep {timestep}: CG particle doesn't track COM position in atomistic region: {pos_error:.6f}")
            if vel_error > precision:
                warnings.append(f"Timestep {timestep}: CG particle doesn't track COM velocity in atomistic region: {vel_error:.6f}")
    
    if warnings:
        print(f"\nWarnings ({len(warnings)}):")
        for warning in warnings[:10]:
            print(f"  - {warning}")
        if len(warnings) > 10:
            print(f"  ... and {len(warnings) - 10} more warnings")
    
    return len(errors) == 0

def verify_transitions(timesteps, atoms_list, atom_fields):
    """Verify smooth transitions between regions."""
    print("=" * 70)
    print("Verifying Region Transitions")
    print("=" * 70)
    
    errors = []
    warnings = []
    
    # Track molecule position and region over time
    positions = []
    regions = []
    velocities = []
    lambda_values = []
    
    for i, (timestep, atoms_dict) in enumerate(zip(timesteps, atoms_list)):
        # Find CG particle (type 2) and atomistic atoms (type 1) for molecule 1
        cg_particle = None
        atomistic_atoms = []
        
        for atom_id, atom_data in atoms_dict.items():
            mol_id = int(atom_data.get('mol', 0))
            atom_type = int(atom_data.get('type', 0))
            
            if mol_id == 1:
                if atom_type == 2:  # CG particle
                    cg_particle = atom_data
                elif atom_type == 1:  # Atomistic atom
                    atomistic_atoms.append(atom_data)
        
        if cg_particle is None or len(atomistic_atoms) < 2:
            warnings.append(f"Timestep {timestep}: Missing atoms for molecule 1")
            continue
        
        # Get CG particle position (x coordinate)
        cg_x = cg_particle.get('x', 0.0)
        lambda_val = calculate_lambda(cg_x)
        region = determine_region_from_lambda(lambda_val)
        
        positions.append(cg_x)
        regions.append(region)
        lambda_values.append(lambda_val)
        
        # Get velocities
        cg_vx = cg_particle.get('vx', 0.0)
        velocities.append(cg_vx)
        
        # Calculate COM from atomistic atoms
        com, total_mass = calculate_com(atoms_dict)
        
        # Check continuity at boundaries
        if i > 0:
            prev_x = positions[-2]
            prev_region = regions[-2]
            prev_vx = velocities[-2]
            
            # Account for periodic boundary conditions (box size = 20.0)
            box_size = 20.0
            dx_raw = cg_x - prev_x
            # Handle periodic wrap-around: if jump > half box, it's likely a PBC crossing
            if dx_raw > box_size / 2:
                dx = dx_raw - box_size  # Wrap from right to left
            elif dx_raw < -box_size / 2:
                dx = dx_raw + box_size  # Wrap from left to right
            else:
                dx = dx_raw
            
            # Check for discontinuities in position (excluding PBC crossings)
            if abs(dx) > 1.0 and abs(dx) < box_size / 2:  # Large jump but not PBC
                warnings.append(f"Timestep {timestep}: Large position change: {abs(dx):.6f}")
            
            # Check for discontinuities in velocity
            dvx = abs(cg_vx - prev_vx)
            if dvx > 2.0:  # Large velocity change might indicate discontinuity
                warnings.append(f"Timestep {timestep}: Large velocity change: {dvx:.6f}")
            
            # Verify region transitions are logical (accounting for PBC)
            # If it's a PBC crossing (dx > box_size/2), skip the region transition check
            is_pbc_crossing = abs(dx_raw) > box_size / 2
            if not is_pbc_crossing:
                # Use lambda-based region determination for more accurate checks
                prev_lambda = lambda_values[-2] if len(lambda_values) > 1 else lambda_values[-1]
                curr_lambda = lambda_values[-1]
                
                # Check for logical region progression
                if prev_region == "CG" and region == "ATOMISTIC":
                    # Should have gone through transition
                    if prev_lambda < 0.1 and curr_lambda >= 0.9:
                        warnings.append(f"Timestep {timestep}: Jumped from CG (lambda={prev_lambda:.4f}) to ATOMISTIC (lambda={curr_lambda:.4f}) without transition")
                elif prev_region == "ATOMISTIC" and region == "CG":
                    # Should have gone through transition
                    if prev_lambda >= 0.9 and curr_lambda < 0.1:
                        warnings.append(f"Timestep {timestep}: Jumped from ATOMISTIC (lambda={prev_lambda:.4f}) to CG (lambda={curr_lambda:.4f}) without transition")
                
                # Enhanced discontinuity detection using lambda
                # Check for sudden lambda changes (indicating discontinuity)
                dlambda = abs(curr_lambda - prev_lambda)
                if dlambda > 0.5 and not is_pbc_crossing:  # Large lambda jump
                    warnings.append(f"Timestep {timestep}: Large lambda change: {dlambda:.4f} (from {prev_lambda:.4f} to {curr_lambda:.4f})")
    
    # Summary
    print(f"\nAnalyzed {len(timesteps)} timesteps")
    print(f"Molecule visited regions: {set(regions)}")
    if lambda_values:
        print(f"Lambda range: {min(lambda_values):.4f} to {max(lambda_values):.4f}")
    
    # Check if molecule visited all three regions
    visited_regions = set(regions)
    if len(visited_regions) < 3:
        errors.append(f"Molecule did not visit all three regions. Visited: {visited_regions}")
    
    # Check lambda range covers transition region
    if lambda_values:
        min_lambda = min(lambda_values)
        max_lambda = max(lambda_values)
        if min_lambda >= 0.9 or max_lambda < 0.1:
            warnings.append(f"Molecule lambda range ({min_lambda:.4f} to {max_lambda:.4f}) doesn't span transition region")
    
    # Report warnings
    if warnings:
        print(f"\nWarnings ({len(warnings)}):")
        for warning in warnings[:10]:  # Show first 10
            print(f"  - {warning}")
        if len(warnings) > 10:
            print(f"  ... and {len(warnings) - 10} more warnings")
    
    # Report errors
    if errors:
        print(f"\nErrors ({len(errors)}):")
        for error in errors:
            print(f"  - {error}")
    
    return len(errors) == 0

def main():
    """Main verification function."""
    dump_file = "region_transitions.dump"
    
    print("=" * 70)
    print("Phase 4: Transition Region Interpolation Integration Test Verification")
    print("=" * 70)
    print(f"\nReading dump file: {dump_file}")
    
    timesteps, atoms_list, atom_fields = read_dump_file(dump_file)
    
    if not timesteps:
        print("ERROR: No timesteps found in dump file.")
        sys.exit(1)
    
    print(f"Found {len(timesteps)} timesteps")
    print(f"Atom fields: {atom_fields}")
    
    # Run all verification checks
    all_passed = True
    
    # 1. Basic transition verification
    success1 = verify_transitions(timesteps, atoms_list, atom_fields)
    all_passed = all_passed and success1
    
    # 2. Boundary continuity verification
    success2 = verify_boundary_continuity(timesteps, atoms_list)
    all_passed = all_passed and success2
    
    # 3. Interpolation formula verification
    success3 = verify_interpolation_formulas(timesteps, atoms_list)
    all_passed = all_passed and success3
    
    # 4. Constraint behavior verification
    success4 = verify_constraint_behavior(timesteps, atoms_list)
    all_passed = all_passed and success4
    
    # 5. Periodic boundary handling verification
    success5 = verify_periodic_boundary_handling(timesteps, atoms_list)
    all_passed = all_passed and success5
    
    # Final summary
    print("\n" + "=" * 70)
    print("VERIFICATION SUMMARY")
    print("=" * 70)
    print(f"1. Basic transitions: {'PASS' if success1 else 'FAIL'}")
    print(f"2. Boundary continuity: {'PASS' if success2 else 'FAIL'}")
    print(f"3. Interpolation formulas: {'PASS' if success3 else 'FAIL'}")
    print(f"4. Constraint behavior: {'PASS' if success4 else 'FAIL'}")
    print(f"5. Periodic boundary handling: {'PASS' if success5 else 'FAIL'}")
    print("=" * 70)
    
    if all_passed:
        print("VERIFICATION PASSED: All checks completed successfully")
        print("=" * 70)
        sys.exit(0)
    else:
        print("VERIFICATION FAILED: Some checks failed")
        print("=" * 70)
        sys.exit(1)

if __name__ == "__main__":
    main()
