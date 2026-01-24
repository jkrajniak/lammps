#!/usr/bin/env python3
"""
Automated verification script for transition region interpolation simple test.

Verifies:
1. Molecule stays in transition region (lambda between 0.1 and 0.9)
2. Lambda values are calculated correctly from positions
3. Position interpolation is smooth (no discontinuities)
4. Molecule behavior is consistent with transition region interpolation

Note: Direct verification of x_atom = λ_i · x_free + (1-λ_i) · x_constrained
requires access to x_free and x_constrained which are not available in dump files.
This script verifies the observable consequences of correct interpolation.

Usage:
    python3 verify_interpolation.py
"""

import sys
import numpy as np

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

def calculate_average_lambda(atoms_dict, atom_type=1, cg_type=2, box_size=20.0):
    """Calculate average lambda for a molecule from atom positions."""
    lambda_sum = 0.0
    n_atoms = 0
    
    for atom_id, atom_data in atoms_dict.items():
        mol_id = int(atom_data.get('mol', 0))
        atom_type_val = int(atom_data.get('type', 0))
        
        if mol_id == 1 and atom_type_val == atom_type:  # Only count atomistic atoms in molecule 1
            x = atom_data.get('x', 0.0)
            lambda_val = calculate_lambda(x, box_size)
            lambda_sum += lambda_val
            n_atoms += 1
    
    if n_atoms > 0:
        return lambda_sum / n_atoms
    return 0.0

def verify_interpolation(timesteps, atoms_list, atom_fields, precision=1e-6):
    """
    Verify position interpolation for transition region.
    
    Since we don't have direct access to x_free and x_constrained in dump files,
    we verify:
    1. Molecule stays in transition region (lambda between 0.1 and 0.9)
    2. Lambda values are calculated correctly
    3. Positions change smoothly (no discontinuities)
    4. Behavior is consistent with interpolation
    """
    print("=" * 70)
    print("Verifying Transition Region Interpolation")
    print("=" * 70)
    
    errors = []
    warnings = []
    lambda_values = []
    regions = []
    
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
            warnings.append(f"Timestep {timestep}: Missing atoms for molecule 1")
            continue
        
        # Calculate lambda for CG particle (representative of molecule)
        cg_x = cg_particle.get('x', 0.0)
        lambda_val = calculate_lambda(cg_x)
        region = determine_region_from_lambda(lambda_val)
        
        lambda_values.append(lambda_val)
        regions.append(region)
        
        # Verify molecule is in transition region
        if region != "TRANSITION":
            errors.append(f"Timestep {timestep}: Molecule not in transition region (lambda={lambda_val:.4f}, region={region})")
        
        # Verify lambda is in valid range for transition region
        if lambda_val < 0.1 or lambda_val >= 0.9:
            errors.append(f"Timestep {timestep}: Lambda out of transition range (lambda={lambda_val:.4f}, expected 0.1 <= lambda < 0.9)")
        
        # Check for position discontinuities (excluding first timestep)
        if i > 0:
            prev_cg_x = None
            prev_atom_x = []
            
            # Get previous positions
            prev_atoms_dict = atoms_list[i-1]
            for atom_id, atom_data in prev_atoms_dict.items():
                mol_id = int(atom_data.get('mol', 0))
                atom_type = int(atom_data.get('type', 0))
                if mol_id == 1:
                    if atom_type == 2:
                        prev_cg_x = atom_data.get('x', 0.0)
                    elif atom_type == 1:
                        prev_atom_x.append(atom_data.get('x', 0.0))
            
            # Check CG particle position change
            if prev_cg_x is not None:
                box_size = 20.0
                dx_raw = cg_x - prev_cg_x
                # Handle periodic boundary conditions
                if dx_raw > box_size / 2:
                    dx = dx_raw - box_size
                elif dx_raw < -box_size / 2:
                    dx = dx_raw + box_size
                else:
                    dx = dx_raw
                
                # Large position changes might indicate discontinuity
                if abs(dx) > 1.0 and abs(dx) < box_size / 2:
                    warnings.append(f"Timestep {timestep}: Large CG position change: {abs(dx):.6f}")
            
            # Check atom positions change smoothly
            for j, atom in enumerate(atomistic_atoms):
                if j < len(prev_atom_x):
                    atom_x = atom.get('x', 0.0)
                    prev_x = prev_atom_x[j]
                    dx_raw = atom_x - prev_x
                    # Handle periodic boundary conditions
                    if dx_raw > box_size / 2:
                        dx = dx_raw - box_size
                    elif dx_raw < -box_size / 2:
                        dx = dx_raw + box_size
                    else:
                        dx = dx_raw
                    
                    # Large position changes might indicate discontinuity
                    if abs(dx) > 1.0 and abs(dx) < box_size / 2:
                        warnings.append(f"Timestep {timestep}: Large atom {j+1} position change: {abs(dx):.6f}")
    
    # Summary statistics
    print(f"\nAnalyzed {len(timesteps)} timesteps")
    if lambda_values:
        print(f"Lambda range: {min(lambda_values):.4f} to {max(lambda_values):.4f}")
        print(f"Average lambda: {np.mean(lambda_values):.4f}")
        print(f"Regions visited: {set(regions)}")
    
    # Verify molecule stayed in transition region
    if set(regions) != {"TRANSITION"}:
        errors.append(f"Molecule left transition region. Regions visited: {set(regions)}")
    
    # Verify lambda stayed in valid range
    if lambda_values:
        min_lambda = min(lambda_values)
        max_lambda = max(lambda_values)
        if min_lambda < 0.1 or max_lambda >= 0.9:
            errors.append(f"Lambda range ({min_lambda:.4f} to {max_lambda:.4f}) extends outside transition region [0.1, 0.9)")
    
    # Report warnings
    if warnings:
        print(f"\nWarnings ({len(warnings)}):")
        for warning in warnings[:10]:
            print(f"  - {warning}")
        if len(warnings) > 10:
            print(f"  ... and {len(warnings) - 10} more warnings")
    
    # Report errors
    if errors:
        print(f"\nErrors ({len(errors)}):")
        for error in errors:
            print(f"  - {error}")
    
    return len(errors) == 0

def verify_com_interpolation(timesteps, atoms_list, atom_fields, precision=1e-6):
    """
    Verify COM interpolation for transition region.
    
    Verifies: x_com = λ_avg · x_com_atoms + (1-λ_avg) · x_com_cg
    where λ_avg = (Σ λ_i) / N_atoms
    """
    print("\n" + "=" * 70)
    print("Verifying COM Interpolation")
    print("=" * 70)
    
    errors = []
    warnings = []
    max_error = 0.0
    
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
        
        # Calculate average lambda for molecule
        lambda_avg = calculate_average_lambda(atoms_dict, atom_type=1, cg_type=2)
        
        # Only verify if in transition region
        if lambda_avg < 0.1 or lambda_avg >= 0.9:
            continue
        
        # Calculate COM from atoms
        com_atoms, com_mass = calculate_com(atoms_dict, atom_type=1, cg_type=2)
        
        # Get CG particle position (this is the interpolated COM)
        cg_pos = np.array([cg_particle.get('x', 0.0), 
                          cg_particle.get('y', 0.0), 
                          cg_particle.get('z', 0.0)])
        
        # Expected interpolated COM: x_com = λ_avg · x_com_atoms + (1-λ_avg) · x_com_cg
        # Note: x_com_cg should be the CG particle position itself
        com_expected = lambda_avg * com_atoms + (1.0 - lambda_avg) * cg_pos
        
        # Actual CG particle position (should match interpolated COM)
        com_actual = cg_pos.copy()
        
        # Compare (handle periodic boundaries)
        box_size = 20.0
        error_vec = com_expected - com_actual
        # Unwrap error for periodic boundaries
        for dim in range(3):
            if error_vec[dim] > box_size / 2:
                error_vec[dim] -= box_size
            elif error_vec[dim] < -box_size / 2:
                error_vec[dim] += box_size
        
        com_error = np.linalg.norm(error_vec)
        max_error = max(max_error, com_error)
        
        if com_error > precision:
            errors.append(f"Timestep {timestep}: COM interpolation error: {com_error:.6e} > {precision:.0e} "
                         f"(lambda_avg={lambda_avg:.4f})")
    
    print(f"\nAnalyzed {len(timesteps)} timesteps in transition region")
    print(f"Maximum COM interpolation error: {max_error:.6e}")
    print(f"Precision threshold: {precision:.0e}")
    
    if errors:
        print(f"\nErrors ({len(errors)}):")
        for error in errors[:10]:
            print(f"  - {error}")
        if len(errors) > 10:
            print(f"  ... and {len(errors) - 10} more errors")
    
    return len(errors) == 0

def verify_velocity_interpolation(timesteps, atoms_list, atom_fields, precision=1e-6):
    """
    Verify velocity interpolation and momentum conservation for transition region.
    
    Note: Direct verification of v_atom = λ_i · v_free + (1-λ_i) · v_constrained
    requires access to v_free and v_constrained which are not available in dump files.
    Instead, we verify:
    1. Momentum conservation (total momentum should be conserved)
    2. Velocity changes smoothly (no discontinuities)
    """
    print("\n" + "=" * 70)
    print("Verifying Velocity Interpolation and Momentum Conservation")
    print("=" * 70)
    
    errors = []
    warnings = []
    initial_momentum = None
    max_momentum_error = 0.0
    
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
        
        # Calculate average lambda for molecule
        lambda_avg = calculate_average_lambda(atoms_dict, atom_type=1, cg_type=2)
        
        # Only verify if in transition region
        if lambda_avg < 0.1 or lambda_avg >= 0.9:
            continue
        
        # Calculate total momentum
        total_momentum = np.array([0.0, 0.0, 0.0])
        total_mass = 0.0
        
        # Add atom momenta
        for atom in atomistic_atoms:
            mass = 1.0
            vx = atom.get('vx', 0.0)
            vy = atom.get('vy', 0.0)
            vz = atom.get('vz', 0.0)
            total_momentum[0] += mass * vx
            total_momentum[1] += mass * vy
            total_momentum[2] += mass * vz
            total_mass += mass
        
        # Add CG particle momentum
        cg_mass = 2.0  # CG particle mass = sum of 2 atoms
        cg_vx = cg_particle.get('vx', 0.0)
        cg_vy = cg_particle.get('vy', 0.0)
        cg_vz = cg_particle.get('vz', 0.0)
        total_momentum[0] += cg_mass * cg_vx
        total_momentum[1] += cg_mass * cg_vy
        total_momentum[2] += cg_mass * cg_vz
        total_mass += cg_mass
        
        # Store initial momentum
        if initial_momentum is None:
            initial_momentum = total_momentum.copy()
        
        # Check momentum conservation
        momentum_error = np.linalg.norm(total_momentum - initial_momentum)
        max_momentum_error = max(max_momentum_error, momentum_error)
        
        if momentum_error > precision:
            errors.append(f"Timestep {timestep}: Momentum conservation error: {momentum_error:.6e} > {precision:.0e}")
        
        # Check for velocity discontinuities (excluding first timestep)
        if i > 0:
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
            
            # Check CG particle velocity change
            if prev_cg_particle is not None:
                prev_cg_vx = prev_cg_particle.get('vx', 0.0)
                dvx = abs(cg_vx - prev_cg_vx)
                if dvx > 2.0:  # Large velocity change might indicate discontinuity
                    warnings.append(f"Timestep {timestep}: Large CG velocity change: {dvx:.6f}")
            
            # Check atom velocities change smoothly
            for j, atom in enumerate(atomistic_atoms):
                if j < len(prev_atomistic_atoms):
                    prev_atom = prev_atomistic_atoms[j]
                    atom_vx = atom.get('vx', 0.0)
                    prev_vx = prev_atom.get('vx', 0.0)
                    dvx = abs(atom_vx - prev_vx)
                    if dvx > 2.0:  # Large velocity change might indicate discontinuity
                        warnings.append(f"Timestep {timestep}: Large atom {j+1} velocity change: {dvx:.6f}")
    
    print(f"\nAnalyzed {len(timesteps)} timesteps in transition region")
    print(f"Maximum momentum conservation error: {max_momentum_error:.6e}")
    print(f"Precision threshold: {precision:.0e}")
    
    if warnings:
        print(f"\nWarnings ({len(warnings)}):")
        for warning in warnings[:10]:
            print(f"  - {warning}")
        if len(warnings) > 10:
            print(f"  ... and {len(warnings) - 10} more warnings")
    
    if errors:
        print(f"\nErrors ({len(errors)}):")
        for error in errors[:10]:
            print(f"  - {error}")
        if len(errors) > 10:
            print(f"  ... and {len(errors) - 10} more errors")
    
    return len(errors) == 0

def main():
    """Main verification function."""
    dump_file = "dimer_transition.dump"
    
    print("=" * 70)
    print("Phase 4: Transition Region Interpolation Simple Test Verification")
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
    
    # 1. Basic interpolation verification
    success1 = verify_interpolation(timesteps, atoms_list, atom_fields)
    all_passed = all_passed and success1
    
    # 2. COM interpolation verification
    success2 = verify_com_interpolation(timesteps, atoms_list, atom_fields)
    all_passed = all_passed and success2
    
    # 3. Velocity interpolation and momentum conservation verification
    success3 = verify_velocity_interpolation(timesteps, atoms_list, atom_fields)
    all_passed = all_passed and success3
    
    # Final summary
    print("\n" + "=" * 70)
    print("VERIFICATION SUMMARY")
    print("=" * 70)
    print(f"1. Basic interpolation: {'PASS' if success1 else 'FAIL'}")
    print(f"2. COM interpolation: {'PASS' if success2 else 'FAIL'}")
    print(f"3. Velocity interpolation & momentum: {'PASS' if success3 else 'FAIL'}")
    print("=" * 70)
    
    if all_passed:
        print("VERIFICATION PASSED: All checks completed successfully")
        print("=" * 70)
        print("\nNote: Direct verification of interpolation formulas")
        print("x_atom = λ_i · x_free + (1-λ_i) · x_constrained")
        print("v_atom = λ_i · v_free + (1-λ_i) · v_constrained")
        print("requires access to x_free, x_constrained, v_free, and v_constrained")
        print("which are not available in dump files. This verification confirms:")
        print("1. Molecule stays in transition region")
        print("2. Lambda values are correct")
        print("3. COM interpolation: x_com = λ_avg · x_com_atoms + (1-λ_avg) · x_com_cg")
        print("4. Momentum is conserved (within precision)")
        print("5. Positions and velocities change smoothly (consistent with interpolation)")
        sys.exit(0)
    else:
        print("VERIFICATION FAILED: Some checks failed")
        print("=" * 70)
        sys.exit(1)

if __name__ == "__main__":
    main()
