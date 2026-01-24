#!/usr/bin/env python3
"""
Generate LAMMPS data file for dimer fluid system.

This script creates a LAMMPS data file with dimer molecules (2 atoms per molecule)
placed on a lattice with random orientations. The system is designed for AdResS
verification with an elongated box geometry.
"""

import numpy as np
import argparse
import sys
import math


def generate_lattice_positions(box_size, num_molecules):
    """
    Generate molecule center positions on a simple cubic lattice.
    
    Args:
        box_size: Tuple (lx, ly, lz) - box dimensions
        num_molecules: Number of molecules to place
        
    Returns:
        Array of shape (num_molecules, 3) with molecule center positions
    """
    lx, ly, lz = box_size
    
    # Calculate grid spacing
    # We want approximately num_molecules molecules
    # For a simple cubic lattice: N = nx * ny * nz
    nx = int(np.ceil((num_molecules * lx / (ly * lz)) ** (1/3)))
    ny = int(np.ceil((num_molecules * ly / (lx * lz)) ** (1/3)))
    nz = int(np.ceil((num_molecules / (nx * ny))))
    
    # Ensure we have enough grid points
    while nx * ny * nz < num_molecules:
        if lx >= ly and lx >= lz:
            nx += 1
        elif ly >= lz:
            ny += 1
        else:
            nz += 1
    
    # Generate grid spacing
    dx = lx / nx
    dy = ly / ny
    dz = lz / nz
    
    # Generate positions
    positions = []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                if len(positions) >= num_molecules:
                    break
                x = (i + 0.5) * dx
                y = (j + 0.5) * dy
                z = (k + 0.5) * dz
                positions.append([x, y, z])
            if len(positions) >= num_molecules:
                break
        if len(positions) >= num_molecules:
            break
    
    return np.array(positions[:num_molecules])


def generate_random_orientations(num_molecules, seed=None):
    """
    Generate random orientation vectors for molecules.
    
    Args:
        num_molecules: Number of molecules
        seed: Random seed (optional)
        
    Returns:
        Array of shape (num_molecules, 3) with normalized orientation vectors
    """
    if seed is not None:
        np.random.seed(seed)
    
    # Generate random unit vectors (uniform on sphere)
    orientations = np.random.randn(num_molecules, 3)
    norms = np.linalg.norm(orientations, axis=1, keepdims=True)
    orientations = orientations / norms
    
    return orientations


def generate_atom_positions(molecule_centers, orientations, bond_length):
    """
    Generate atom positions for dimer molecules.
    
    Args:
        molecule_centers: Array (num_molecules, 3) - molecule center positions
        orientations: Array (num_molecules, 3) - normalized orientation vectors
        bond_length: Bond length (equilibrium distance)
        
    Returns:
        Array of shape (num_molecules * 2, 3) with all atom positions
    """
    num_molecules = len(molecule_centers)
    atom_positions = np.zeros((num_molecules * 2, 3))
    
    for i in range(num_molecules):
        center = molecule_centers[i]
        orientation = orientations[i]
        
        # Place atoms along orientation vector
        offset = 0.5 * bond_length * orientation
        atom_positions[2*i] = center - offset
        atom_positions[2*i + 1] = center + offset
    
    return atom_positions


def write_lammps_data(filename, atom_positions, box_size, bond_length, mass, num_molecules):
    """
    Write LAMMPS data file.
    
    Args:
        filename: Output filename
        atom_positions: Array (num_atoms, 3) - atom positions
        box_size: Tuple (lx, ly, lz) - box dimensions
        bond_length: Bond length
        mass: Mass per atom
        num_molecules: Number of molecules
    """
    num_atoms = len(atom_positions)
    num_bonds = num_molecules
    lx, ly, lz = box_size
    
    with open(filename, 'w') as f:
        # Header
        f.write("LAMMPS data file - Dimer fluid system\n")
        f.write("\n")
        f.write(f"{num_atoms} atoms\n")
        f.write(f"{num_bonds} bonds\n")
        f.write("1 atom types\n")
        f.write("1 bond types\n")
        f.write("\n")
        
        # Box dimensions
        f.write(f"0.0 {lx} xlo xhi\n")
        f.write(f"0.0 {ly} ylo yhi\n")
        f.write(f"0.0 {lz} zlo zhi\n")
        f.write("\n")
        
        # Masses
        f.write("Masses\n")
        f.write("\n")
        f.write(f"1 {mass}\n")
        f.write("\n")
        
        # Atoms
        f.write("Atoms # molecular\n")
        f.write("\n")
        for i in range(num_atoms):
            molecule_id = i // 2 + 1
            atom_type = 1
            x, y, z = atom_positions[i]
            # Format: atom-ID molecule-ID atom-type x y z ix iy iz
            # Image flags are 0 0 0 for initial configuration
            f.write(f"{i+1} {molecule_id} {atom_type} {x:.6f} {y:.6f} {z:.6f} 0 0 0\n")
        
        f.write("\n")
        
        # Bonds
        f.write("Bonds\n")
        f.write("\n")
        for i in range(num_molecules):
            bond_type = 1
            atom1_id = 2 * i + 1
            atom2_id = 2 * i + 2
            f.write(f"{i+1} {bond_type} {atom1_id} {atom2_id}\n")


def main():
    parser = argparse.ArgumentParser(description='Generate LAMMPS data file for dimer fluid system')
    parser.add_argument('--box-size', type=float, nargs=3, default=[30.0, 10.0, 10.0],
                        metavar=('Lx', 'Ly', 'Lz'),
                        help='Box dimensions (default: 30 10 10)')
    parser.add_argument('--num-molecules', type=int, default=100,
                        help='Number of dimer molecules (default: 100)')
    parser.add_argument('--bond-length', type=float, default=0.97,
                        help='Equilibrium bond length (default: 0.97)')
    parser.add_argument('--mass', type=float, default=1.0,
                        help='Mass per atom (default: 1.0)')
    parser.add_argument('--output', type=str, default='system.data',
                        help='Output data file name (default: system.data)')
    parser.add_argument('--seed', type=int, default=None,
                        help='Random seed for orientations (optional)')
    
    args = parser.parse_args()
    
    # Validate inputs
    if args.num_molecules <= 0:
        print("Error: num-molecules must be positive", file=sys.stderr)
        sys.exit(1)
    
    if any(s <= 0 for s in args.box_size):
        print("Error: box dimensions must be positive", file=sys.stderr)
        sys.exit(1)
    
    if args.bond_length <= 0:
        print("Error: bond-length must be positive", file=sys.stderr)
        sys.exit(1)
    
    if args.mass <= 0:
        print("Error: mass must be positive", file=sys.stderr)
        sys.exit(1)
    
    try:
        # Generate system
        print(f"Generating dimer system with {args.num_molecules} molecules...")
        molecule_centers = generate_lattice_positions(args.box_size, args.num_molecules)
        orientations = generate_random_orientations(args.num_molecules, seed=args.seed)
        atom_positions = generate_atom_positions(molecule_centers, orientations, args.bond_length)
        
        # Write LAMMPS data file
        print(f"Writing LAMMPS data file: {args.output}")
        write_lammps_data(args.output, atom_positions, args.box_size, 
                         args.bond_length, args.mass, args.num_molecules)
        
        print(f"Success! Generated {args.output}")
        print(f"  Atoms: {len(atom_positions)}")
        print(f"  Bonds: {args.num_molecules}")
        print(f"  Box size: {args.box_size[0]} x {args.box_size[1]} x {args.box_size[2]}")
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
