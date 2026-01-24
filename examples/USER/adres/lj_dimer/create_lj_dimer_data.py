#!/usr/bin/env python3
"""
Create initial configuration for LJ dimer AdResS simulation.

Atomistic: Two LJ particles connected by harmonic bond (dimer)
CG: Single LJ bead representing the dimer
"""

import numpy as np
import random

# System parameters
n_dimmers = 100  # Number of dimers
box_size = 20.0  # Box size in Angstrom (cubic)
bond_length = 2.0  # Equilibrium bond length in Angstrom

# LJ parameters (same for both atomistic and CG)
sigma = 3.4  # Angstrom (typical LJ sigma)
epsilon = 0.238  # kcal/mol (typical LJ epsilon)

# Atom types
# Type 1: First particle in dimer
# Type 2: Second particle in dimer
# Type 3: CG bead (represents the dimer)

print("Creating LJ dimer system...")
print(f"  Number of dimers: {n_dimmers}")
print(f"  Box size: {box_size} x {box_size} x {box_size} Angstrom")
print(f"  Bond length: {bond_length} Angstrom")

# Generate random positions for dimers
random.seed(12345)
np.random.seed(12345)

# Atomistic data file
atoms_atomistic = []
bonds = []
atom_id = 1
bond_id = 1

for i in range(n_dimmers):
    # Random center position for dimer
    x_center = random.uniform(0, box_size)
    y_center = random.uniform(0, box_size)
    z_center = random.uniform(0, box_size)
    
    # Random orientation vector (unit vector)
    theta = random.uniform(0, 2 * np.pi)
    phi = random.uniform(0, np.pi)
    dx = np.sin(phi) * np.cos(theta)
    dy = np.sin(phi) * np.sin(theta)
    dz = np.cos(phi)
    
    # Position of first particle
    x1 = x_center - 0.5 * bond_length * dx
    y1 = y_center - 0.5 * bond_length * dy
    z1 = z_center - 0.5 * bond_length * dz
    
    # Position of second particle
    x2 = x_center + 0.5 * bond_length * dx
    y2 = y_center + 0.5 * bond_length * dy
    z2 = z_center + 0.5 * bond_length * dz
    
    # Apply periodic boundary conditions
    coords = [x1, x2, y1, y2, z1, z2]
    for i in range(len(coords)):
        if coords[i] < 0:
            coords[i] += box_size
        elif coords[i] >= box_size:
            coords[i] -= box_size
    x1, x2, y1, y2, z1, z2 = coords
    
    # Add atoms
    atoms_atomistic.append((atom_id, 1, x1, y1, z1))  # Type 1
    atom_id += 1
    atoms_atomistic.append((atom_id, 2, x2, y2, z2))  # Type 2
    atom_id += 1
    
    # Add bond
    bonds.append((bond_id, atom_id - 2, atom_id - 1))
    bond_id += 1

# Write atomistic data file
print("\nWriting atomistic data file...")
with open('data.lj_dimer.atomistic', 'w') as f:
    f.write("# LJ dimer system - atomistic only\n\n")
    f.write(f"{len(atoms_atomistic)} atoms\n")
    f.write(f"{len(bonds)} bonds\n")
    f.write("0 angles\n")
    f.write("0 dihedrals\n")
    f.write("0 impropers\n\n")
    
    f.write("2 atom types\n")
    f.write("1 bond types\n\n")
    
    f.write(f"0.0 {box_size} xlo xhi\n")
    f.write(f"0.0 {box_size} ylo yhi\n")
    f.write(f"0.0 {box_size} zlo zhi\n\n")
    
    f.write("Masses\n\n")
    f.write("1 1.0  # First particle in dimer\n")
    f.write("2 1.0  # Second particle in dimer\n\n")
    
    f.write("Atoms # full\n\n")
    for atom_id, atom_type, x, y, z in atoms_atomistic:
        # LAMMPS full format: atom-ID molecule-ID atom-type q x y z
        molecule_id = (atom_id + 1) // 2  # Each dimer is one molecule
        f.write(f"{atom_id} {molecule_id} {atom_type} 0.0 {x:.6f} {y:.6f} {z:.6f}\n")
    
    f.write("\nBonds\n\n")
    for bond_id, atom1, atom2 in bonds:
        f.write(f"{bond_id} 1 {atom1} {atom2}\n")

print(f"  Created data.lj_dimer.atomistic with {len(atoms_atomistic)} atoms and {len(bonds)} bonds")

# Create AdResS data file (atomistic + CG)
print("\nCreating AdResS data file...")
atoms_adress = []
bonds_adress = []
atom_id = 1
bond_id = 1

for i in range(n_dimmers):
    # Get positions from atomistic configuration
    idx1 = i * 2
    idx2 = i * 2 + 1
    _, _, x1, y1, z1 = atoms_atomistic[idx1]
    _, _, x2, y2, z2 = atoms_atomistic[idx2]
    
    # Center of mass for CG bead
    x_cg = 0.5 * (x1 + x2)
    y_cg = 0.5 * (y1 + y2)
    z_cg = 0.5 * (z1 + z2)
    
    # Apply periodic boundary conditions
    coords_cg = [x_cg, y_cg, z_cg]
    for i in range(len(coords_cg)):
        if coords_cg[i] < 0:
            coords_cg[i] += box_size
        elif coords_cg[i] >= box_size:
            coords_cg[i] -= box_size
    x_cg, y_cg, z_cg = coords_cg
    
    # Add atomistic particles
    atoms_adress.append((atom_id, 1, x1, y1, z1))  # Type 1
    atom_id += 1
    atoms_adress.append((atom_id, 2, x2, y2, z2))  # Type 2
    atom_id += 1
    
    # Add CG bead
    atoms_adress.append((atom_id, 3, x_cg, y_cg, z_cg))  # Type 3 (CG)
    atom_id += 1
    
    # Add bond (only between atomistic particles)
    bonds_adress.append((bond_id, atom_id - 3, atom_id - 2))
    bond_id += 1

# Write AdResS data file
with open('data.lj_dimer.adress', 'w') as f:
    f.write("# LJ dimer system - AdResS (atomistic + CG)\n\n")
    f.write(f"{len(atoms_adress)} atoms\n")
    f.write(f"{len(bonds_adress)} bonds\n")
    f.write("0 angles\n")
    f.write("0 dihedrals\n")
    f.write("0 impropers\n\n")
    
    f.write("3 atom types\n")
    f.write("1 bond types\n\n")
    
    f.write(f"0.0 {box_size} xlo xhi\n")
    f.write(f"0.0 {box_size} ylo yhi\n")
    f.write(f"0.0 {box_size} zlo zhi\n\n")
    
    f.write("Masses\n\n")
    f.write("1 1.0  # First particle in dimer\n")
    f.write("2 1.0  # Second particle in dimer\n")
    f.write("3 2.0  # CG bead (mass = sum of dimer particles)\n\n")
    
    f.write("Atoms # full\n\n")
    for atom_id, atom_type, x, y, z in atoms_adress:
        # LAMMPS full format: atom-ID molecule-ID atom-type q x y z
        # Each dimer has 3 atoms: 2 atomistic + 1 CG, all in same molecule
        molecule_id = (atom_id + 2) // 3  # Each group of 3 atoms is one molecule
        f.write(f"{atom_id} {molecule_id} {atom_type} 0.0 {x:.6f} {y:.6f} {z:.6f}\n")
    
    f.write("\nBonds\n\n")
    for bond_id, atom1, atom2 in bonds_adress:
        f.write(f"{bond_id} 1 {atom1} {atom2}\n")

print(f"  Created data.lj_dimer.adress with {len(atoms_adress)} atoms ({2*n_dimmers} atomistic + {n_dimmers} CG) and {len(bonds_adress)} bonds")
print("\nDone!")

