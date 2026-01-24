#!/usr/bin/env python3
"""
Analyze density profile from AdResS simulation.

This script reads the density profile data and creates plots to validate
that the density is uniform across atomistic, transition, and CG regions.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os

def read_density_file(filename):
    """Read LAMMPS density profile file.
    
    Format from fix ave/chunk:
    # Header comments (3 lines)
    # Timestep Number-of-chunks Total-count
    timestep nchunks total_count
    chunk_id coord1 ncount density/number
    ...
    """
    if not os.path.exists(filename):
        print(f"Warning: {filename} not found")
        return None, None
    
    # Read file and skip comment lines and timestep line
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Find first non-comment, non-empty line (should be timestep line)
    data_start = 0
    for i, line in enumerate(lines):
        if line.strip() and not line.strip().startswith('#'):
            data_start = i
            break
    
    # Skip timestep line (has 3 columns) and read chunk data (has 4 columns)
    # Format: chunk_id coord1 ncount density
    data_lines = []
    for line in lines[data_start + 1:]:  # Skip timestep line
        line = line.strip()
        if line and not line.startswith('#'):
            parts = line.split()
            if len(parts) >= 4:  # chunk_id coord1 ncount density
                data_lines.append([float(parts[1]), float(parts[3])])  # coord1, density
    
    if not data_lines:
        print(f"Warning: {filename} contains no data")
        return None, None
    
    data = np.array(data_lines)
    x = data[:, 0]  # coord1 (bin center)
    density = data[:, 1]  # density/number
    
    return x, density

def plot_density_profile(density_files, output_file='density_profile.png'):
    """Plot density profiles from LAMMPS output files."""
    
    fig, axes = plt.subplots(2, 1, figsize=(10, 8))
    
    # Top plot: Total density
    ax1 = axes[0]
    
    # Try to read total density
    if 'density_profile.dat' in density_files:
        x_total, rho_total = read_density_file('density_profile.dat')
        if x_total is not None:
            ax1.plot(x_total, rho_total, 'b-', linewidth=2, label='Total density')
    
    # Try to read atomistic and CG densities separately
    x_at = None
    rho_at = None
    x_cg = None
    rho_cg = None
    
    if 'density_atomistic.dat' in density_files:
        x_at, rho_at = read_density_file('density_atomistic.dat')
        if x_at is not None:
            ax1.plot(x_at, rho_at, 'r--', linewidth=1.5, label='Atomistic (types 1-2)')
    
    if 'density_cg.dat' in density_files:
        x_cg, rho_cg = read_density_file('density_cg.dat')
        if x_cg is not None:
            ax1.plot(x_cg, rho_cg, 'g--', linewidth=1.5, label='CG (type 3)')
    
    # Add region boundaries (LJ dimer system)
    box_x = 20.0
    at_center = box_x / 2.0
    at_width = 8.0
    trans_width = 2.0
    
    at_xlo = at_center - at_width / 2.0
    at_xhi = at_center + at_width / 2.0
    trans1_xlo = at_xlo - trans_width
    trans1_xhi = at_xlo
    trans2_xlo = at_xhi
    trans2_xhi = at_xhi + trans_width
    
    # Draw vertical lines for region boundaries
    ax1.axvline(trans1_xlo, color='gray', linestyle=':', alpha=0.7, label='Region boundaries')
    ax1.axvline(at_xlo, color='gray', linestyle=':', alpha=0.7)
    ax1.axvline(at_xhi, color='gray', linestyle=':', alpha=0.7)
    ax1.axvline(trans2_xhi, color='gray', linestyle=':', alpha=0.7)
    
    # Add text labels for regions
    ax1.text(trans1_xlo/2, ax1.get_ylim()[1]*0.9, 'CG', ha='center', fontsize=10)
    ax1.text((trans1_xlo + at_xlo)/2, ax1.get_ylim()[1]*0.9, 'Transition', ha='center', fontsize=10)
    ax1.text((at_xlo + at_xhi)/2, ax1.get_ylim()[1]*0.9, 'Atomistic', ha='center', fontsize=10)
    ax1.text((at_xhi + trans2_xhi)/2, ax1.get_ylim()[1]*0.9, 'Transition', ha='center', fontsize=10)
    ax1.text((trans2_xhi + box_x)/2, ax1.get_ylim()[1]*0.9, 'CG', ha='center', fontsize=10)
    
    ax1.set_xlabel('x (Angstrom)')
    ax1.set_ylabel('Number density (1/Angstrom³)')
    ax1.set_title('AdResS Density Profile Validation')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Bottom plot: Relative density (normalized by average)
    ax2 = axes[1]
    
    if x_total is not None and rho_total is not None:
        rho_avg = np.mean(rho_total)
        rho_rel = rho_total / rho_avg
        ax2.plot(x_total, rho_rel, 'b-', linewidth=2, label='Total (normalized)')
        ax2.axhline(1.0, color='k', linestyle='--', alpha=0.5, label='Uniform density')
    
    # Draw region boundaries again
    ax2.axvline(trans1_xlo, color='gray', linestyle=':', alpha=0.7)
    ax2.axvline(at_xlo, color='gray', linestyle=':', alpha=0.7)
    ax2.axvline(at_xhi, color='gray', linestyle=':', alpha=0.7)
    ax2.axvline(trans2_xhi, color='gray', linestyle=':', alpha=0.7)
    
    ax2.set_xlabel('x (Angstrom)')
    ax2.set_ylabel('Relative density (normalized)')
    ax2.set_title('Normalized Density Profile (should be ~1.0 everywhere)')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim([0.8, 1.2])  # Focus on small deviations
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    print(f"Density profile plot saved to {output_file}")
    
    # Print statistics
    if x_total is not None and rho_total is not None:
        print("\nDensity Statistics:")
        print(f"  Mean density: {np.mean(rho_total):.6f} 1/Angstrom³")
        print(f"  Std deviation: {np.std(rho_total):.6f} 1/Angstrom³")
        print(f"  Relative std: {np.std(rho_total)/np.mean(rho_total)*100:.2f}%")
        
        # Check uniformity
        rho_rel = rho_total / np.mean(rho_total)
        max_dev = np.max(np.abs(rho_rel - 1.0)) * 100
        print(f"  Max deviation from uniform: {max_dev:.2f}%")
        
        if max_dev < 5.0:
            print("  ✓ Density is uniform (good!)")
        else:
            print("  ⚠ Density variations detected - check thermodynamic force")

if __name__ == '__main__':
    density_files = ['density_profile.dat', 'density_atomistic.dat', 'density_cg.dat']
    plot_density_profile(density_files)
    
    # Also create a simple text summary
    if os.path.exists('density_profile.dat'):
        x, rho = read_density_file('density_profile.dat')
        if x is not None:
            print(f"\nDensity profile summary written to density_summary.txt")
            with open('density_summary.txt', 'w') as f:
                f.write("# AdResS Density Profile Summary\n")
                f.write("# x (Angstrom)  density (1/Angstrom³)  relative_density\n")
                rho_avg = np.mean(rho)
                for xi, rhoi in zip(x, rho):
                    f.write(f"{xi:.3f}  {rhoi:.6e}  {rhoi/rho_avg:.6f}\n")

