#!/usr/bin/env python3
"""
Iterative verification manager for Force-based AdResS.

This script implements Iterative Boltzmann Inversion (IBI) to iteratively
tune the thermodynamic force until a flat density profile (target 0.844)
is achieved across all AdResS zones.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for cloud environments
import matplotlib.pyplot as plt
import subprocess
import argparse
import sys
import os
from pathlib import Path


def create_initial_force_table(filename, num_points=100, x_min=0.0, x_max=30.0):
    """
    Create initial thermodynamic force table with all zeros.
    
    Args:
        filename: Output filename
        num_points: Number of grid points
        x_min: Minimum x coordinate
        x_max: Maximum x coordinate
    """
    dx = (x_max - x_min) / num_points
    
    with open(filename, 'w') as f:
        f.write("# Thermodynamic force table\n")
        f.write(f"# {num_points} points from {x_min} to {x_max}\n")
        f.write(f"{num_points}\n")
        
        for i in range(num_points):
            x = x_min + (i + 0.5) * dx
            f.write(f"{x:.6f} 0.0\n")


def read_force_table(filename):
    """
    Read thermodynamic force table.
    
    Args:
        filename: Input filename
        
    Returns:
        x_values: Array of x coordinates
        force_values: Array of force values
    """
    x_values = []
    force_values = []
    
    with open(filename, 'r') as f:
        # Skip header comments
        line = f.readline()
        while line.strip().startswith('#'):
            line = f.readline()
        
        # Read number of points (if present)
        try:
            num_points = int(line.strip())
            line = f.readline()
        except ValueError:
            # No number of points line, just read data
            pass
        
        # Read data
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 2:
                x_values.append(float(parts[0]))
                force_values.append(float(parts[1]))
    
    return np.array(x_values), np.array(force_values)


def read_density_profile(filename, box_size=(30.0, 10.0, 10.0)):
    """
    Read density profile from LAMMPS output.
    
    LAMMPS fix ave/chunk output format:
    # Chunk Coord1 Ncount density/number
    1000 3000 400  <-- timestep line
    1 0.005 0 0
    2 0.015 0 0
    ...
    Where columns are: Chunk_ID Coord1 Ncount density/number
    
    Note: density/number is the COUNT of atoms per bin, not density.
    We need to convert to actual density by dividing by bin volume.
    
    Args:
        filename: Input filename (density_profile.dat)
        box_size: Tuple (lx, ly, lz) - box dimensions for density conversion
        
    Returns:
        x_values: Array of bin center coordinates
        density_values: Array of density values (atoms/volume)
    """
    x_values = []
    density_values = []
    
    if not os.path.exists(filename):
        raise FileNotFoundError(f"Density profile file not found: {filename}")
    
    lx, ly, lz = box_size
    
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Find the LAST timestep line (most recent data)
    last_timestep_idx = None
    for i, line in enumerate(lines):
        line = line.strip()
        if line.startswith('#'):
            continue
        parts = line.split()
        if len(parts) == 3:
            # This is a timestep line
            last_timestep_idx = i
    
    # Read data from the last timestep only
    start_idx = last_timestep_idx + 1 if last_timestep_idx is not None else 0
    
    # Count bins to calculate bin volume
    num_bins = 0
    for line in lines[start_idx:]:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        parts = line.split()
        if len(parts) >= 4:
            num_bins += 1
    
    # Calculate bin volume
    if num_bins > 0:
        bin_volume = (lx / num_bins) * ly * lz
    else:
        bin_volume = 1.0  # fallback
    
    # Read data from last timestep
    for line in lines[start_idx:]:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        parts = line.split()
        # Skip timestep lines
        if len(parts) == 3:
            continue
        # Read data lines: Chunk_ID Coord1 Ncount density/number
        if len(parts) >= 4:
            try:
                chunk_id = int(parts[0])
                x = float(parts[1])  # Coord1
                ncount = int(parts[2])  # Ncount
                # density/number is already density (atoms/volume), not count!
                density = float(parts[3])  # density/number (already density)
                
                x_values.append(x)
                density_values.append(density)
            except (ValueError, IndexError):
                continue
    
    if len(x_values) == 0:
        raise ValueError(f"No valid data found in density profile file: {filename}")
    
    return np.array(x_values), np.array(density_values)


def interpolate_to_grid(x_source, y_source, x_target):
    """
    Interpolate values from source grid to target grid.
    
    Args:
        x_source: Source x coordinates
        y_source: Source y values
        x_target: Target x coordinates
        
    Returns:
        y_target: Interpolated y values on target grid
    """
    return np.interp(x_target, x_source, y_source)


def update_force_table(input_file, output_file, density_x, density_values, 
                      target_density, prefactor, force_x):
    """
    Update thermodynamic force table using IBI formula.
    
    Args:
        input_file: Current force table file
        output_file: Output force table file
        density_x: Density profile x coordinates
        density_values: Density profile values
        target_density: Target density value
        prefactor: IBI update prefactor
        force_x: Force table x coordinates
    """
    # Read current force table
    force_x_current, force_values = read_force_table(input_file)
    
    # Interpolate density to force grid points
    density_interp = interpolate_to_grid(density_x, density_values, force_x_current)
    
    # Update forces using IBI formula: Force_new = Force_old - PreFactor * (Density - Target)
    force_new = force_values - prefactor * (density_interp - target_density)
    
    # Write updated table
    num_points = len(force_x_current)
    with open(output_file, 'w') as f:
        f.write("# Thermodynamic force table (updated)\n")
        f.write(f"# {num_points} points\n")
        f.write(f"{num_points}\n")
        
        for i in range(num_points):
            f.write(f"{force_x_current[i]:.6f} {force_new[i]:.6f}\n")


def plot_density_profile(x_values, density_values, target_density, iteration, output_dir):
    """
    Plot density profile and save to file.
    
    Args:
        x_values: Bin center coordinates
        density_values: Density values
        target_density: Target density (horizontal line)
        iteration: Iteration number
        output_dir: Output directory for plots
    """
    plt.figure(figsize=(10, 6))
    plt.plot(x_values, density_values, 'b-', linewidth=2, label='Measured Density')
    plt.axhline(y=target_density, color='r', linestyle='--', linewidth=2, label=f'Target ({target_density})')
    plt.xlabel('x position', fontsize=12)
    plt.ylabel('Density', fontsize=12)
    plt.title(f'Density Profile - Iteration {iteration}', fontsize=14)
    plt.legend(fontsize=10)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    # Save plot
    os.makedirs(output_dir, exist_ok=True)
    plot_filename = os.path.join(output_dir, f'iteration_{iteration}.png')
    plt.savefig(plot_filename, dpi=150)
    plt.close()
    
    print(f"  Plot saved: {plot_filename}")


def run_lammps(lammps_binary, input_script, data_file):
    """
    Run LAMMPS simulation.
    
    Args:
        lammps_binary: Path to LAMMPS binary
        input_script: LAMMPS input script
        data_file: LAMMPS data file
        
    Returns:
        exit_code: Exit code from LAMMPS
        stdout: Standard output
        stderr: Standard error
    """
    if not os.path.exists(lammps_binary):
        raise FileNotFoundError(f"LAMMPS binary not found: {lammps_binary}")
    
    if not os.path.exists(input_script):
        raise FileNotFoundError(f"Input script not found: {input_script}")
    
    if not os.path.exists(data_file):
        raise FileNotFoundError(f"Data file not found: {data_file}")
    
    # Run LAMMPS
    cmd = [lammps_binary, '-in', input_script]
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    return result.returncode, result.stdout, result.stderr


def main():
    parser = argparse.ArgumentParser(
        description='Iterative verification manager for Force-based AdResS'
    )
    parser.add_argument('--lammps-binary', type=str, default='./lmp_serial',
                        help='Path to LAMMPS binary (default: ./lmp_serial)')
    parser.add_argument('--input-script', type=str, default='in.adress_dimer',
                        help='LAMMPS input script (default: in.adress_dimer)')
    parser.add_argument('--data-file', type=str, default='system.data',
                        help='LAMMPS data file (default: system.data)')
    parser.add_argument('--max-iterations', type=int, default=20,
                        help='Maximum iterations (default: 20)')
    parser.add_argument('--prefactor', type=float, default=0.5,
                        help='IBI update prefactor (default: 0.5)')
    parser.add_argument('--target-density', type=float, default=0.844,
                        help='Target density (default: 0.844)')
    parser.add_argument('--table-file', type=str, default='thermo_force.table',
                        help='Thermodynamic force table file (default: thermo_force.table)')
    parser.add_argument('--plot-dir', type=str, default='plots',
                        help='Directory for plots (default: plots/)')
    parser.add_argument('--convergence-threshold', type=float, default=0.01,
                        help='Convergence threshold (default: 0.01)')
    
    args = parser.parse_args()
    
    # Validate inputs
    if args.prefactor <= 0 or args.prefactor > 1.0:
        print("Error: prefactor must be in (0, 1.0]", file=sys.stderr)
        sys.exit(1)
    
    if args.max_iterations <= 0:
        print("Error: max-iterations must be positive", file=sys.stderr)
        sys.exit(1)
    
    if args.target_density <= 0:
        print("Error: target-density must be positive", file=sys.stderr)
        sys.exit(1)
    
    # Step A: Initialization
    print("=" * 60)
    print("Force-based AdResS Verification")
    print("=" * 60)
    print(f"LAMMPS binary: {args.lammps_binary}")
    print(f"Input script: {args.input_script}")
    print(f"Data file: {args.data_file}")
    print(f"Max iterations: {args.max_iterations}")
    print(f"Prefactor: {args.prefactor}")
    print(f"Target density: {args.target_density}")
    print("=" * 60)
    
    # Create initial force table if it doesn't exist
    if not os.path.exists(args.table_file):
        print(f"\nCreating initial force table: {args.table_file}")
        create_initial_force_table(args.table_file)
    else:
        print(f"\nUsing existing force table: {args.table_file}")
    
    # Step B: Iterative Loop
    density_profile_file = 'density_profile.dat'
    
    for iteration in range(args.max_iterations):
        print(f"\n--- Iteration {iteration + 1}/{args.max_iterations} ---")
        
        # 1. Run LAMMPS simulation
        print("Running LAMMPS simulation...")
        exit_code, stdout, stderr = run_lammps(
            args.lammps_binary, args.input_script, args.data_file
        )
        
        if exit_code != 0:
            print(f"Error: LAMMPS simulation failed with exit code {exit_code}", file=sys.stderr)
            print("STDERR:", file=sys.stderr)
            print(stderr, file=sys.stderr)
            sys.exit(1)
        
        # 2. Analyze: Read density profile
        if not os.path.exists(density_profile_file):
            print(f"Error: Density profile file not found: {density_profile_file}", file=sys.stderr)
            sys.exit(1)
        
        print("Reading density profile...")
        try:
            # Read box size from data file or use defaults
            # For now, use default box size (30 x 10 x 10)
            density_x, density_values = read_density_profile(density_profile_file, box_size=(30.0, 10.0, 10.0))
        except Exception as e:
            print(f"Error reading density profile: {e}", file=sys.stderr)
            sys.exit(1)
        
        # 3. Calculate: Compare measured vs target
        density_error = density_values - args.target_density
        max_error = np.max(np.abs(density_error))
        mean_error = np.mean(np.abs(density_error))
        
        print(f"  Density range: [{np.min(density_values):.4f}, {np.max(density_values):.4f}]")
        print(f"  Max error: {max_error:.6f}")
        print(f"  Mean error: {mean_error:.6f}")
        
        # 4. Update: Apply IBI formula
        print("Updating thermodynamic force table...")
        force_x, _ = read_force_table(args.table_file)
        
        update_force_table(
            args.table_file, args.table_file,
            density_x, density_values,
            args.target_density, args.prefactor, force_x
        )
        
        # 5. Plot: Generate visualization
        print("Generating plot...")
        plot_density_profile(
            density_x, density_values,
            args.target_density, iteration, args.plot_dir
        )
        
        # Check convergence
        if max_error < args.convergence_threshold:
            print(f"\n✓ Converged! Max error ({max_error:.6f}) < threshold ({args.convergence_threshold})")
            break
    
    # Final summary
    print("\n" + "=" * 60)
    print("Verification complete")
    print("=" * 60)
    print(f"Final density profile plot: {args.plot_dir}/iteration_{iteration}.png")
    print(f"Final force table: {args.table_file}")


if __name__ == '__main__':
    main()
