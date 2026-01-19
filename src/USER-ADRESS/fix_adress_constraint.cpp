/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_adress_constraint.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix_adress_region.h"
#include "memory.h"
#include "modify.h"
#include "update.h"
#include "utils.h"

#include <cmath>
#include <cstring>
#include <map>
#include <set>
#include <vector>
#include <mpi.h>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixAdResSConstraint::FixAdResSConstraint(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), id_fix_region(nullptr), fix_region(nullptr), cg_type(0),
    molecule_map(nullptr), cg_atom_index(nullptr), mol_id_list(nullptr), nmolecules(0), maxmolecule(0),
    x_com(nullptr), v_com(nullptr), f_com(nullptr), x_cg_stored(nullptr),
    v_cg_stored(nullptr), mass_com(nullptr), displace(nullptr), maxatom(0),
    lambda_cg_threshold(0.1), lambda_at_threshold(0.9)
{
  if (narg < 4) utils::missing_cmd_args(FLERR, "fix adress/constraint", error);

  // Parse arguments: fix ID group-ID adress/constraint fix-region-ID cg-type
  id_fix_region = utils::strdup(arg[3]);
  cg_type = utils::numeric(FLERR, arg[4], false, lmp);

  if (cg_type < 1 || cg_type > atom->ntypes)
    error->all(FLERR, "Invalid CG atom type {} for fix adress/constraint", cg_type);

  // Parse optional thresholds
  int iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "lambda_cg") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix adress/constraint lambda_cg", error);
      lambda_cg_threshold = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if (lambda_cg_threshold < 0.0 || lambda_cg_threshold > 1.0)
        error->all(FLERR, "lambda_cg threshold must be between 0.0 and 1.0");
      iarg += 2;
    } else if (strcmp(arg[iarg], "lambda_at") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix adress/constraint lambda_at", error);
      lambda_at_threshold = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if (lambda_at_threshold < 0.0 || lambda_at_threshold > 1.0)
        error->all(FLERR, "lambda_at threshold must be between 0.0 and 1.0");
      iarg += 2;
    } else {
      error->all(FLERR, "Unknown fix adress/constraint keyword: {}", arg[iarg]);
    }
  }

  if (lambda_cg_threshold >= lambda_at_threshold)
    error->all(FLERR, "lambda_cg threshold ({}) must be less than lambda_at threshold ({})",
               lambda_cg_threshold, lambda_at_threshold);

  // Allocate arrays
  maxatom = atom->nmax;
  memory->create(molecule_map, maxatom, "adress/constraint:molecule_map");
  memory->create(displace, maxatom, 3, "adress/constraint:displace");

  // Initialize molecule_map to -1 (unmapped)
  for (int i = 0; i < maxatom; i++) {
    molecule_map[i] = -1;
    displace[i][0] = 0.0;
    displace[i][1] = 0.0;
    displace[i][2] = 0.0;
  }

  // Restart support
  restart_peratom = 0;
  restart_global = 0;
}

/* ---------------------------------------------------------------------- */

FixAdResSConstraint::~FixAdResSConstraint()
{
  delete[] id_fix_region;
  memory->destroy(molecule_map);
  memory->destroy(cg_atom_index);
  memory->destroy(mol_id_list);
  memory->destroy(x_com);
  memory->destroy(v_com);
  memory->destroy(f_com);
  memory->destroy(x_cg_stored);
  memory->destroy(v_cg_stored);
  memory->destroy(mass_com);
  memory->destroy(displace);
}

/* ---------------------------------------------------------------------- */

int FixAdResSConstraint::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= POST_FORCE;
  mask |= FINAL_INTEGRATE;
  mask |= POST_NEIGHBOR;  // To remap molecules if atoms migrate
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAdResSConstraint::init()
{
  // Get pointer to fix adress/region
  fix_region = dynamic_cast<FixAdResSRegion *>(modify->get_fix_by_id(id_fix_region));
  if (!fix_region)
    error->all(FLERR, "Fix {} for fix adress/constraint does not exist", id_fix_region);

  // Check that molecule IDs are available
  if (!atom->molecule)
    error->all(FLERR, "fix adress/constraint requires molecule IDs to be defined");

  // Map molecules
  map_molecules();
}

/* ---------------------------------------------------------------------- */

void FixAdResSConstraint::setup(int /*vflag*/)
{
  // Remap molecules in case system changed
  map_molecules();
}

/* ---------------------------------------------------------------------- */

void FixAdResSConstraint::post_neighbor()
{
  // Remap molecules after neighbor list rebuild (atoms may have migrated)
  map_molecules();
}

/* ---------------------------------------------------------------------- */

void FixAdResSConstraint::map_molecules()
{
  // This method maps atoms to molecules and finds corresponding CG particles
  // Called during init(), setup(), and post_neighbor()

  int nlocal = atom->nlocal;
  tagint *molecule = atom->molecule;
  int *type = atom->type;

  // First pass: count unique molecules and find CG particles
  std::map<tagint, int> mol_id_to_index;  // Map molecule ID to molecule index
  std::map<tagint, int> mol_id_to_cg;     // Map molecule ID to CG atom local index

  nmolecules = 0;

  // Find all unique molecule IDs and their CG particles
  for (int i = 0; i < nlocal; i++) {
    if (molecule[i] == 0) continue;  // Skip atoms without molecule ID

    // Check if this is a CG particle
    if (type[i] == cg_type) {
      if (mol_id_to_cg.find(molecule[i]) != mol_id_to_cg.end()) {
        error->warning(FLERR, "Multiple CG particles found for molecule ID {}. Using first one.",
                       molecule[i]);
      } else {
        mol_id_to_cg[molecule[i]] = i;
      }
    }

    // Track unique molecule IDs
    if (mol_id_to_index.find(molecule[i]) == mol_id_to_index.end()) {
      mol_id_to_index[molecule[i]] = nmolecules;
      nmolecules++;
    }
  }

  // MPI communication to get global molecule count and CG particle locations
  // Gather all molecule IDs from all processors to build global molecule list
  
  // Step 1: Collect local molecule IDs
  std::vector<tagint> local_mol_ids;
  for (auto &pair : mol_id_to_index) {
    local_mol_ids.push_back(pair.first);
  }
  int nlocal_mols = local_mol_ids.size();
  
  // Step 2: Gather counts from all processors
  int *mol_counts = nullptr;
  int *mol_displs = nullptr;
  if (comm->me == 0) {
    mol_counts = new int[comm->nprocs];
    mol_displs = new int[comm->nprocs];
  }
  MPI_Gather(&nlocal_mols, 1, MPI_INT, mol_counts, 1, MPI_INT, 0, world);
  
  // Step 3: Calculate displacements and total count
  int total_mols = 0;
  if (comm->me == 0) {
    mol_displs[0] = 0;
    for (int i = 0; i < comm->nprocs; i++) {
      if (i > 0) mol_displs[i] = mol_displs[i-1] + mol_counts[i-1];
      total_mols += mol_counts[i];
    }
  }
  
  // Step 4: Gather all molecule IDs
  std::vector<tagint> all_mol_ids;
  if (comm->me == 0) {
    all_mol_ids.resize(total_mols);
  }
  MPI_Gatherv(local_mol_ids.data(), nlocal_mols, MPI_LMP_TAGINT,
               comm->me == 0 ? all_mol_ids.data() : nullptr,
               mol_counts, mol_displs, MPI_LMP_TAGINT, 0, world);
  
  // Step 5: Build unique global molecule list (on proc 0)
  std::set<tagint> unique_mols;
  if (comm->me == 0) {
    for (tagint mol_id : all_mol_ids) {
      if (mol_id > 0) unique_mols.insert(mol_id);
    }
    nmolecules = unique_mols.size();
  }
  
  // Broadcast global molecule count
  MPI_Bcast(&nmolecules, 1, MPI_INT, 0, world);
  
  // Step 6: Broadcast unique molecule IDs to all processors
  std::vector<tagint> global_mol_list;
  if (comm->me == 0) {
    global_mol_list.assign(unique_mols.begin(), unique_mols.end());
  } else {
    global_mol_list.resize(nmolecules);
  }
  if (nmolecules > 0) {
    MPI_Bcast(global_mol_list.data(), nmolecules, MPI_LMP_TAGINT, 0, world);
  }
  
  // Step 7: Rebuild local molecule mapping using global list
  mol_id_to_index.clear();
  for (int m = 0; m < nmolecules; m++) {
    mol_id_to_index[global_mol_list[m]] = m;
  }
  
  // Step 8: Determine which processor owns each molecule's CG particle
  // For each molecule, find which processor has the CG particle locally
  // Initialize with large value, set to comm->me if we have CG particle, then use MPI_MIN
  std::vector<int> cg_owner(nmolecules, comm->nprocs);
  for (auto &pair : mol_id_to_cg) {
    tagint mol_id = pair.first;
    auto it = mol_id_to_index.find(mol_id);
    if (it != mol_id_to_index.end()) {
      cg_owner[it->second] = comm->me;  // This processor owns the CG particle
    }
  }
  
  // Reduce to find global owner (processor with lowest rank that has CG particle)
  // Use MPI_MIN to find lowest rank processor
  std::vector<int> global_cg_owner(nmolecules);
  MPI_Allreduce(cg_owner.data(), global_cg_owner.data(), nmolecules, 
                MPI_INT, MPI_MIN, world);
  
  // Clean up temporary arrays
  if (comm->me == 0) {
    delete[] mol_counts;
    delete[] mol_displs;
  }

  // Allocate molecule arrays if needed
  if (nmolecules > maxmolecule) {
    maxmolecule = nmolecules + 100;  // Add some padding
    memory->destroy(cg_atom_index);
    memory->destroy(mol_id_list);
    memory->destroy(x_com);
    memory->destroy(v_com);
    memory->destroy(f_com);
    memory->destroy(x_cg_stored);
    memory->destroy(v_cg_stored);
    memory->destroy(mass_com);

    memory->create(cg_atom_index, maxmolecule, "adress/constraint:cg_atom_index");
    memory->create(mol_id_list, maxmolecule, "adress/constraint:mol_id_list");
    memory->create(x_com, maxmolecule, 3, "adress/constraint:x_com");
    memory->create(v_com, maxmolecule, 3, "adress/constraint:v_com");
    memory->create(f_com, maxmolecule, 3, "adress/constraint:f_com");
    memory->create(x_cg_stored, maxmolecule, 3, "adress/constraint:x_cg_stored");
    memory->create(v_cg_stored, maxmolecule, 3, "adress/constraint:v_cg_stored");
    memory->create(mass_com, maxmolecule, "adress/constraint:mass_com");

    // Initialize
    for (int m = 0; m < maxmolecule; m++) {
      cg_atom_index[m] = -1;
      mol_id_list[m] = 0;
      x_com[m][0] = x_com[m][1] = x_com[m][2] = 0.0;
      v_com[m][0] = v_com[m][1] = v_com[m][2] = 0.0;
      f_com[m][0] = f_com[m][1] = f_com[m][2] = 0.0;
      x_cg_stored[m][0] = x_cg_stored[m][1] = x_cg_stored[m][2] = 0.0;
      v_cg_stored[m][0] = v_cg_stored[m][1] = v_cg_stored[m][2] = 0.0;
      mass_com[m] = 0.0;
    }
  }
  
  // Store molecule ID list
  for (int m = 0; m < nmolecules; m++) {
    if (m < (int)global_mol_list.size()) {
      mol_id_list[m] = global_mol_list[m];
    } else {
      mol_id_list[m] = 0;
    }
  }

  // Second pass: build molecule_map and cg_atom_index
  for (int i = 0; i < nlocal; i++) {
    if (molecule[i] == 0) {
      molecule_map[i] = -1;  // No molecule
      continue;
    }

    auto it = mol_id_to_index.find(molecule[i]);
    if (it != mol_id_to_index.end()) {
      molecule_map[i] = it->second;
    } else {
      molecule_map[i] = -1;
    }
  }

  // Build cg_atom_index
  // Only set cg_atom_index if this processor owns the CG particle
  for (auto &pair : mol_id_to_cg) {
    tagint mol_id = pair.first;
    int cg_local = pair.second;
    auto it = mol_id_to_index.find(mol_id);
    if (it != mol_id_to_index.end()) {
      int mol_idx = it->second;
      // Only set if this processor owns the CG particle
      if (global_cg_owner[mol_idx] == comm->me) {
        cg_atom_index[mol_idx] = cg_local;
      }
    }
  }

  // Warn about molecules without CG particles
  int nmol_no_cg = 0;
  for (int m = 0; m < nmolecules; m++) {
    if (cg_atom_index[m] == -1) {
      nmol_no_cg++;
    }
  }
  if (nmol_no_cg > 0 && comm->me == 0) {
    error->warning(FLERR, "fix adress/constraint: {} molecules have no CG particles. "
                   "These molecules will not be constrained.",
                   nmol_no_cg);
  }
}

/* ---------------------------------------------------------------------- */

int FixAdResSConstraint::find_cg_particle_for_molecule(tagint mol_id)
{
  // Find CG particle for given molecule ID
  // Returns local atom index, or -1 if not found
  int nlocal = atom->nlocal;
  tagint *molecule = atom->molecule;
  int *type = atom->type;

  for (int i = 0; i < nlocal; i++) {
    if (molecule[i] == mol_id && type[i] == cg_type) {
      return i;
    }
  }
  return -1;
}

/* ---------------------------------------------------------------------- */

void FixAdResSConstraint::initial_integrate(int /*vflag*/)
{
  // Called BEFORE normal integrator updates positions/velocities
  // For CG region molecules: constrain atoms to follow COM before integration
  // The normal integrator will update both CG particle and atoms, but we'll
  // re-constrain atoms in post_force() after forces are computed
  
  if (!fix_region) return;
  
  // Process each molecule
  for (int mol_idx = 0; mol_idx < nmolecules; mol_idx++) {
    // Check if molecule is in CG region
    if (!is_molecule_in_cg_region(mol_idx)) continue;
    
    // Get COM from CG particle (current position, will be updated by integrator)
    calculate_com_from_cg(mol_idx);
    
    // Constrain atom positions and velocities to follow COM
    // This ensures atoms start from constrained positions before integrator runs
    constrain_atoms_to_com(mol_idx);
  }
}

/* ---------------------------------------------------------------------- */

void FixAdResSConstraint::post_force(int /*vflag*/)
{
  // Called after force computation, before final_integrate()
  // For CG region molecules: 
  // 1. Calculate COM from CG particle (which has been updated by integrator)
  // 2. Store relative positions if needed
  // 3. Re-constrain atoms to follow updated COM (after integrator has moved them)
  // 4. Zero out forces on constrained atoms
  // For atomistic region molecules (Phase 3):
  // 1. Calculate COM from atom positions (atoms have moved from atomistic forces)
  // 2. Update CG particle position and velocity from COM
  // Transition region molecules are handled in Phase 4
  
  if (!fix_region) return;
  
  double **x = atom->x;
  double **f = atom->f;
  tagint *molecule = atom->molecule;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  
  // Process each molecule
  for (int mol_idx = 0; mol_idx < nmolecules; mol_idx++) {
    // Handle CG region molecules (Phase 2)
    if (is_molecule_in_cg_region(mol_idx)) {
      // Calculate COM from CG particle (now updated by integrator)
      calculate_com_from_cg(mol_idx);
      
      // Get COM position
      double *x_com_mol = x_com[mol_idx];
      
      // Check if relative positions need to be stored
      // Store if displace is effectively zero (hasn't been set yet)
      bool need_store = false;
      
      // Check if any atom in this molecule has zero displace (indicating not stored yet)
      for (int i = 0; i < nlocal; i++) {
        if (molecule_map[i] == mol_idx && type[i] != cg_type) {
          double disp_mag = displace[i][0]*displace[i][0] + 
                            displace[i][1]*displace[i][1] + 
                            displace[i][2]*displace[i][2];
          if (disp_mag < 1.0e-10) {
            need_store = true;
            break;
          }
        }
      }
      
      // Store relative positions if needed
      if (need_store) {
        for (int i = 0; i < nlocal; i++) {
          if (molecule_map[i] != mol_idx) continue;
          if (type[i] == cg_type) continue;  // Skip CG particles
          
          // Calculate relative position: displace = x_atom - x_com
          // Handle periodic boundaries using minimum image convention
          double dx[3];
          dx[0] = x[i][0] - x_com_mol[0];
          dx[1] = x[i][1] - x_com_mol[1];
          dx[2] = x[i][2] - x_com_mol[2];
          
          // Apply minimum image convention to get shortest distance
          domain->minimum_image(FLERR, dx[0], dx[1], dx[2]);
          
          // Store relative position
          displace[i][0] = dx[0];
          displace[i][1] = dx[1];
          displace[i][2] = dx[2];
        }
      }
      
      // Re-constrain atoms to follow updated COM (after integrator has moved CG particle)
      constrain_atoms_to_com(mol_idx, false);
      
      // Zero out forces on constrained atoms to prevent energy non-conservation
      // Constrained atoms move as rigid bodies with the CG particle
      // Forces on individual atoms would do work and break energy conservation
      for (int i = 0; i < nlocal; i++) {
        if (molecule_map[i] == mol_idx && type[i] != cg_type) {
          f[i][0] = 0.0;
          f[i][1] = 0.0;
          f[i][2] = 0.0;
        }
      }
    }
    // Handle atomistic region molecules (Phase 3)
    else if (is_molecule_in_atomistic_region(mol_idx)) {
      // Calculate COM from atom positions (atoms have moved from atomistic forces)
      calculate_com_from_atoms(mol_idx);
      
      // Update CG particle position and velocity from COM
      update_cg_from_com(mol_idx);
    }
    // Transition region molecules are handled in Phase 4
  }
}

/* ---------------------------------------------------------------------- */

void FixAdResSConstraint::final_integrate()
{
  // Called AFTER normal integrator updates velocities
  // For CG region molecules: constrain atom velocities to match COM velocity
  
  if (!fix_region) return;
  
  // Process each molecule
  for (int mol_idx = 0; mol_idx < nmolecules; mol_idx++) {
    // Check if molecule is in CG region
    if (!is_molecule_in_cg_region(mol_idx)) continue;
    
    // Get COM velocity from CG particle (updated by integrator)
    calculate_com_from_cg(mol_idx);
    
    // Constrain atom velocities to match COM velocity
    double **v = atom->v;
    int *type = atom->type;
    int nlocal = atom->nlocal;
    double *v_com_mol = v_com[mol_idx];
    
    // Update velocities for all atoms in this molecule
    for (int i = 0; i < nlocal; i++) {
      if (molecule_map[i] != mol_idx) continue;
      if (type[i] == cg_type) continue;  // Skip CG particles
      
      // Constrain atom velocity: v_atom = v_com
      v[i][0] = v_com_mol[0];
      v[i][1] = v_com_mol[1];
      v[i][2] = v_com_mol[2];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixAdResSConstraint::calculate_com_from_atoms(int mol_idx)
{
  // Calculate COM from atom positions for given molecule
  // Uses MPI to aggregate across processors for distributed molecules
  
  if (mol_idx < 0 || mol_idx >= nmolecules) return;
  
  tagint mol_id = mol_id_list[mol_idx];
  if (mol_id == 0) return;  // Invalid molecule ID
  
  double **x = atom->x;
  double **v = atom->v;
  imageint *image = atom->image;
  tagint *molecule = atom->molecule;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  
  // Partial sums for COM calculation
  double sum_mx[3] = {0.0, 0.0, 0.0};  // Σ(m_i * x_i)
  double sum_mv[3] = {0.0, 0.0, 0.0};  // Σ(m_i * v_i)
  double sum_m = 0.0;                   // Σ(m_i)
  
  // Unwrap coordinates relative to first atom in molecule (or CG particle if available)
  // Find reference atom (first atom or CG particle)
  double x_ref[3] = {0.0, 0.0, 0.0};
  int found_ref = 0;
  
  // Try to use CG particle as reference if available
  if (cg_atom_index[mol_idx] >= 0 && cg_atom_index[mol_idx] < nall) {
    int cg_idx = cg_atom_index[mol_idx];
    x_ref[0] = x[cg_idx][0];
    x_ref[1] = x[cg_idx][1];
    x_ref[2] = x[cg_idx][2];
    found_ref = 1;
  } else {
    // Use first atom in molecule as reference
    for (int i = 0; i < nlocal; i++) {
      if (molecule[i] == mol_id) {
        x_ref[0] = x[i][0];
        x_ref[1] = x[i][1];
        x_ref[2] = x[i][2];
        found_ref = 1;
        break;
      }
    }
  }
  
  // If no reference found, cannot calculate COM
  if (!found_ref) {
    x_com[mol_idx][0] = x_com[mol_idx][1] = x_com[mol_idx][2] = 0.0;
    v_com[mol_idx][0] = v_com[mol_idx][1] = v_com[mol_idx][2] = 0.0;
    mass_com[mol_idx] = 0.0;
    return;
  }
  
  // Calculate partial sums for all atoms in molecule
  for (int i = 0; i < nall; i++) {
    if (molecule[i] != mol_id) continue;
    
    // Get mass
    double massone;
    if (rmass) {
      massone = rmass[i];
    } else {
      massone = mass[type[i]];
    }
    
    // Unwrap coordinates relative to reference
    double dx[3], xunwrap[3];
    dx[0] = x[i][0] - x_ref[0];
    dx[1] = x[i][1] - x_ref[1];
    dx[2] = x[i][2] - x_ref[2];
    
    // Apply minimum image convention
    domain->minimum_image(FLERR, dx[0], dx[1], dx[2]);
    
    xunwrap[0] = x_ref[0] + dx[0];
    xunwrap[1] = x_ref[1] + dx[1];
    xunwrap[2] = x_ref[2] + dx[2];
    
    // Accumulate sums
    sum_mx[0] += xunwrap[0] * massone;
    sum_mx[1] += xunwrap[1] * massone;
    sum_mx[2] += xunwrap[2] * massone;
    
    sum_mv[0] += v[i][0] * massone;
    sum_mv[1] += v[i][1] * massone;
    sum_mv[2] += v[i][2] * massone;
    
    sum_m += massone;
  }
  
  // Aggregate across all processors using MPI
  double all_sum_mx[3], all_sum_mv[3], all_sum_m;
  MPI_Allreduce(sum_mx, all_sum_mx, 3, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(sum_mv, all_sum_mv, 3, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(&sum_m, &all_sum_m, 1, MPI_DOUBLE, MPI_SUM, world);
  
  // Calculate final COM
  if (all_sum_m > 0.0) {
    x_com[mol_idx][0] = all_sum_mx[0] / all_sum_m;
    x_com[mol_idx][1] = all_sum_mx[1] / all_sum_m;
    x_com[mol_idx][2] = all_sum_mx[2] / all_sum_m;
    
    v_com[mol_idx][0] = all_sum_mv[0] / all_sum_m;
    v_com[mol_idx][1] = all_sum_mv[1] / all_sum_m;
    v_com[mol_idx][2] = all_sum_mv[2] / all_sum_m;
    
    mass_com[mol_idx] = all_sum_m;
  } else {
    x_com[mol_idx][0] = x_com[mol_idx][1] = x_com[mol_idx][2] = 0.0;
    v_com[mol_idx][0] = v_com[mol_idx][1] = v_com[mol_idx][2] = 0.0;
    mass_com[mol_idx] = 0.0;
  }
}

/* ---------------------------------------------------------------------- */

void FixAdResSConstraint::calculate_com_from_cg(int mol_idx)
{
  // Get COM from CG particle position
  // Simple case: COM = CG particle position (CG particle represents the molecule)
  
  if (mol_idx < 0 || mol_idx >= nmolecules) return;
  
  int cg_idx = cg_atom_index[mol_idx];
  if (cg_idx < 0) {
    // No CG particle - cannot calculate COM from CG
    x_com[mol_idx][0] = x_com[mol_idx][1] = x_com[mol_idx][2] = 0.0;
    v_com[mol_idx][0] = v_com[mol_idx][1] = v_com[mol_idx][2] = 0.0;
    mass_com[mol_idx] = 0.0;
    return;
  }
  
  double **x = atom->x;
  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int nall = atom->nlocal + atom->nghost;
  
  if (cg_idx >= nall) {
    // CG particle index out of range
    x_com[mol_idx][0] = x_com[mol_idx][1] = x_com[mol_idx][2] = 0.0;
    v_com[mol_idx][0] = v_com[mol_idx][1] = v_com[mol_idx][2] = 0.0;
    mass_com[mol_idx] = 0.0;
    return;
  }
  
  // COM = CG particle position
  x_com[mol_idx][0] = x[cg_idx][0];
  x_com[mol_idx][1] = x[cg_idx][1];
  x_com[mol_idx][2] = x[cg_idx][2];
  
  // COM velocity = CG particle velocity
  v_com[mol_idx][0] = v[cg_idx][0];
  v_com[mol_idx][1] = v[cg_idx][1];
  v_com[mol_idx][2] = v[cg_idx][2];
  
  // Mass = CG particle mass
  if (rmass) {
    mass_com[mol_idx] = rmass[cg_idx];
  } else {
    mass_com[mol_idx] = mass[type[cg_idx]];
  }
}

/* ---------------------------------------------------------------------- */

void FixAdResSConstraint::constrain_atoms_to_com(int mol_idx, bool constrain_velocities)
{
  // Constrain atoms in molecule to follow COM position and velocity
  // Used for CG region where molecules move as rigid bodies
  
  if (mol_idx < 0 || mol_idx >= nmolecules) return;
  
  double **x = atom->x;
  double **v = atom->v;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  
  // Get COM position and velocity for this molecule
  double *x_com_mol = x_com[mol_idx];
  double *v_com_mol = v_com[mol_idx];
  
  // Iterate through all local atoms and constrain those in this molecule
  for (int i = 0; i < nlocal; i++) {
    // Check if atom belongs to this molecule
    if (molecule_map[i] != mol_idx) continue;
    
    // Skip CG particles themselves (they are not constrained)
    if (type[i] == cg_type) continue;
    
    // Constrain atom position: x_atom = x_com + displace
    x[i][0] = x_com_mol[0] + displace[i][0];
    x[i][1] = x_com_mol[1] + displace[i][1];
    x[i][2] = x_com_mol[2] + displace[i][2];
    
    // Constrain atom velocity: v_atom = v_com
    v[i][0] = v_com_mol[0];
    v[i][1] = v_com_mol[1];
    v[i][2] = v_com_mol[2];
  }
}

/* ---------------------------------------------------------------------- */

bool FixAdResSConstraint::is_molecule_in_cg_region(int mol_idx)
{
  // Check if molecule is in CG region by checking lambda of representative atom
  // Returns true if lambda < lambda_cg_threshold
  
  if (mol_idx < 0 || mol_idx >= nmolecules) return false;
  if (!fix_region) return false;
  
  // Try to use CG particle as representative if available
  int cg_idx = cg_atom_index[mol_idx];
  if (cg_idx >= 0 && cg_idx < atom->nlocal + atom->nghost) {
    double lambda = fix_region->get_lambda(cg_idx);
    return (lambda < lambda_cg_threshold);
  }
  
  // Otherwise, find first atom in molecule and check its lambda
  tagint mol_id = mol_id_list[mol_idx];
  if (mol_id == 0) return false;
  
  tagint *molecule = atom->molecule;
  int nlocal = atom->nlocal;
  
  for (int i = 0; i < nlocal; i++) {
    if (molecule[i] == mol_id) {
      double lambda = fix_region->get_lambda(i);
      return (lambda < lambda_cg_threshold);
    }
  }
  
  return false;
}

/* ---------------------------------------------------------------------- */

bool FixAdResSConstraint::is_molecule_in_atomistic_region(int mol_idx)
{
  // Check if molecule is in atomistic region by checking lambda of representative atom
  // Returns true if lambda > lambda_at_threshold
  
  if (mol_idx < 0 || mol_idx >= nmolecules) return false;
  if (!fix_region) return false;
  
  // Try to use CG particle as representative if available
  int cg_idx = cg_atom_index[mol_idx];
  if (cg_idx >= 0 && cg_idx < atom->nlocal + atom->nghost) {
    double lambda = fix_region->get_lambda(cg_idx);
    return (lambda > lambda_at_threshold);
  }
  
  // Otherwise, find first atom in molecule and check its lambda
  tagint mol_id = mol_id_list[mol_idx];
  if (mol_id == 0) return false;
  
  tagint *molecule = atom->molecule;
  int nlocal = atom->nlocal;
  
  for (int i = 0; i < nlocal; i++) {
    if (molecule[i] == mol_id) {
      double lambda = fix_region->get_lambda(i);
      return (lambda > lambda_at_threshold);
    }
  }
  
  return false;
}

/* ---------------------------------------------------------------------- */

void FixAdResSConstraint::update_cg_from_com(int mol_idx)
{
  // Update CG particle position and velocity from calculated COM
  // Used for atomistic region where CG particle tracks COM of atomistic atoms
  
  if (mol_idx < 0 || mol_idx >= nmolecules) return;
  
  int cg_idx = cg_atom_index[mol_idx];
  if (cg_idx < 0) {
    // No CG particle - cannot update
    if (comm->me == 0) {
      error->warning(FLERR, "fix adress/constraint: Cannot update CG particle for molecule {}: no CG particle found", mol_idx);
    }
    return;
  }
  
  int nall = atom->nlocal + atom->nghost;
  if (cg_idx >= nall) {
    // CG particle index out of range
    return;
  }
  
  double **x = atom->x;
  double **v = atom->v;
  imageint *image = atom->image;
  double *x_com_mol = x_com[mol_idx];
  double *v_com_mol = v_com[mol_idx];
  
  // Update CG particle position: x_cg = x_com
  x[cg_idx][0] = x_com_mol[0];
  x[cg_idx][1] = x_com_mol[1];
  x[cg_idx][2] = x_com_mol[2];
  
  // Update CG particle velocity: v_cg = v_com
  v[cg_idx][0] = v_com_mol[0];
  v[cg_idx][1] = v_com_mol[1];
  v[cg_idx][2] = v_com_mol[2];
  
  // Handle image flags for periodic boundaries
  domain->remap(x[cg_idx], image[cg_idx]);
}

/* ---------------------------------------------------------------------- */

void FixAdResSConstraint::grow_arrays(int nmax)
{
  memory->grow(molecule_map, nmax, "adress/constraint:molecule_map");
  memory->grow(displace, nmax, 3, "adress/constraint:displace");
  maxatom = nmax;
}

/* ---------------------------------------------------------------------- */

void FixAdResSConstraint::copy_arrays(int i, int j, int /*delflag*/)
{
  molecule_map[j] = molecule_map[i];
  displace[j][0] = displace[i][0];
  displace[j][1] = displace[i][1];
  displace[j][2] = displace[i][2];
}

/* ---------------------------------------------------------------------- */

int FixAdResSConstraint::pack_exchange(int i, double *buf)
{
  int n = 0;
  buf[n++] = molecule_map[i];
  buf[n++] = displace[i][0];
  buf[n++] = displace[i][1];
  buf[n++] = displace[i][2];
  return n;
}

/* ---------------------------------------------------------------------- */

int FixAdResSConstraint::unpack_exchange(int nlocal, double *buf)
{
  int n = 0;
  molecule_map[nlocal] = (int) buf[n++];
  displace[nlocal][0] = buf[n++];
  displace[nlocal][1] = buf[n++];
  displace[nlocal][2] = buf[n++];
  return n;
}

/* ---------------------------------------------------------------------- */

int FixAdResSConstraint::pack_restart(int i, double *buf)
{
  int n = 0;
  buf[n++] = molecule_map[i];
  buf[n++] = displace[i][0];
  buf[n++] = displace[i][1];
  buf[n++] = displace[i][2];
  return n;
}

/* ---------------------------------------------------------------------- */

void FixAdResSConstraint::unpack_restart(int nlocal, int nth)
{
  double **extra = atom->extra;
  int n = (nth - 1) * 4;
  molecule_map[nlocal] = (int) extra[nlocal][n++];
  displace[nlocal][0] = extra[nlocal][n++];
  displace[nlocal][1] = extra[nlocal][n++];
  displace[nlocal][2] = extra[nlocal][n++];
}

/* ---------------------------------------------------------------------- */

int FixAdResSConstraint::size_restart(int /*nlocal*/)
{
  return 4;
}

/* ---------------------------------------------------------------------- */

int FixAdResSConstraint::maxsize_restart()
{
  return 4;
}

/* ---------------------------------------------------------------------- */

double FixAdResSConstraint::memory_usage()
{
  double bytes = 0.0;
  bytes += (double) maxatom * sizeof(int);           // molecule_map
  bytes += (double) maxatom * 3 * sizeof(double);    // displace
  bytes += (double) maxmolecule * sizeof(int);        // cg_atom_index
  bytes += (double) maxmolecule * sizeof(tagint);     // mol_id_list
  bytes += (double) maxmolecule * 3 * sizeof(double); // x_com, v_com, f_com, x_cg_stored, v_cg_stored
  bytes += (double) maxmolecule * sizeof(double);     // mass_com
  return bytes;
}

/* ---------------------------------------------------------------------- */

void *FixAdResSConstraint::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str, "cg_type") == 0) {
    dim = 0;
    return &cg_type;
  }
  return nullptr;
}

