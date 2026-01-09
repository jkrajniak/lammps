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
#include "error.h"
#include "fix_adress_region.h"
#include "memory.h"
#include "modify.h"
#include "update.h"
#include "utils.h"

#include <cmath>
#include <cstring>
#include <map>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixAdResSConstraint::FixAdResSConstraint(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), id_fix_region(nullptr), fix_region(nullptr), cg_type(0),
    molecule_map(nullptr), cg_atom_index(nullptr), nmolecules(0), maxmolecule(0),
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
  // For now, assume each processor has complete molecule information
  // TODO: Handle distributed molecules properly

  // Allocate molecule arrays if needed
  if (nmolecules > maxmolecule) {
    maxmolecule = nmolecules + 100;  // Add some padding
    memory->destroy(cg_atom_index);
    memory->destroy(x_com);
    memory->destroy(v_com);
    memory->destroy(f_com);
    memory->destroy(x_cg_stored);
    memory->destroy(v_cg_stored);
    memory->destroy(mass_com);

    memory->create(cg_atom_index, maxmolecule, "adress/constraint:cg_atom_index");
    memory->create(x_com, maxmolecule, 3, "adress/constraint:x_com");
    memory->create(v_com, maxmolecule, 3, "adress/constraint:v_com");
    memory->create(f_com, maxmolecule, 3, "adress/constraint:f_com");
    memory->create(x_cg_stored, maxmolecule, 3, "adress/constraint:x_cg_stored");
    memory->create(v_cg_stored, maxmolecule, 3, "adress/constraint:v_cg_stored");
    memory->create(mass_com, maxmolecule, "adress/constraint:mass_com");

    // Initialize
    for (int m = 0; m < maxmolecule; m++) {
      cg_atom_index[m] = -1;
      x_com[m][0] = x_com[m][1] = x_com[m][2] = 0.0;
      v_com[m][0] = v_com[m][1] = v_com[m][2] = 0.0;
      f_com[m][0] = f_com[m][1] = f_com[m][2] = 0.0;
      x_cg_stored[m][0] = x_cg_stored[m][1] = x_cg_stored[m][2] = 0.0;
      v_cg_stored[m][0] = v_cg_stored[m][1] = v_cg_stored[m][2] = 0.0;
      mass_com[m] = 0.0;
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
  for (auto &pair : mol_id_to_cg) {
    tagint mol_id = pair.first;
    int cg_local = pair.second;
    auto it = mol_id_to_index.find(mol_id);
    if (it != mol_id_to_index.end()) {
      cg_atom_index[it->second] = cg_local;
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
  // Placeholder - will be implemented in Phase 2
  // For now, just verify fix_region is available
  if (!fix_region) return;
}

/* ---------------------------------------------------------------------- */

void FixAdResSConstraint::post_force(int /*vflag*/)
{
  // Placeholder - will be implemented in Phase 2-4
  // For now, just verify fix_region is available
  if (!fix_region) return;
}

/* ---------------------------------------------------------------------- */

void FixAdResSConstraint::final_integrate()
{
  // Placeholder - will be implemented in Phase 2
  // For now, just verify fix_region is available
  if (!fix_region) return;
}

/* ---------------------------------------------------------------------- */

void FixAdResSConstraint::calculate_com_from_atoms(int mol_idx)
{
  // Placeholder - will be implemented in Phase 3
  // Calculate COM from atom positions for given molecule
}

/* ---------------------------------------------------------------------- */

void FixAdResSConstraint::calculate_com_from_cg(int mol_idx)
{
  // Placeholder - will be implemented in Phase 2
  // Get COM from CG particle position
}

/* ---------------------------------------------------------------------- */

void FixAdResSConstraint::constrain_atoms_to_com(int mol_idx)
{
  // Placeholder - will be implemented in Phase 2
  // Constrain atoms to follow COM
}

/* ---------------------------------------------------------------------- */

void FixAdResSConstraint::update_cg_from_com(int mol_idx)
{
  // Placeholder - will be implemented in Phase 3
  // Update CG particle position from COM
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

