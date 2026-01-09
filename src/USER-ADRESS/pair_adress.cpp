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

#include "pair_adress.h"

#include <cmath>

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "fix_adress_region.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "update.h"
#include "utils.h"

#include <cmath>
#include <cstring>
#include <mpi.h>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairAdResS::PairAdResS(LAMMPS *lmp) : Pair(lmp), cut_global(0.0), cut(nullptr),
    id_fix_region(nullptr), fix_region(nullptr), pair_atomistic(nullptr),
    pair_cg(nullptr), style_atomistic(nullptr), style_cg(nullptr),
    f_atomistic(nullptr), f_cg(nullptr), nmax_force(0), cg_type(0)
{
  restartinfo = 1;
  writedata = 1;
  allocated = 0;
}

/* ---------------------------------------------------------------------- */

PairAdResS::~PairAdResS()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
  }
  delete[] id_fix_region;
  delete[] style_atomistic;
  delete[] style_cg;
  if (pair_atomistic) delete pair_atomistic;
  if (pair_cg) delete pair_cg;
  memory->destroy(f_atomistic);
  memory->destroy(f_cg);
}

/* ---------------------------------------------------------------------- */

void PairAdResS::compute(int eflag, int vflag)
{
  // Check that sub-styles are initialized
  if (!pair_atomistic || !pair_cg)
    error->all(FLERR, "Pair adress sub-styles not initialized");

  // Get fix adress/region if available
  if (id_fix_region && !fix_region) {
    fix_region = dynamic_cast<FixAdResSRegion *>(modify->get_fix_by_id(id_fix_region));
    if (!fix_region) error->all(FLERR, "Fix {} for pair adress does not exist", id_fix_region);
  }

  if (!fix_region) {
    error->all(FLERR, "Pair adress requires fix adress/region for filter-during-compute approach");
  }

  // Allocate temporary force arrays if needed
  const int nall = atom->nlocal + atom->nghost;
  if (nall > nmax_force) {
    memory->destroy(f_atomistic);
    memory->destroy(f_cg);
    nmax_force = atom->nmax;
    memory->create(f_atomistic, nmax_force, 3, "pair_adress:f_atomistic");
    memory->create(f_cg, nmax_force, 3, "pair_adress:f_cg");
  }

  // Zero force arrays
  for (int i = 0; i < nall; i++) {
    f_atomistic[i][0] = 0.0;
    f_atomistic[i][1] = 0.0;
    f_atomistic[i][2] = 0.0;
    f_cg[i][0] = 0.0;
    f_cg[i][1] = 0.0;
    f_cg[i][2] = 0.0;
  }

  // Get CG type if not already set
  if (cg_type == 0) {
    cg_type = get_cg_type();
    if (cg_type == 0) {
      error->all(FLERR, "Pair adress requires CG type to be set via pair_coeff cg_type or fix adress/constraint");
    }
  }

  // Initialize energy and virial
  ev_init(eflag, vflag);
  
  // Zero energy and virial for sub-styles (we'll accumulate manually)
  pair_atomistic->eng_vdwl = 0.0;
  pair_atomistic->eng_coul = 0.0;
  pair_cg->eng_vdwl = 0.0;
  pair_cg->eng_coul = 0.0;
  for (int n = 0; n < 6; n++) {
    pair_atomistic->virial[n] = 0.0;
    pair_cg->virial[n] = 0.0;
  }

  // Get neighbor list (normal list, all atoms within cutoff)
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int inum = list->inum;
  
  double **x = atom->x;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  // Iterate over neighbor list and filter during computation
  for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii];
    
    int *jlist = firstneigh[i];
    int jnum = numneigh[i];
    
    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj] & NEIGHMASK;
      
      // Skip if j < i and newton_pair is off (half list)
      if (!newton_pair && j < i) continue;
      
      // Check distance (cutoff check)
      double delx = x[i][0] - x[j][0];
      double dely = x[i][1] - x[j][1];
      double delz = x[i][2] - x[j][2];
      double rsq = delx * delx + dely * dely + delz * delz;
      
      // Check if pair is within cutoff
      if (rsq > cutsq[type[i]][type[j]]) continue;
      
      // Filter: Check if pair should have AT forces
      if (should_compute_AT_force(i, j)) {
        compute_AT_force_pair(i, j, rsq, delx, dely, delz, eflag, vflag);
      }
      
      // Filter: Check if pair should have CG forces
      if (should_compute_CG_force(i, j)) {
        compute_CG_force_pair(i, j, rsq, delx, dely, delz, eflag, vflag);
      }
      
      // If neither, skip entirely (computational savings!)
    }
  }

  // Interpolate forces based on lambda
  interpolate_forces();

  // Interpolate energy and virial
  if (eflag) interpolate_energy();
  if (vflag) interpolate_virial();

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

void PairAdResS::settings(int narg, char **arg)
{
  // Syntax: pair_style adress cut atomistic_style atomistic_args ... cg_style cg_args ...
  // Minimum: cut + 2 style names
  if (narg < 3) error->all(FLERR, "Illegal pair_style adress command");

  cut_global = utils::numeric(FLERR, arg[0], false, lmp);
  if (cut_global <= 0.0) error->all(FLERR, "Illegal pair_style adress command");

  // Find where CG style starts (look for second pair style name)
  // Start from arg[2] to skip the atomistic style name at arg[1]
  int iarg = 2;
  while (iarg < narg && !force->pair_map->count(arg[iarg]) &&
         !lmp->match_style("pair", arg[iarg])) {
    iarg++;
  }

  if (iarg >= narg) error->all(FLERR, "Illegal pair_style adress command: missing CG pair style");

  // Create atomistic pair style
  style_atomistic = utils::strdup(arg[1]);
  int dummy = 0;
  pair_atomistic = force->new_pair(arg[1], 1, dummy);

  // Determine arguments for atomistic style (everything between style names)
  // arg[0] = cut_global
  // arg[1] = atomistic style name
  // arg[2] ... arg[iarg-1] = atomistic arguments
  // arg[iarg] = CG style name
  // So the number of atomistic arguments is (iarg - 2)
  int narg_at = iarg - 2;
  if (narg_at > 0) {
    // Pass narg_at arguments starting from arg[2]
    pair_atomistic->settings(narg_at, &arg[2]);
  }

  // Create CG pair style
  style_cg = utils::strdup(arg[iarg]);
  pair_cg = force->new_pair(arg[iarg], 1, dummy);

  // Determine arguments for CG style (everything after CG style name)
  if (iarg + 1 < narg) {
    pair_cg->settings(narg - iarg - 1, &arg[iarg + 1]);
  }
}

/* ---------------------------------------------------------------------- */

void PairAdResS::coeff(int narg, char **arg)
{
  // New syntax: pair_coeff i j atomistic <atomistic_args> cg <cg_args> [fix fix_id]
  // Old syntax: pair_coeff i j <args> [fix fix_id] (backward compatible)
  if (narg < 2) error->all(FLERR, "Incorrect args for pair coefficients");

  if (!allocated) allocate();

  // Find keywords: atomistic, cg, fix
  int atomistic_arg = -1;
  int cg_arg = -1;
  int fix_arg = -1;
  
  for (int i = 0; i < narg; i++) {
    if (strcmp(arg[i], "atomistic") == 0) {
      atomistic_arg = i;
    } else if (strcmp(arg[i], "cg") == 0) {
      cg_arg = i;
    } else if (strcmp(arg[i], "fix") == 0) {
      fix_arg = i;
    }
  }

  // Determine if using new syntax (both keywords present) or old syntax
  bool new_syntax = (atomistic_arg >= 0 && cg_arg >= 0);
  
  if (new_syntax) {
    // NEW SYNTAX: Separate coefficients for atomistic and CG
    
    // Validate keyword order
    if (atomistic_arg >= cg_arg) {
      error->all(FLERR, "pair_coeff: 'atomistic' keyword must come before 'cg' keyword");
    }
    if (fix_arg >= 0 && cg_arg >= fix_arg) {
      error->all(FLERR, "pair_coeff: 'cg' keyword must come before 'fix' keyword");
    }
    
    // Extract atomistic arguments (between "atomistic" and "cg")
    int narg_at = cg_arg - atomistic_arg - 1;
    if (narg_at < 0) {
      error->all(FLERR, "pair_coeff: No arguments provided for atomistic sub-style");
    }
    
    // Extract CG arguments (between "cg" and "fix" or end)
    int end_arg = (fix_arg >= 0) ? fix_arg : narg;
    int narg_cg = end_arg - cg_arg - 1;
    if (narg_cg < 0) {
      error->all(FLERR, "pair_coeff: No arguments provided for CG sub-style");
    }
    
    // Build argument arrays for each sub-style
    // Atomistic: [i, j, arg[atomistic_arg+1], ..., arg[cg_arg-1]]
    // CG: [i, j, arg[cg_arg+1], ..., arg[end_arg-1]]
    
    // For atomistic sub-style
    if (pair_atomistic && narg_at > 0) {
      // Create temporary array: [arg[0], arg[1], atomistic_args...]
      char **arg_at = new char*[narg_at + 2];
      arg_at[0] = arg[0];  // i
      arg_at[1] = arg[1];  // j
      for (int i = 0; i < narg_at; i++) {
        arg_at[i + 2] = arg[atomistic_arg + 1 + i];
      }
      pair_atomistic->coeff(narg_at + 2, arg_at);
      delete[] arg_at;
    }
    
    // For CG sub-style
    if (pair_cg && narg_cg > 0) {
      // Create temporary array: [arg[0], arg[1], cg_args...]
      char **arg_cg = new char*[narg_cg + 2];
      arg_cg[0] = arg[0];  // i
      arg_cg[1] = arg[1];  // j
      for (int i = 0; i < narg_cg; i++) {
        arg_cg[i + 2] = arg[cg_arg + 1 + i];
      }
      pair_cg->coeff(narg_cg + 2, arg_cg);
      delete[] arg_cg;
    }
    
  } else {
    // OLD SYNTAX: Same coefficients for both (backward compatibility)
    // Check if only one keyword is present (error case)
    if (atomistic_arg >= 0 || cg_arg >= 0) {
      error->all(FLERR, "pair_coeff: Both 'atomistic' and 'cg' keywords must be present, or neither");
    }
    
    // Original behavior: pass same coefficients to both sub-styles
    int narg_sub = fix_arg >= 0 ? fix_arg : narg;
    if (pair_atomistic && narg_sub > 0) {
      pair_atomistic->coeff(narg_sub, arg);
    }
    if (pair_cg && narg_sub > 0) {
      pair_cg->coeff(narg_sub, arg);
    }
  }

  // Extract fix ID if present
  if (fix_arg >= 0 && fix_arg + 1 < narg) {
    delete[] id_fix_region;
    id_fix_region = utils::strdup(arg[fix_arg + 1]);
  }

  // Look for CG type specification: "cg_type <type>"
  // This can appear anywhere in the arguments
  for (int i = 0; i < narg; i++) {
    if (strcmp(arg[i], "cg_type") == 0 && i + 1 < narg) {
      cg_type = utils::numeric(FLERR, arg[i + 1], false, lmp);
      if (cg_type < 1 || cg_type > atom->ntypes)
        error->all(FLERR, "Invalid CG atom type {} in pair_coeff", cg_type);
      break;
    }
  }

  // Set cutoffs and flags (use maximum of both styles)
  int ilo, ihi, jlo, jhi;
  if (strcmp(arg[0], "*") == 0 && strcmp(arg[1], "*") == 0) {
    ilo = 1; ihi = atom->ntypes;
    jlo = 1; jhi = atom->ntypes;
  } else {
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);
  }

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      // Use maximum cutoff from both styles
      // Base Pair class has cutforce (max cutoff) and cutsq[i][j] (cutoff squared)
      double cut_at = cut_global;
      double cut_cg = cut_global;
      if (pair_atomistic) {
        if (pair_atomistic->cutsq && pair_atomistic->cutsq[i][j] > 0.0)
          cut_at = sqrt(pair_atomistic->cutsq[i][j]);
        else
          cut_at = pair_atomistic->cutforce;
      }
      if (pair_cg) {
        if (pair_cg->cutsq && pair_cg->cutsq[i][j] > 0.0)
          cut_cg = sqrt(pair_cg->cutsq[i][j]);
        else
          cut_cg = pair_cg->cutforce;
      }
      cut[i][j] = MAX(cut_at, cut_cg);
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

/* ---------------------------------------------------------------------- */

void PairAdResS::init_style()
{
  // Initialize sub-styles
  if (pair_atomistic) pair_atomistic->init_style();
  if (pair_cg) pair_cg->init_style();

  // Validate that both sub-styles support pairwise potentials only
  // AdResS is limited to pairwise potentials (three-body, many-body not supported)
  if (pair_atomistic) {
    // Check for known many-body styles
    const char *many_body_styles[] = {
      "tersoff", "sw", "eam", "meam", "pace", "threebody", "fourbody",
      "airebo", "rebo", "airebo/morse", "reax", "reax/c"
    };
    for (int i = 0; i < sizeof(many_body_styles)/sizeof(many_body_styles[0]); i++) {
      if (strcmp(style_atomistic, many_body_styles[i]) == 0) {
        error->all(FLERR,
          "pair_style adress: Atomistic pair style '{}' is a many-body potential.\n"
          "AdResS only works with pairwise potentials that support the single() method.\n"
          "Supported styles include: lj/cut, lj/cut/coul/cut, table, coul/cut, etc.\n"
          "Three-body, four-body, or many-body potentials are not supported.",
          style_atomistic);
      }
    }
  }

  if (pair_cg) {
    // Check for known many-body styles
    const char *many_body_styles[] = {
      "tersoff", "sw", "eam", "meam", "pace", "threebody", "fourbody",
      "airebo", "rebo", "airebo/morse", "reax", "reax/c"
    };
    for (int i = 0; i < sizeof(many_body_styles)/sizeof(many_body_styles[0]); i++) {
      if (strcmp(style_cg, many_body_styles[i]) == 0) {
        error->all(FLERR,
          "pair_style adress: CG pair style '{}' is a many-body potential.\n"
          "AdResS only works with pairwise potentials that support the single() method.\n"
          "Supported styles include: lj/cut, lj/cut/coul/cut, table, coul/cut, etc.\n"
          "Three-body, four-body, or many-body potentials are not supported.",
          style_cg);
      }
    }
  }

  // request regular neighbor list
  neighbor->add_request(this);

  // look for fix adres/region
  if (id_fix_region) {
    fix_region = dynamic_cast<FixAdResSRegion *>(modify->get_fix_by_id(id_fix_region));
    if (!fix_region) error->all(FLERR, "Fix {} for pair adress does not exist", id_fix_region);
  }
}

/* ---------------------------------------------------------------------- */

double PairAdResS::init_one(int i, int j)
{
  // Initialize sub-styles first
  if (pair_atomistic) pair_atomistic->init_one(i, j);
  if (pair_cg) pair_cg->init_one(i, j);

  if (setflag[i][j] == 0) {
    // Use maximum cutoff from sub-styles
    // Base Pair class has cutforce (max cutoff) and cutsq[i][j] (cutoff squared)
    double cut_at = cut_global;
    double cut_cg = cut_global;
    if (pair_atomistic) {
      if (pair_atomistic->cutsq && pair_atomistic->cutsq[i][j] > 0.0)
        cut_at = sqrt(pair_atomistic->cutsq[i][j]);
      else
        cut_at = pair_atomistic->cutforce;
    }
    if (pair_cg) {
      if (pair_cg->cutsq && pair_cg->cutsq[i][j] > 0.0)
        cut_cg = sqrt(pair_cg->cutsq[i][j]);
      else
        cut_cg = pair_cg->cutforce;
    }
    cut[i][j] = MAX(cut_at, cut_cg);
  }

  if (cut[i][j] > 0.0)
    cutsq[i][j] = cut[i][j] * cut[i][j];
  else
    cutsq[i][j] = 0.0;

  return cut[i][j];
}

/* ---------------------------------------------------------------------- */

void PairAdResS::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");
  memory->create(cut, n + 1, n + 1, "pair:cut");

  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++) setflag[i][j] = 0;
}

/* ---------------------------------------------------------------------- */

double PairAdResS::switching_function(double r, double rcut, double lambda)
{
  // Linear switching function based on lambda
  // lambda = 1: full atomistic, lambda = 0: full CG
  if (lambda >= 1.0) return 1.0;
  if (lambda <= 0.0) return 0.0;
  return lambda;
}

/* ---------------------------------------------------------------------- */

void PairAdResS::interpolate_forces()
{
  // Interpolate forces based on per-atom lambda values
  // F = lambda * F_at + (1 - lambda) * F_cg
  
  double **f = atom->f;
  const int nlocal = atom->nlocal;
  const int nall = nlocal + atom->nghost;

  if (!fix_region) {
    // No lambda values - use equal weighting
    for (int i = 0; i < nall; i++) {
      f[i][0] = 0.5 * (f_atomistic[i][0] + f_cg[i][0]);
      f[i][1] = 0.5 * (f_atomistic[i][1] + f_cg[i][1]);
      f[i][2] = 0.5 * (f_atomistic[i][2] + f_cg[i][2]);
    }
    return;
  }

  // Interpolate based on per-atom lambda
  for (int i = 0; i < nall; i++) {
    double lambda_i = fix_region->get_lambda(i);
    
    f[i][0] = lambda_i * f_atomistic[i][0] + (1.0 - lambda_i) * f_cg[i][0];
    f[i][1] = lambda_i * f_atomistic[i][1] + (1.0 - lambda_i) * f_cg[i][1];
    f[i][2] = lambda_i * f_atomistic[i][2] + (1.0 - lambda_i) * f_cg[i][2];
  }
}

/* ---------------------------------------------------------------------- */

int PairAdResS::get_cg_type()
{
  // If already set via pair_coeff, return it
  if (cg_type > 0) return cg_type;

  // Try to get CG type from fix adress/constraint
  // Search through all fixes to find fix adress/constraint
  if (modify) {
    for (int i = 0; i < modify->nfix; i++) {
      Fix *fix = modify->fix[i];
      if (strcmp(fix->style, "adress/constraint") == 0) {
        // Try to extract cg_type using extract method
        int dim = 0;
        int *cg_type_ptr = (int *) fix->extract("cg_type", dim);
        if (cg_type_ptr) {
          cg_type = *cg_type_ptr;
          return cg_type;
        }
      }
    }
  }
  return 0;  // Not found, user must specify via pair_coeff cg_type
}

/* ---------------------------------------------------------------------- */

void PairAdResS::filter_forces_by_region()
{
  // Filter forces based on region to ensure correct physics
  // NOTE: This does NOT provide computational savings - forces are already computed.
  // It only zeros forces that shouldn't exist based on region.
  // 
  // For true computational efficiency, we need to skip pair calculations
  // entirely (see REGION_FILTERING_CORRECT_APPROACH.md)
  //
  // This filtering ensures:
  // - AT forces: Only for AT atoms in Transition + Atomistic regions
  // - CG forces: Only for CG particles in CG region

  if (!fix_region) {
    // No region information - can't filter
    // This means forces are computed everywhere (inefficient but works)
    return;
  }

  const int nall = atom->nlocal + atom->nghost;
  int *type = atom->type;
  const int RES_ATOMISTIC = 0;
  const int RES_TRANSITION = 1;
  const int RES_CG = 2;

  // Filter AT forces:
  // - Zero AT forces for atoms in CG region
  // - Zero AT forces for CG particles (they don't have AT interactions)
  for (int i = 0; i < nall; i++) {
    int res_i = fix_region->get_resolution(i);
    bool is_cg = (cg_type > 0 && type[i] == cg_type);

    // Zero AT forces if:
    // 1. Atom is in CG region, OR
    // 2. Atom is a CG particle
    if (res_i == RES_CG || is_cg) {
      f_atomistic[i][0] = 0.0;
      f_atomistic[i][1] = 0.0;
      f_atomistic[i][2] = 0.0;
    }
  }

  // Filter CG forces:
  // - Zero CG forces for AT atoms (CG forces only act on CG particles)
  // - Zero CG forces for CG particles not in CG region
  for (int i = 0; i < nall; i++) {
    int res_i = fix_region->get_resolution(i);
    bool is_cg = (cg_type > 0 && type[i] == cg_type);

    // Zero CG forces if:
    // 1. Atom is not a CG particle, OR
    // 2. CG particle is not in CG region
    if (!is_cg || res_i != RES_CG) {
      f_cg[i][0] = 0.0;
      f_cg[i][1] = 0.0;
      f_cg[i][2] = 0.0;
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairAdResS::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i, j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j], sizeof(int), 1, fp);
      if (setflag[i][j]) fwrite(&cut[i][j], sizeof(double), 1, fp);
    }
}

/* ---------------------------------------------------------------------- */

void PairAdResS::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i, j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR, &setflag[i][j], sizeof(int), 1, fp, nullptr, error);
      MPI_Bcast(&setflag[i][j], 1, MPI_INT, 0, world);
      if (setflag[i][j]) {
        if (me == 0) utils::sfread(FLERR, &cut[i][j], sizeof(double), 1, fp, nullptr, error);
        MPI_Bcast(&cut[i][j], 1, MPI_DOUBLE, 0, world);
      }
    }
}

/* ---------------------------------------------------------------------- */

void PairAdResS::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global, sizeof(double), 1, fp);
  fwrite(&mix_flag, sizeof(int), 1, fp);
}

/* ---------------------------------------------------------------------- */

void PairAdResS::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR, &cut_global, sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &mix_flag, sizeof(int), 1, fp, nullptr, error);
  }
  MPI_Bcast(&cut_global, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&mix_flag, 1, MPI_INT, 0, world);
}

/* ---------------------------------------------------------------------- */

double PairAdResS::single(int /*i*/, int /*j*/, int itype, int jtype, double rsq,
                         double /*factor_coul*/, double factor_lj, double &fforce)
{
  double r = sqrt(rsq);
  double lambda_ij = 0.5;    // default
  double sw = switching_function(r, cut[itype][jtype], lambda_ij);
  fforce = 0.0;
  return sw * factor_lj;
}

/* ---------------------------------------------------------------------- */

void *PairAdResS::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str, "cut") == 0) {
    dim = 2;
    return (void *) cut;
  }
  return nullptr;
}

/* ----------------------------------------------------------------------
   Check if pair (i, j) should have AT forces
   Returns true if:
   - Both atoms are NOT CG particles
   - Both atoms are in Transition OR Atomistic region (not CG region)
------------------------------------------------------------------------- */

bool PairAdResS::should_compute_AT_force(int i, int j)
{
  if (!fix_region) return false;
  
  int *type = atom->type;
  int res_i = fix_region->get_resolution(i);
  int res_j = fix_region->get_resolution(j);
  
  // Both atoms must NOT be CG particles
  bool is_cg_i = (cg_type > 0 && type[i] == cg_type);
  bool is_cg_j = (cg_type > 0 && type[j] == cg_type);
  if (is_cg_i || is_cg_j) return false;
  
  // Both atoms must be in Transition OR Atomistic region (not CG region)
  const int RES_ATOMISTIC = 0;
  const int RES_TRANSITION = 1;
  const int RES_CG = 2;
  
  if (res_i == RES_CG || res_j == RES_CG) return false;
  
  return true;
}

/* ----------------------------------------------------------------------
   Check if pair (i, j) should have CG forces
   Returns true if:
   - Both atoms ARE CG particles
   - Both atoms are in CG region
------------------------------------------------------------------------- */

bool PairAdResS::should_compute_CG_force(int i, int j)
{
  if (!fix_region) return false;
  
  int *type = atom->type;
  int res_i = fix_region->get_resolution(i);
  int res_j = fix_region->get_resolution(j);
  
  // Both atoms must BE CG particles
  bool is_cg_i = (cg_type > 0 && type[i] == cg_type);
  bool is_cg_j = (cg_type > 0 && type[j] == cg_type);
  if (!is_cg_i || !is_cg_j) return false;
  
  // Both atoms must be in CG region
  const int RES_CG = 2;
  if (res_i != RES_CG || res_j != RES_CG) return false;
  
  return true;
}

/* ----------------------------------------------------------------------
   Compute AT force for a single pair (i, j)
   Uses single() method if available, otherwise falls back to full computation
------------------------------------------------------------------------- */

void PairAdResS::compute_AT_force_pair(int i, int j, double rsq, 
                                       double delx, double dely, double delz,
                                       int eflag, int vflag)
{
  int *type = atom->type;
  int itype = type[i];
  int jtype = type[j];
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  
  // Check if single() is available
  if (pair_atomistic->single_enable == 0) {
    // Fall back: This shouldn't happen often, but if it does, we need full computation
    // For now, skip this pair (could implement fallback later)
    static int warned = 0;
    if (!warned && comm->me == 0) {
      error->warning(FLERR, 
        "Pair style {} does not support single() method. "
        "Some AT forces may not be computed correctly.",
        pair_atomistic->style);
      warned = 1;
    }
    return;
  }
  
  // Get special bond factors
  double *special_lj = force->special_lj;
  double *special_coul = force->special_coul;
  
  // Check for special bonds (1-2, 1-3, 1-4 neighbors)
  tagint *tag = atom->tag;
  tagint **special = atom->special;
  int **nspecial = atom->nspecial;
  
  int which = 0;
  if (atom->molecular != Atom::ATOMIC) {
    // Find if j is in special list of i
    for (int k = 0; k < nspecial[i][0]; k++) {
      if (special[i][k] == tag[j]) {
        if (k < nspecial[i][1]) which = 1;  // 1-2 neighbor
        else if (k < nspecial[i][2]) which = 2;  // 1-3 neighbor
        else which = 3;  // 1-4 neighbor
        break;
      }
    }
  }
  
  double factor_lj = special_lj[which];
  double factor_coul = special_coul[which];
  
  // Compute force using single() method
  double fpair;
  double energy = pair_atomistic->single(i, j, itype, jtype, rsq, 
                                          factor_coul, factor_lj, fpair);
  
  // Compute force vector
  double r = sqrt(rsq);
  if (r > 0.0) {
    double fforce = fpair / r;
    
    f_atomistic[i][0] += delx * fforce;
    f_atomistic[i][1] += dely * fforce;
    f_atomistic[i][2] += delz * fforce;
    
    if (newton_pair || j < nlocal) {
      f_atomistic[j][0] -= delx * fforce;
      f_atomistic[j][1] -= dely * fforce;
      f_atomistic[j][2] -= delz * fforce;
    }
    
    // Tally energy and virial
    if (eflag) {
      pair_atomistic->eng_vdwl += energy;
      // Note: eng_coul would be added here if coulombic interactions exist
    }
    
    if (vflag) {
      // Tally virial manually (ev_tally would tally to this pair's virial, not sub-style's)
      // Virial = r_ij * F_ij (force times distance)
      if (vflag_global) {
        pair_atomistic->virial[0] += delx * delx * fpair;
        pair_atomistic->virial[1] += dely * dely * fpair;
        pair_atomistic->virial[2] += delz * delz * fpair;
        pair_atomistic->virial[3] += delx * dely * fpair;
        pair_atomistic->virial[4] += delx * delz * fpair;
        pair_atomistic->virial[5] += dely * delz * fpair;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   Compute CG force for a single pair (i, j)
   Uses single() method if available, otherwise falls back to full computation
------------------------------------------------------------------------- */

void PairAdResS::compute_CG_force_pair(int i, int j, double rsq,
                                        double delx, double dely, double delz,
                                        int eflag, int vflag)
{
  int *type = atom->type;
  int itype = type[i];
  int jtype = type[j];
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  
  // Check if single() is available
  if (pair_cg->single_enable == 0) {
    // Fall back: This shouldn't happen often, but if it does, we need full computation
    static int warned = 0;
    if (!warned && comm->me == 0) {
      error->warning(FLERR,
        "Pair style {} does not support single() method. "
        "Some CG forces may not be computed correctly.",
        pair_cg->style);
      warned = 1;
    }
    return;
  }
  
  // Get special bond factors
  double *special_lj = force->special_lj;
  double *special_coul = force->special_coul;
  
  // Check for special bonds (1-2, 1-3, 1-4 neighbors)
  tagint *tag = atom->tag;
  tagint **special = atom->special;
  int **nspecial = atom->nspecial;
  
  int which = 0;
  if (atom->molecular != Atom::ATOMIC) {
    // Find if j is in special list of i
    for (int k = 0; k < nspecial[i][0]; k++) {
      if (special[i][k] == tag[j]) {
        if (k < nspecial[i][1]) which = 1;  // 1-2 neighbor
        else if (k < nspecial[i][2]) which = 2;  // 1-3 neighbor
        else which = 3;  // 1-4 neighbor
        break;
      }
    }
  }
  
  double factor_lj = special_lj[which];
  double factor_coul = special_coul[which];
  
  // Compute force using single() method
  double fpair;
  double energy = pair_cg->single(i, j, itype, jtype, rsq,
                                   factor_coul, factor_lj, fpair);
  
  // Compute force vector
  double r = sqrt(rsq);
  if (r > 0.0) {
    double fforce = fpair / r;
    
    f_cg[i][0] += delx * fforce;
    f_cg[i][1] += dely * fforce;
    f_cg[i][2] += delz * fforce;
    
    if (newton_pair || j < nlocal) {
      f_cg[j][0] -= delx * fforce;
      f_cg[j][1] -= dely * fforce;
      f_cg[j][2] -= delz * fforce;
    }
    
    // Tally energy and virial
    if (eflag) {
      pair_cg->eng_vdwl += energy;
      // Note: eng_coul would be added here if coulombic interactions exist
    }
    
    if (vflag) {
      // Tally virial manually (ev_tally would tally to this pair's virial, not sub-style's)
      // Virial = r_ij * F_ij (force times distance)
      if (vflag_global) {
        pair_cg->virial[0] += delx * delx * fpair;
        pair_cg->virial[1] += dely * dely * fpair;
        pair_cg->virial[2] += delz * delz * fpair;
        pair_cg->virial[3] += delx * dely * fpair;
        pair_cg->virial[4] += delx * delz * fpair;
        pair_cg->virial[5] += dely * delz * fpair;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   Interpolate energy based on per-atom lambda values
------------------------------------------------------------------------- */

void PairAdResS::interpolate_energy()
{
  if (!fix_region) {
    // No lambda values - use equal weighting
    eng_vdwl = 0.5 * (pair_atomistic->eng_vdwl + pair_cg->eng_vdwl);
    eng_coul = 0.5 * (pair_atomistic->eng_coul + pair_cg->eng_coul);
    return;
  }
  
  // For now, use average lambda for energy interpolation
  // Could be improved with per-atom lambda weighting
  double lambda_avg = 0.5;
  eng_vdwl = lambda_avg * pair_atomistic->eng_vdwl + (1.0 - lambda_avg) * pair_cg->eng_vdwl;
  eng_coul = lambda_avg * pair_atomistic->eng_coul + (1.0 - lambda_avg) * pair_cg->eng_coul;
}

/* ----------------------------------------------------------------------
   Interpolate virial based on per-atom lambda values
------------------------------------------------------------------------- */

void PairAdResS::interpolate_virial()
{
  if (!fix_region) {
    // No lambda values - use equal weighting
    for (int n = 0; n < 6; n++) {
      virial[n] = 0.5 * (pair_atomistic->virial[n] + pair_cg->virial[n]);
    }
    return;
  }
  
  // For now, use average lambda for virial interpolation
  // Could be improved with per-atom lambda weighting
  double lambda_avg = 0.5;
  for (int n = 0; n < 6; n++) {
    virial[n] = lambda_avg * pair_atomistic->virial[n] + (1.0 - lambda_avg) * pair_cg->virial[n];
  }
}
