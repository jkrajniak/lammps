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
    f_atomistic(nullptr), f_cg(nullptr), nmax_force(0)
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

  // get fix adres/region if available
  if (id_fix_region && !fix_region) {
    fix_region = dynamic_cast<FixAdResSRegion *>(modify->get_fix_by_id(id_fix_region));
    if (!fix_region) error->all(FLERR, "Fix {} for pair adress does not exist", id_fix_region);
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

  // Save current forces
  double **f = atom->f;
  double **f_saved;
  memory->create(f_saved, nall, 3, "pair_adress:f_saved");
  for (int i = 0; i < nall; i++) {
    f_saved[i][0] = f[i][0];
    f_saved[i][1] = f[i][1];
    f_saved[i][2] = f[i][2];
  }

  // Zero forces for atomistic computation
  for (int i = 0; i < nall; i++) {
    f[i][0] = 0.0;
    f[i][1] = 0.0;
    f[i][2] = 0.0;
  }

  // Compute atomistic forces
  pair_atomistic->compute(eflag, vflag);
  
  // Save atomistic forces
  for (int i = 0; i < nall; i++) {
    f_atomistic[i][0] = f[i][0];
    f_atomistic[i][1] = f[i][1];
    f_atomistic[i][2] = f[i][2];
  }

  // Zero forces for CG computation
  for (int i = 0; i < nall; i++) {
    f[i][0] = 0.0;
    f[i][1] = 0.0;
    f[i][2] = 0.0;
  }

  // Compute CG forces
  pair_cg->compute(eflag, vflag);

  // Save CG forces
  for (int i = 0; i < nall; i++) {
    f_cg[i][0] = f[i][0];
    f_cg[i][1] = f[i][1];
    f_cg[i][2] = f[i][2];
  }

  // Interpolate forces based on lambda
  interpolate_forces();

  // Restore saved forces and add interpolated forces
  for (int i = 0; i < nall; i++) {
    f[i][0] = f_saved[i][0] + f[i][0];
    f[i][1] = f_saved[i][1] + f[i][1];
    f[i][2] = f_saved[i][2] + f[i][2];
  }

  memory->destroy(f_saved);

  // Interpolate energies and virials
        if (eflag) {
    eng_vdwl = 0.0;
    eng_coul = 0.0;
    // Note: Energy interpolation would need per-atom lambda, simplified here
    // For now, use average of atomistic and CG energies
    if (fix_region) {
      // Simplified: use average lambda for energy
      double lambda_avg = 0.5;  // Could be improved with per-atom lambda
      eng_vdwl = lambda_avg * pair_atomistic->eng_vdwl + (1.0 - lambda_avg) * pair_cg->eng_vdwl;
      eng_coul = lambda_avg * pair_atomistic->eng_coul + (1.0 - lambda_avg) * pair_cg->eng_coul;
    } else {
      eng_vdwl = 0.5 * (pair_atomistic->eng_vdwl + pair_cg->eng_vdwl);
      eng_coul = 0.5 * (pair_atomistic->eng_coul + pair_cg->eng_coul);
    }
  }

  if (vflag) {
    for (int n = 0; n < 6; n++) {
      if (fix_region) {
        double lambda_avg = 0.5;
        virial[n] = lambda_avg * pair_atomistic->virial[n] + (1.0 - lambda_avg) * pair_cg->virial[n];
      } else {
        virial[n] = 0.5 * (pair_atomistic->virial[n] + pair_cg->virial[n]);
      }
    }
  }

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

