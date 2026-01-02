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

#include "pair_adres.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "fix_adress_region.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "utils.h"

#include <cmath>
#include <cstring>
#include <mpi.h>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairAdResS::PairAdRes(LAMMPS *lmp) : Pair(lmp), cut_global(0.0), cut(nullptr),
    id_fix_region(nullptr), fix_region(nullptr)
{
  restartinfo = 1;
  writedata = 1;
  allocated = 0;
}

/* ---------------------------------------------------------------------- */

PairAdResS::~PairAdRes()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
  }
  delete[] id_fix_region;
}

/* ---------------------------------------------------------------------- */

void PairAdResS::compute(int eflag, int vflag)
{
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz, evdwl, fpair;
  double rsq, r, factor_lj;
  int *ilist, *jlist, *numneigh, **firstneigh;

  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // get fix adres/region if available
  if (id_fix_region && !fix_region) {
    fix_region = dynamic_cast<FixAdResSRegion *>(modify->get_fix_by_id(id_fix_region));
    if (!fix_region) error->all(FLERR, "Fix {} for pair adress does not exist", id_fix_region);
  }

  // loop over neighbors of my atoms
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r = sqrt(rsq);

        // Get lambda values for switching
        double lambda_i = 1.0;
        double lambda_j = 1.0;
        if (fix_region) {
          lambda_i = fix_region->get_lambda(i);
          lambda_j = fix_region->get_lambda(j);
        }

        // Calculate effective lambda for pair interaction
        double lambda_ij = 0.5 * (lambda_i + lambda_j);

        // Apply switching function
        double sw = switching_function(r, cut[itype][jtype], lambda_ij);

        // Calculate force (simplified - would use actual pair potential)
        // This is a placeholder that applies switching
        fpair = 0.0;
        if (r > 0.0) {
          // Simplified force calculation
          double invr = 1.0 / r;
          fpair = sw * invr * invr;    // placeholder force
        }

        f[i][0] += delx * fpair;
        f[i][1] += dely * fpair;
        f[i][2] += delz * fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx * fpair;
          f[j][1] -= dely * fpair;
          f[j][2] -= delz * fpair;
        }

        if (eflag) {
          evdwl = sw * 1.0;    // placeholder energy
          evdwl *= factor_lj;
        }

        if (evflag) ev_tally(i, j, nlocal, newton_pair, evdwl, 0.0, fpair, delx, dely, delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

void PairAdResS::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR, "Illegal pair_style command");

  cut_global = utils::numeric(FLERR, arg[0], false, lmp);
  if (cut_global <= 0.0) error->all(FLERR, "Illegal pair_style command");
}

/* ---------------------------------------------------------------------- */

void PairAdResS::coeff(int narg, char **arg)
{
  if (narg < 2 || narg > 4) error->all(FLERR, "Incorrect args for pair coefficients");

  if (!allocated) allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  double cut_one = cut_global;
  if (narg >= 3) cut_one = utils::numeric(FLERR, arg[2], false, lmp);

  // optional fix adres/region ID
  if (narg == 4) {
    delete[] id_fix_region;
    id_fix_region = utils::strdup(arg[3]);
  }

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

/* ---------------------------------------------------------------------- */

void PairAdResS::init_style()
{
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
  if (setflag[i][j] == 0) {
    if (cut[i][i] > 0.0 && cut[j][j] > 0.0)
      cut[i][j] = 0.5 * (cut[i][i] + cut[j][j]);
    else
      cut[i][j] = cut_global;
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

