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

#include "fix_adres_region.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "memory.h"
#include "region.h"
#include "update.h"
#include "utils.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixAdResRegion::FixAdResRegion(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), id_region_atomistic(nullptr), id_region_transition(nullptr),
    id_region_cg(nullptr), region_atomistic(nullptr), region_transition(nullptr),
    region_cg(nullptr), resolution(nullptr), lambda(nullptr), maxatom(0)
{
  if (narg < 6) utils::missing_cmd_args(FLERR, "fix adres/region", error);

  restart_peratom = 1;
  restart_global = 0;

  id_region_atomistic = utils::strdup(arg[3]);
  id_region_transition = utils::strdup(arg[4]);
  id_region_cg = utils::strdup(arg[5]);

  // allocate per-atom arrays
  maxatom = atom->nmax;
  memory->create(resolution, maxatom, "adres/region:resolution");
  memory->create(lambda, maxatom, "adres/region:lambda");

  // initialize arrays
  for (int i = 0; i < maxatom; i++) {
    resolution[i] = RES_CG;    // default to CG
    lambda[i] = 0.0;
  }
}

/* ---------------------------------------------------------------------- */

FixAdResRegion::~FixAdResRegion()
{
  delete[] id_region_atomistic;
  delete[] id_region_transition;
  delete[] id_region_cg;
  memory->destroy(resolution);
  memory->destroy(lambda);
}

/* ---------------------------------------------------------------------- */

int FixAdResRegion::setmask()
{
  int mask = 0;
  mask |= POST_NEIGHBOR;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAdResRegion::init()
{
  // get region pointers
  region_atomistic = domain->get_region_by_id(id_region_atomistic);
  if (!region_atomistic)
    error->all(FLERR, "Region {} for fix adres/region does not exist", id_region_atomistic);

  region_transition = domain->get_region_by_id(id_region_transition);
  if (!region_transition)
    error->all(FLERR, "Region {} for fix adres/region does not exist", id_region_transition);

  region_cg = domain->get_region_by_id(id_region_cg);
  if (!region_cg)
    error->all(FLERR, "Region {} for fix adres/region does not exist", id_region_cg);
}

/* ---------------------------------------------------------------------- */

void FixAdResRegion::setup(int /*vflag*/)
{
  update_resolution();
}

/* ---------------------------------------------------------------------- */

void FixAdResRegion::post_neighbor()
{
  update_resolution();
}

/* ---------------------------------------------------------------------- */

void FixAdResRegion::update_resolution()
{
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // update regions before matching
  if (region_atomistic) region_atomistic->prematch();
  if (region_transition) region_transition->prematch();
  if (region_cg) region_cg->prematch();

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      int reg = determine_region(x[i][0], x[i][1], x[i][2]);
      resolution[i] = reg;

      if (reg == RES_TRANSITION) {
        lambda[i] = calculate_lambda(x[i][0], x[i][1], x[i][2]);
      } else if (reg == RES_ATOMISTIC) {
        lambda[i] = 1.0;
      } else {
        lambda[i] = 0.0;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

int FixAdResRegion::determine_region(double x, double y, double z)
{
  // Check regions in order: atomistic, transition, CG
  // Transition region takes precedence if atom is in both transition and atomistic/CG
  if (region_transition && region_transition->match(x, y, z)) {
    return RES_TRANSITION;
  }
  if (region_atomistic && region_atomistic->match(x, y, z)) {
    return RES_ATOMISTIC;
  }
  if (region_cg && region_cg->match(x, y, z)) {
    return RES_CG;
  }
  // Default to CG if not in any region
  return RES_CG;
}

/* ---------------------------------------------------------------------- */

double FixAdResRegion::calculate_lambda(double x, double y, double z)
{
  // Simple linear interpolation based on distance from atomistic region
  // This is a placeholder - actual implementation would depend on region geometry
  // For now, return 0.5 as a default transition value
  if (region_atomistic && region_atomistic->match(x, y, z)) {
    return 1.0;
  }
  if (region_cg && region_cg->match(x, y, z)) {
    return 0.0;
  }
  // In transition region, interpolate
  return 0.5;
}

/* ---------------------------------------------------------------------- */

int FixAdResRegion::get_resolution(int i) const
{
  if (i < 0 || i >= atom->nmax) return RES_CG;
  return resolution[i];
}

/* ---------------------------------------------------------------------- */

double FixAdResRegion::get_lambda(int i) const
{
  if (i < 0 || i >= atom->nmax) return 0.0;
  return lambda[i];
}

/* ---------------------------------------------------------------------- */

void FixAdResRegion::grow_arrays(int nmax)
{
  memory->grow(resolution, nmax, "adres/region:resolution");
  memory->grow(lambda, nmax, "adres/region:lambda");
  maxatom = nmax;
}

/* ---------------------------------------------------------------------- */

void FixAdResRegion::copy_arrays(int i, int j, int delflag)
{
  resolution[j] = resolution[i];
  lambda[j] = lambda[i];
}

/* ---------------------------------------------------------------------- */

int FixAdResRegion::pack_exchange(int i, double *buf)
{
  int n = 0;
  buf[n++] = resolution[i];
  buf[n++] = lambda[i];
  return n;
}

/* ---------------------------------------------------------------------- */

int FixAdResRegion::unpack_exchange(int nlocal, double *buf)
{
  int n = 0;
  resolution[nlocal] = (int) buf[n++];
  lambda[nlocal] = buf[n++];
  return n;
}

/* ---------------------------------------------------------------------- */

int FixAdResRegion::pack_restart(int i, double *buf)
{
  int n = 0;
  buf[n++] = resolution[i];
  buf[n++] = lambda[i];
  return n;
}

/* ---------------------------------------------------------------------- */

void FixAdResRegion::unpack_restart(int nlocal, int nth)
{
  double **extra = atom->extra;
  int n = (nth - 1) * 2;
  resolution[nlocal] = (int) extra[nlocal][n++];
  lambda[nlocal] = extra[nlocal][n++];
}

/* ---------------------------------------------------------------------- */

int FixAdResRegion::size_restart(int /*nlocal*/)
{
  return 2;
}

/* ---------------------------------------------------------------------- */

int FixAdResRegion::maxsize_restart()
{
  return 2;
}

/* ---------------------------------------------------------------------- */

double FixAdResRegion::memory_usage()
{
  double bytes = 0.0;
  bytes += (double) maxatom * sizeof(int);      // resolution
  bytes += (double) maxatom * sizeof(double);    // lambda
  return bytes;
}

