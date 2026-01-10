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

#include "fix_adress_region.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "memory.h"
#include "modify.h"
#include "region.h"
#include "region_block.h"
#include "region_sphere.h"
#include "region_cylinder.h"
#include "update.h"
#include "utils.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixAdResSRegion::FixAdResSRegion(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), id_region_atomistic(nullptr), id_region_transition(nullptr),
    id_region_cg(nullptr), region_atomistic(nullptr), region_transition(nullptr),
    region_cg(nullptr), resolution(nullptr), lambda(nullptr), maxatom(0),
    transition_axis(0), cg_boundary(0.0), at_boundary(0.0), transition_width(0.0)
{
  if (narg < 6) utils::missing_cmd_args(FLERR, "fix adress/region", error);

  restart_peratom = 1;
  restart_global = 0;

  id_region_atomistic = utils::strdup(arg[3]);
  id_region_transition = utils::strdup(arg[4]);
  id_region_cg = utils::strdup(arg[5]);

  // Parse optional axis parameter (default: x-axis)
  transition_axis = 0;  // default to x-axis
  if (narg >= 7) {
    if (strcmp(arg[6], "x") == 0) {
      transition_axis = 0;
    } else if (strcmp(arg[6], "y") == 0) {
      transition_axis = 1;
    } else if (strcmp(arg[6], "z") == 0) {
      transition_axis = 2;
    } else {
      error->all(FLERR, "Invalid axis for fix adress/region: {}. Must be x, y, or z", arg[6]);
    }
  }

  // allocate per-atom arrays
  maxatom = atom->nmax;
  memory->create(resolution, maxatom, "adress/region:resolution");
  memory->create(lambda, maxatom, "adress/region:lambda");

  // initialize arrays
  for (int i = 0; i < maxatom; i++) {
    resolution[i] = RES_CG;    // default to CG
    lambda[i] = 0.0;
  }
}

/* ---------------------------------------------------------------------- */

FixAdResSRegion::~FixAdResSRegion()
{
  delete[] id_region_atomistic;
  delete[] id_region_transition;
  delete[] id_region_cg;
  memory->destroy(resolution);
  memory->destroy(lambda);
}

/* ---------------------------------------------------------------------- */

int FixAdResSRegion::setmask()
{
  int mask = 0;
  mask |= POST_NEIGHBOR;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAdResSRegion::init()
{
  // get region pointers
  region_atomistic = domain->get_region_by_id(id_region_atomistic);
  if (!region_atomistic)
    error->all(FLERR, "Region {} for fix adress/region does not exist", id_region_atomistic);

  region_transition = domain->get_region_by_id(id_region_transition);
  if (!region_transition)
    error->all(FLERR, "Region {} for fix adress/region does not exist", id_region_transition);

  region_cg = domain->get_region_by_id(id_region_cg);
  if (!region_cg)
    error->all(FLERR, "Region {} for fix adress/region does not exist", id_region_cg);

  // Calculate boundaries along transition axis
  calculate_boundaries();
}

/* ---------------------------------------------------------------------- */

void FixAdResSRegion::setup(int /*vflag*/)
{
  // Recalculate boundaries in case regions have changed
  calculate_boundaries();
  update_resolution();
}

/* ---------------------------------------------------------------------- */

void FixAdResSRegion::post_neighbor()
{
  // Recalculate boundaries in case regions are dynamic
  calculate_boundaries();
  update_resolution();
}

/* ---------------------------------------------------------------------- */

void FixAdResSRegion::update_resolution()
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

int FixAdResSRegion::determine_region(double x, double y, double z)
{
  // Get position along transition axis
  double pos = (transition_axis == 0) ? x : (transition_axis == 1) ? y : z;

  // Use axis-based detection: check position relative to boundaries
  if (pos < cg_boundary) {
    // Before CG boundary - check if in CG region
    if (region_cg && region_cg->match(x, y, z)) {
      return RES_CG;
    }
  } else if (pos > at_boundary) {
    // After atomistic boundary - check if in atomistic region
    if (region_atomistic && region_atomistic->match(x, y, z)) {
      return RES_ATOMISTIC;
    }
  } else {
    // In transition region - check if actually in transition region
    if (region_transition && region_transition->match(x, y, z)) {
      return RES_TRANSITION;
    }
  }

  // Fallback: check regions in order (for cases where regions don't align with axis)
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

double FixAdResSRegion::calculate_lambda(double x, double y, double z)
{
  // Get position along transition axis
  double pos = (transition_axis == 0) ? x : (transition_axis == 1) ? y : z;

  // If outside transition region boundaries, return fixed values
  if (pos <= cg_boundary) {
    return 0.0;  // CG region
  }
  if (pos >= at_boundary) {
    return 1.0;  // Atomistic region
  }

  // Calculate normalized distance from CG boundary
  // d = (pos - cg_boundary) / transition_width
  double d = (pos - cg_boundary) / transition_width;

  // Clamp to [0, 1] range (safety check)
  d = fmax(0.0, fmin(1.0, d));

  // Square for smooth transition: λ = d²
  return d * d;
}

/* ---------------------------------------------------------------------- */

int FixAdResSRegion::get_resolution(int i) const
{
  if (i < 0 || i >= atom->nmax) return RES_CG;
  return resolution[i];
}

/* ---------------------------------------------------------------------- */

double FixAdResSRegion::get_lambda(int i) const
{
  if (i < 0 || i >= atom->nmax) return 0.0;
  return lambda[i];
}

/* ---------------------------------------------------------------------- */

double FixAdResSRegion::get_lambda_cg(int i) const
{
  if (i < 0 || i >= atom->nmax) return 0.0;
  
  double lambda_at = get_lambda(i);
  
  // Get CG type from fix adress/constraint if available
  int cg_type = 0;
  if (modify) {
    for (int j = 0; j < modify->nfix; j++) {
      Fix *fix = modify->fix[j];
      if (strcmp(fix->style, "adress/constraint") == 0) {
        int dim = 0;
        int *cg_type_ptr = (int *) fix->extract("cg_type", dim);
        if (cg_type_ptr) {
          cg_type = *cg_type_ptr;
          break;
        }
      }
    }
  }
  
  // If this is a CG particle, invert lambda: λ_CG = 1 - λ_AT
  if (cg_type > 0 && atom->type[i] == cg_type) {
    return 1.0 - lambda_at;
  }
  
  // Atomistic particles use normal lambda
  return lambda_at;
}

/* ---------------------------------------------------------------------- */

void FixAdResSRegion::grow_arrays(int nmax)
{
  memory->grow(resolution, nmax, "adress/region:resolution");
  memory->grow(lambda, nmax, "adress/region:lambda");
  maxatom = nmax;
}

/* ---------------------------------------------------------------------- */

void FixAdResSRegion::copy_arrays(int i, int j, int delflag)
{
  resolution[j] = resolution[i];
  lambda[j] = lambda[i];
}

/* ---------------------------------------------------------------------- */

int FixAdResSRegion::pack_exchange(int i, double *buf)
{
  int n = 0;
  buf[n++] = resolution[i];
  buf[n++] = lambda[i];
  return n;
}

/* ---------------------------------------------------------------------- */

int FixAdResSRegion::unpack_exchange(int nlocal, double *buf)
{
  int n = 0;
  resolution[nlocal] = (int) buf[n++];
  lambda[nlocal] = buf[n++];
  return n;
}

/* ---------------------------------------------------------------------- */

int FixAdResSRegion::pack_restart(int i, double *buf)
{
  int n = 0;
  buf[n++] = resolution[i];
  buf[n++] = lambda[i];
  return n;
}

/* ---------------------------------------------------------------------- */

void FixAdResSRegion::unpack_restart(int nlocal, int nth)
{
  double **extra = atom->extra;
  int n = (nth - 1) * 2;
  resolution[nlocal] = (int) extra[nlocal][n++];
  lambda[nlocal] = extra[nlocal][n++];
}

/* ---------------------------------------------------------------------- */

int FixAdResSRegion::size_restart(int /*nlocal*/)
{
  return 2;
}

/* ---------------------------------------------------------------------- */

int FixAdResSRegion::maxsize_restart()
{
  return 2;
}

/* ---------------------------------------------------------------------- */

double FixAdResSRegion::memory_usage()
{
  double bytes = 0.0;
  bytes += (double) maxatom * sizeof(int);      // resolution
  bytes += (double) maxatom * sizeof(double);    // lambda
  return bytes;
}

/* ---------------------------------------------------------------------- */

void FixAdResSRegion::calculate_boundaries()
{
  // Calculate boundaries along transition axis for each region
  double cg_lo, cg_hi, at_lo, at_hi;

  // Get CG region boundaries
  cg_lo = get_boundary_coordinate(region_cg, transition_axis);
  cg_hi = get_boundary_coordinate(region_cg, transition_axis + 3);  // +3 for hi boundaries

  // Get atomistic region boundaries
  at_lo = get_boundary_coordinate(region_atomistic, transition_axis);
  at_hi = get_boundary_coordinate(region_atomistic, transition_axis + 3);

  // Determine which boundaries to use based on region layout
  // Assume CG region comes before atomistic region along transition axis
  // CG boundary is where CG region ends (closest to atomistic)
  // Atomistic boundary is where atomistic region starts (closest to CG)
  
  if (cg_hi < at_lo) {
    // Normal layout: CG -> Transition -> Atomistic
    cg_boundary = cg_hi;
    at_boundary = at_lo;
  } else if (at_hi < cg_lo) {
    // Reversed layout: Atomistic -> Transition -> CG
    cg_boundary = cg_lo;
    at_boundary = at_hi;
  } else {
    // Regions overlap or are misaligned - use extent values
    cg_boundary = (transition_axis == 0) ? region_cg->extent_xhi :
                   (transition_axis == 1) ? region_cg->extent_yhi : region_cg->extent_zhi;
    at_boundary = (transition_axis == 0) ? region_atomistic->extent_xlo :
                   (transition_axis == 1) ? region_atomistic->extent_ylo : region_atomistic->extent_zlo;
  }

  // Calculate transition width
  transition_width = fabs(at_boundary - cg_boundary);

  // Validate
  if (transition_width <= 0.0) {
    error->warning(FLERR, "Fix adress/region: Transition width is zero or negative. "
                   "Regions may not be properly aligned along transition axis.");
  }
}

/* ---------------------------------------------------------------------- */

double FixAdResSRegion::get_boundary_coordinate(class Region *reg, int boundary_type)
{
  // boundary_type: 0=xlo, 1=ylo, 2=zlo, 3=xhi, 4=yhi, 5=zhi
  if (!reg) return 0.0;

  // Try to get boundary from region extent (works for all region types)
  if (boundary_type == 0) return reg->extent_xlo;
  if (boundary_type == 1) return reg->extent_ylo;
  if (boundary_type == 2) return reg->extent_zlo;
  if (boundary_type == 3) return reg->extent_xhi;
  if (boundary_type == 4) return reg->extent_yhi;
  if (boundary_type == 5) return reg->extent_zhi;

  return 0.0;
}

