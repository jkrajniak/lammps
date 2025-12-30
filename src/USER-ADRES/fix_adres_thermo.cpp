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

#include "fix_adres_thermo.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "fix_adres_region.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "region.h"
#include "respa.h"
#include "update.h"
#include "utils.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixAdResThermo::FixAdResThermo(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), id_region(nullptr), id_fix_region(nullptr), region(nullptr),
    fix_region(nullptr), transition_width(0.0), kT(0.0), thermo_force(nullptr), maxatom(0)
{
  if (narg < 5) utils::missing_cmd_args(FLERR, "fix adres/thermo", error);

  id_region = utils::strdup(arg[3]);
  transition_width = utils::numeric(FLERR, arg[4], false, lmp);

  if (transition_width <= 0.0)
    error->all(FLERR, "Fix adres/thermo transition-width must be > 0.0");

  // optional fix adres/region ID
  int iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "fix") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix adres/thermo fix", error);
      id_fix_region = utils::strdup(arg[iarg + 1]);
      iarg += 2;
    } else
      error->all(FLERR, "Unknown fix adres/thermo keyword: {}", arg[iarg]);
  }

  // allocate per-atom array for thermodynamic force
  maxatom = atom->nmax;
  memory->create(thermo_force, maxatom, 3, "adres/thermo:thermo_force");
}

/* ---------------------------------------------------------------------- */

FixAdResThermo::~FixAdResThermo()
{
  delete[] id_region;
  delete[] id_fix_region;
  memory->destroy(thermo_force);
}

/* ---------------------------------------------------------------------- */

int FixAdResThermo::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAdResThermo::init()
{
  // get region pointer
  region = domain->get_region_by_id(id_region);
  if (!region) error->all(FLERR, "Region {} for fix adres/thermo does not exist", id_region);

  // get fix adres/region pointer if specified
  if (id_fix_region) {
    fix_region = dynamic_cast<FixAdResRegion *>(modify->get_fix_by_id(id_fix_region));
    if (!fix_region)
      error->all(FLERR, "Fix {} for fix adres/thermo does not exist", id_fix_region);
  }

  // get temperature from compute
  if (force->boltz == 0.0)
    error->all(FLERR, "Fix adres/thermo requires temperature to be defined");
}

/* ---------------------------------------------------------------------- */

void FixAdResThermo::setup(int vflag)
{
  if (utils::strmatch(update->integrate_style, "^verlet"))
    post_force(vflag);
  else {
    int nlevels_respa = (dynamic_cast<Respa *>(update->integrate))->nlevels;
    for (int ilevel = 0; ilevel < nlevels_respa; ilevel++) {
      (dynamic_cast<Respa *>(update->integrate))->copy_flevel_f(ilevel);
      post_force_respa(vflag, ilevel, 0);
      (dynamic_cast<Respa *>(update->integrate))->copy_f_flevel(ilevel);
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixAdResThermo::post_force(int /*vflag*/)
{
  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // update region before matching
  if (region) region->prematch();

  // reallocate arrays if necessary
  if (atom->nmax > maxatom) {
    maxatom = atom->nmax;
    memory->destroy(thermo_force);
    memory->create(thermo_force, maxatom, 3, "adres/thermo:thermo_force");
  }

  // zero thermodynamic force array
  for (int i = 0; i < nlocal; i++) {
    thermo_force[i][0] = 0.0;
    thermo_force[i][1] = 0.0;
    thermo_force[i][2] = 0.0;
  }

  // calculate thermodynamic force
  calculate_thermodynamic_force();

  // apply thermodynamic force correction
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (region && region->match(x[i][0], x[i][1], x[i][2])) {
        f[i][0] += thermo_force[i][0];
        f[i][1] += thermo_force[i][1];
        f[i][2] += thermo_force[i][2];
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixAdResThermo::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == 0) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixAdResThermo::calculate_thermodynamic_force()
{
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // Get temperature from thermo compute
  // This is a simplified implementation
  // In a full implementation, we would calculate density gradients
  // and apply the thermodynamic force: F_th = -kT * grad(ln(rho))

  if (!fix_region) return;

  // For each atom in transition region, calculate thermodynamic force
  // based on density gradient
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (region && region->match(x[i][0], x[i][1], x[i][2])) {
        double lambda = fix_region->get_lambda(i);
        if (lambda > 0.0 && lambda < 1.0) {
          // Calculate density gradient (simplified)
          double grad_rho = calculate_density_gradient(i);

          // Thermodynamic force: F_th = -kT * grad(ln(rho)) = -kT/rho * grad(rho)
          // Simplified: assume constant density, use lambda gradient
          double lambda_grad = (1.0 - lambda) / transition_width;
          double force_mag = force->boltz * force->boltz * lambda_grad;

          // Apply force in direction of gradient (simplified: x-direction)
          // In full implementation, would use actual gradient direction
          thermo_force[i][0] = force_mag;
          thermo_force[i][1] = 0.0;
          thermo_force[i][2] = 0.0;
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

double FixAdResThermo::calculate_density_gradient(int i)
{
  // Placeholder for density gradient calculation
  // In full implementation, would calculate local density
  // and its gradient using neighbor lists
  return 1.0;
}

/* ---------------------------------------------------------------------- */

double FixAdResThermo::memory_usage()
{
  double bytes = 0.0;
  bytes += (double) maxatom * 3 * sizeof(double);    // thermo_force
  return bytes;
}

