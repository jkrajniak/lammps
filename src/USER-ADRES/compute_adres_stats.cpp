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

#include "compute_adres_stats.h"

#include "atom.h"
#include "error.h"
#include "fix_adres_region.h"
#include "memory.h"
#include "modify.h"
#include "update.h"

#include <cstring>
#include <mpi.h>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeAdResStats::ComputeAdResStats(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), id_fix_region(nullptr), fix_region(nullptr), nregions(3),
    natoms_region(nullptr), nlocal_region(nullptr)
{
  if (narg < 4) utils::missing_cmd_args(FLERR, "compute adres/stats", error);

  vector_flag = 1;
  size_vector = 3;    // atomistic, transition, CG
  extvector = 0;
  extscalar = 0;

  vector = new double[size_vector];

  id_fix_region = utils::strdup(arg[3]);

  memory->create(natoms_region, nregions, "adres/stats:natoms_region");
  memory->create(nlocal_region, nregions, "adres/stats:nlocal_region");
}

/* ---------------------------------------------------------------------- */

ComputeAdResStats::~ComputeAdResStats()
{
  delete[] id_fix_region;
  memory->destroy(natoms_region);
  memory->destroy(nlocal_region);
  if (!copymode) delete[] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeAdResStats::init()
{
  // get fix adres/region pointer
  fix_region = dynamic_cast<FixAdResRegion *>(modify->get_fix_by_id(id_fix_region));
  if (!fix_region) error->all(FLERR, "Fix {} for compute adres/stats does not exist", id_fix_region);
}

/* ---------------------------------------------------------------------- */

void ComputeAdResStats::compute_vector()
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // zero counts
  for (int i = 0; i < nregions; i++) {
    nlocal_region[i] = 0;
    natoms_region[i] = 0;
  }

  // count atoms in each region
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      int res = fix_region->get_resolution(i);
      if (res >= 0 && res < nregions) {
        nlocal_region[res]++;
      }
    }
  }

  // sum across all processors
  MPI_Allreduce(nlocal_region, natoms_region, nregions, MPI_INT, MPI_SUM, world);

  // store in vector
  vector[0] = natoms_region[0];    // atomistic
  vector[1] = natoms_region[1];    // transition
  vector[2] = natoms_region[2];    // CG
}

