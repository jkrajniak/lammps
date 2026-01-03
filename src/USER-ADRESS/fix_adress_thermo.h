/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(adress/thermo,FixAdResSThermo);
// clang-format on
#else

#ifndef LMP_FIX_ADRESS_THERMO_H
#define LMP_FIX_ADRESS_THERMO_H

#include "fix.h"

namespace LAMMPS_NS {

class FixAdResSThermo : public Fix {
 public:
  FixAdResSThermo(class LAMMPS *, int, char **);
  ~FixAdResSThermo() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  double memory_usage() override;

 private:
  char *id_region;
  char *id_fix_region;
  class Region *region;
  class FixAdResSRegion *fix_region;

  double transition_width;
  double kT;                    // temperature in energy units
  double **thermo_force;        // per-atom thermodynamic force
  int maxatom;

  void calculate_thermodynamic_force();
  double calculate_density_gradient(int);
};

}    // namespace LAMMPS_NS

#endif
#endif

