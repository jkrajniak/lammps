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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(adress,PairAdResS);
// clang-format on
#else

#ifndef LMP_PAIR_ADRESS_H
#define LMP_PAIR_ADRESS_H

#include "pair.h"

namespace LAMMPS_NS {

class PairAdResS : public Pair {
 public:
  PairAdResS(class LAMMPS *);
  ~PairAdResS() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;
  double single(int, int, int, int, double, double, double, double &) override;
  void *extract(const char *, int &) override;

 protected:
  double cut_global;
  double **cut;
  char *id_fix_region;
  class FixAdResSRegion *fix_region;

  // Sub-styles for atomistic and CG interactions
  class Pair *pair_atomistic;
  class Pair *pair_cg;
  char *style_atomistic;
  char *style_cg;

  // Temporary force arrays for interpolation
  double **f_atomistic;
  double **f_cg;
  int nmax_force;

  virtual void allocate();
  double switching_function(double, double, double);
  void interpolate_forces();
};

}    // namespace LAMMPS_NS

#endif
#endif

