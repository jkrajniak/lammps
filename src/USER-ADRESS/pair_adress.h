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

  // CG particle type for region-based filtering
  int cg_type;  // Atom type of CG particles (0 = not set)

  virtual void allocate();
  double switching_function(double, double, double);
  void interpolate_forces();
  void filter_forces_by_region();  // Filter forces based on region (zeros invalid forces)
  int get_cg_type();  // Get CG type from fix_adress_constraint or return 0
  
  // Methods for filter-during-compute approach
  bool should_compute_AT_force(int i, int j);  // Check if pair should have AT force
  bool should_compute_CG_force(int i, int j);  // Check if pair should have CG force
  void compute_AT_force_pair(int i, int j, int eflag, int vflag);  // Compute AT force for single pair
  void compute_CG_force_pair(int i, int j, int eflag, int vflag);  // Compute CG force for single pair
  void interpolate_energy();  // Interpolate energy based on lambda
  void interpolate_virial();  // Interpolate virial based on lambda
};

}    // namespace LAMMPS_NS

#endif
#endif

