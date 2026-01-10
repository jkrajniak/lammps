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
FixStyle(adress/region,FixAdResSRegion);
// clang-format on
#else

#ifndef LMP_FIX_ADRESS_REGION_H
#define LMP_FIX_ADRESS_REGION_H

#include "fix.h"

namespace LAMMPS_NS {

class FixAdResSRegion : public Fix {
 public:
  FixAdResSRegion(class LAMMPS *, int, char **);
  ~FixAdResSRegion() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void post_neighbor() override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;
  int pack_restart(int, double *) override;
  void unpack_restart(int, int) override;
  int size_restart(int) override;
  int maxsize_restart() override;
  double memory_usage() override;

  // Access methods for other AdRes components
  int get_resolution(int) const;
  double get_lambda(int) const;
  double get_lambda_cg(int) const;

 private:
  char *id_region_atomistic;
  char *id_region_transition;
  char *id_region_cg;
  class Region *region_atomistic;
  class Region *region_transition;
  class Region *region_cg;

  int *resolution;        // per-atom resolution state: 0=AT, 1=TR, 2=CG
  double *lambda;         // per-atom lambda value [0,1] for transition region
  int maxatom;

  int transition_axis;    // 0=x, 1=y, 2=z
  double cg_boundary;     // CG region boundary along transition axis
  double at_boundary;     // Atomistic region boundary along transition axis
  double transition_width; // Width of transition region along axis

  enum { RES_ATOMISTIC = 0, RES_TRANSITION = 1, RES_CG = 2 };

  void update_resolution();
  int determine_region(double, double, double);
  double calculate_lambda(double, double, double);
  void calculate_boundaries();
  double get_boundary_coordinate(class Region *reg, int axis);
};

}    // namespace LAMMPS_NS

#endif
#endif

