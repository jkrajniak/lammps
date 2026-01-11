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
FixStyle(adress/constraint,FixAdResSConstraint);
// clang-format on
#else

#ifndef LMP_FIX_ADRESS_CONSTRAINT_H
#define LMP_FIX_ADRESS_CONSTRAINT_H

#include "fix.h"

namespace LAMMPS_NS {

class FixAdResSConstraint : public Fix {
 public:
  FixAdResSConstraint(class LAMMPS *, int, char **);
  ~FixAdResSConstraint() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void initial_integrate(int) override;
  void post_force(int) override;
  void final_integrate() override;
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
  void *extract(const char *, int &) override;

 private:
  char *id_fix_region;              // ID of fix adress/region
  class FixAdResSRegion *fix_region; // Pointer to fix adress/region
  int cg_type;                       // Atom type of CG particles (e.g., 4 for WCG)

  // Molecule mapping
  int *molecule_map;                 // molecule_map[atom_index] = molecule_index (-1 if not mapped)
  int *cg_atom_index;                // cg_atom_index[molecule_index] = CG atom local index (-1 if none)
  tagint *mol_id_list;                // mol_id_list[molecule_index] = molecule ID
  int nmolecules;                    // Number of unique molecules
  int maxmolecule;                   // Maximum number of molecules (for array allocation)

  // COM and CG particle data
  double **x_com;                    // COM positions [nmolecules][3]
  double **v_com;                    // COM velocities [nmolecules][3]
  double **f_com;                    // COM forces [nmolecules][3]
  double **x_cg_stored;              // Stored CG particle positions [nmolecules][3]
  double **v_cg_stored;              // Stored CG particle velocities [nmolecules][3]
  double *mass_com;                  // Total mass of each molecule [nmolecules]

  // Relative positions (for CG region constraints)
  double **displace;                 // Relative positions in body frame [natoms][3]
  int maxatom;                       // Maximum number of atoms

  // Lambda thresholds
  double lambda_cg_threshold;        // Below this: CG region (default 0.1)
  double lambda_at_threshold;        // Above this: Atomistic region (default 0.9)

  void map_molecules();
  void calculate_com_from_atoms(int mol_idx);
  void calculate_com_from_cg(int mol_idx);
  void constrain_atoms_to_com(int mol_idx);
  void update_cg_from_com(int mol_idx);
  int find_cg_particle_for_molecule(tagint mol_id);
  bool is_molecule_in_cg_region(int mol_idx);
};

}    // namespace LAMMPS_NS

#endif
#endif

