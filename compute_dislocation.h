/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2015) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Jianfeng Huang Glasgow Caledonian University
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(dislocation,ComputeDislocation)

#else

#ifndef LMP_COMPUTE_DISLOCATION_H
#define LMP_COMPUTE_DISLOCATION_H

#include "compute.h"
#include <vector>
using namespace std;

namespace LAMMPS_NS {

class ComputeDislocation : public Compute {
 public:
  ComputeDislocation(class LAMMPS *, int, char **);
  ~ComputeDislocation();
  void init();
  void init_list(int, class NeighList *);
  void compute_peratom();
  double memory_usage();
  int add_atom_to_dislocation(tagint, class Dislocation *);
  void add_atoms_to_dislocation(int, tagint*, class Dislocation *);
  int is_atom_next_to_dislocation(tagint, Dislocation *);
  void find_neighbor_dislocation_atom(int seedAtom,
	  std::vector<tagint>&,
      std::vector<tagint>&);
  void flat_dislocation_atoms(std::vector<tagint>&);
  void find_nearest_neighbor(tagint, std::vector<tagint>&);
  void Burgers(class Dislocation*);
  void Length(class Dislocation*);
  void CalcLengthAtoms(tagint ,
	  std::vector<tagint> ,
	  std::vector<tagint>& ,
	  class Dislocation* );
 private:
  int nmax;
  double cutsq;
  class NeighList *list;
  int **nearest;
  int *nnearest;
  double *pattern;
  double **atomSlipVector;
  int numDislocation;
  double **displacements;
  class Dislocation **dislocations;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute cna/atom requires a pair style be defined

Self-explantory.

E: Compute cna/atom cutoff is longer than pairwise cutoff

Self-explantory.

W: Compute cna/atom cutoff may be too large to find ghost atom neighbors

The neighbor cutoff used may not encompass enough ghost atoms
to perform this operation correctly.

W: More than one compute cna/atom defined

It is not efficient to use compute cna/atom  more than once.

W: Too many neighbors in CNA for %d atoms

More than the maximum # of neighbors was found multiple times.  This
was unexpected.

W: Too many common neighbors in CNA %d times

More than the maximum # of neighbors was found multiple times.  This
was unexpected.

*/
