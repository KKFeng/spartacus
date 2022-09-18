/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#ifdef COMMAND_CLASS

CommandStyle(adapt_dt_weight,AdaptDtWeight)

#else

#ifndef SPARTA_ADAPT_DT_WEIGHT_H
#define SPARTA_ADAPT_DT_WEIGHT_H

#include "pointers.h"
#include "hash3.h"

#define NVALUEMAX 5

namespace SPARTA_NS {

class AdaptDtWeight : protected Pointers {
 public:


  AdaptDtWeight(class SPARTA *);
  ~AdaptDtWeight();
  void command(int, char **);
  void process_args(int, char **);
  void check_args();

 private:
  int me,nprocs;
  int groupbit, style;
  int mod; //max or min or none(default)

  // style = surf or near_surf
  int sgroupbit, surf_ndt;
  double surf_dist; // for style = near_surf. distance to cell center
  struct MyRegion
  {
	  double x[3], radius;
  };

  MyRegion *regionlist;
  int nregion, maxregion;

  // style = value
  int valuewhich, valindex,icompute,ifix;
  char* computeID, * valueID;
  class Compute* compute;
  class Fix* fix;
  double thresh;
  int max_dt; //also used in style = value_heatflux

  // style = value_heatflux, for each arrary, 0, 1, 2~4=temperature, mass, heat flux_x~z
  int valuewhich_arr[NVALUEMAX], valindex_arr[NVALUEMAX], icompute_arr[NVALUEMAX], ifix_arr[NVALUEMAX];
  int nvalue;
  double coef;
  char* computeID_arr[NVALUEMAX], * valueID_arr[NVALUEMAX];
  class Compute* compute_arr[NVALUEMAX];
  class Fix* fix_arr[NVALUEMAX];

  //style = same
  int same_dt;

  // method

  void set_weight_surf();
  void set_weight_nearsurf();
  void set_weight_value();
  void set_weight_value_heatflux();
  void set_weight_same();
  double value_compute(int, int ico = -1);
  double value_fix(int, int ifi = -1);

  void add_region(double*, double);
  bool in_region(MyRegion&, double*);
  void gather_allregion();

};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
