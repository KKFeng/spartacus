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
/* ----------------------------------------------------------------------
*  SPARTACUS - SPARTA Combined with USp method,
*  a unified stochastic particle (USP) method sovler implement within
*  the framework of SPARTA.
*
*  Developer¡¯s repository:
*  [GitHub]  https://github.com/KKFeng/spartacus
*  or [Gitee]   https://gitee.com/kaikfeng/sparta-usp-para
*
*  Feng Kai, kfeng@buaa.edu.cn
*  Beihang University
------------------------------------------------------------------------- */

#ifdef COMMAND_CLASS

CommandStyle(adapt_grad_compute,AdaptGradCompute)

#else

#ifndef SPARTA_ADAPT_GRAD_COMPUTE_H
#define SPARTA_ADAPT_GRAD_COMPUTE_H

#include "pointers.h"

#define NVALUEMAX 5

namespace SPARTA_NS {

struct AllGrad
{
    cellint id;
    double dt, l;
};
class AdaptGradCompute : protected Pointers {
 public:

  AdaptGradCompute(class SPARTA *);
  ~AdaptGradCompute();
  void command(int, char **);
  void process_args(int, char **);
  void check_args();

 private:
  int me,nprocs;
  int groupbit, style;
  int mod; //max or min or none(default)

  // style = value, for each array, 0, 1, 2~4 = temperature, mass, q_1~3
  int valuewhich_arr[NVALUEMAX], valindex_arr[NVALUEMAX], icompute_arr[NVALUEMAX], ifix_arr[NVALUEMAX];
  char* computeID_arr[NVALUEMAX], * valueID_arr[NVALUEMAX];
  class Compute* compute_arr[NVALUEMAX];
  class Fix* fix_arr[NVALUEMAX];
  // extra buffer, each length = nghost + nlocal
  double* q;
  int* exist_q;

  // if > 0, all cells will get the biggest grad within this range
  // otherwise, they will just get their own grad
  double range;

  //// method
  void compute_grad_value();
  double value_compute(int, int ico);
  double value_fix(int, int ifi);

  void run_comm();
  double cal_grad(int icell);
  void print_grad();

};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
