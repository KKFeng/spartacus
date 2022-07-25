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

#ifndef SPARTA_GRID_COMM_MACRO_H
#define SPARTA_GRID_COMM_MACRO_H

#include "pointers.h"

namespace SPARTA_NS {

class GridCommMacro : protected Pointers {
public:

    GridCommMacro(class SPARTA*);
    ~GridCommMacro();
    void runComm();
    void acquire_macro_comm_list_near();
    const CommMacro* interpolation(class Particle::OnePart*);
    
    int nprocs, me;
    class RanPark* random;
    int rand_flag; //init random when first used
    // sending plan
    int nsendproc;
    int * proclist, // size = nsendproc 
        * nsendeachproc, //size = nprocs
        * sizelist;// size = nsendproc
    int ncellsendall, 
        * sendcelllist, // size = ncellsendall
        * sendfirst;    // size = nprocs, don't change it when runComm()!!!

    // receiving plan
    int nrecvproc, nrecvcell, recvsize, * recvicelllist;

    // buffer & Irregular
    char* rbuf, * sbuf;
    class Irregular* irregular;

    // status
    bigint count_sumInter, count_surfInter, count_originInter, count_neighInter,
        count_boundInter, count_outInter, count_warningInter;

private:
    double xnew[3], xhold[3];
    class Particle::OnePart* ipart;
    typedef const CommMacro* (GridCommMacro::* FnPtr)();
    FnPtr interptr;             // ptr to move method
    const CommMacro* interpolation_2d();
    const CommMacro* interpolation_axisym();
    const CommMacro* interpolation_3d();
};


}
#endif