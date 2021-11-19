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

#ifndef SPARTA_GRID_COMM_MACRO_H
#define SPARTA_GRID_COMM_MACRO_H

//#include "stdio.h"
#include "pointers.h"
//#include "hash3.h"
//#include "my_page.h"
//#include "surf.h"

namespace SPARTA_NS {

class GridCommMacro : protected Pointers {
public:
    GridCommMacro(class SPARTA*);
    ~GridCommMacro();
    void runComm();
    void acquire_macro_comm_list_near();
    
    int nprocs, me;
        
    // sending plan
    int nsendproc;
    int * proclist, // size = nsendproc 
        * nsendeachproc, //size = nprocs
        * sizelist;// size = nsendproc
    int ncellsendall, 
        * sendcelllist, // size = ncellsendall
        * sendfirst;    // size = nprocs, don't change it each runComm()

    // receiving plan
    int nrecvproc, nrecvcell, recvsize, * recvicelllist;

    // buffer & Irregular
    char* rbuf, * sbuf;
    class Irregular* irregular;






};


}
#endif