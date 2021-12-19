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

#include "string.h"
#include "grid.h"
#include "geometry.h"
#include "domain.h"
#include "region.h"
#include "surf.h"
#include "comm.h"
#include "modify.h"
#include "fix.h"
#include "compute.h"
#include "output.h"
#include "dump.h"
#include "irregular.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "grid_comm_macro.h"
#include "collide.h"
#include "particle.h"
#include "random_park.h"
#include "update.h"
#include "domain.h"

using namespace SPARTA_NS;
using namespace MathConst;

#define DELTA 8192
#define DELTAPARENT 1024
#define LARGE 256000
#define BIG 1.0e20
#define MAXGROUP 32
#define MAXLEVEL 32

// default values, can be overridden by global command

#define MAXSURFPERCELL  100
#define MAXSPLITPERCELL 10

//enum { XLO, XHI, YLO, YHI, ZLO, ZHI, INTERIOR };         // same as Domain
enum { PERIODIC, OUTFLOW, REFLECT, SURFACE, AXISYM };  // same as Domain
//enum { REGION_ALL, REGION_ONE, REGION_CENTER };      // same as Surf
//enum { PERAUTO, PERCELL, PERSURF };                  // several files
//
//// cell type = OUTSIDE/INSIDE/OVERLAP if entirely outside/inside surfs
////   or has any overlap with surfs including grazing or touching
//// corner point = OUTSIDE/INSIDE (not OVERLAP) if outside/inside
////   if exactly on surf, is still marked OUTSIDE or INSIDE by cut2d/cut3d
////   corner pts are only set if cell type = OVERLAP
//
//enum { UNKNOWN, OUTSIDE, INSIDE, OVERLAP };           // several files
//enum { NCHILD, NPARENT, NUNKNOWN, NPBCHILD, NPBPARENT, NPBUNKNOWN, NBOUND };  // Update
//enum { NOWEIGHT, VOLWEIGHT, RADWEIGHT };

GridCommMacro::GridCommMacro(SPARTA* sparta) : Pointers(sparta) {
    me = comm->me;
    nprocs = comm->nprocs;
    random = new RanPark(update->ranmaster->uniform());
    nsendproc = 0;
    proclist = new int[nprocs];
    nsendeachproc = new int[nprocs];
    sizelist = new int[nprocs];
    sendfirst = new int[nprocs];
    sendcelllist = nullptr;
    ncellsendall = 0;

    recvsize = 0;
    nrecvproc = 0;
    nrecvcell = 0;
    recvicelllist = nullptr;

    rbuf = nullptr;
    sbuf = nullptr;
    irregular = new Irregular(sparta);

}

GridCommMacro::~GridCommMacro(){
    delete[] proclist;
    delete[] nsendeachproc;
    delete[] sizelist;
    delete[] sendfirst;
    delete[] sendcelllist;
    delete[] recvicelllist;
    delete[] rbuf;
    delete[] sbuf;
    delete irregular;
    delete random;

}


void GridCommMacro::acquire_macro_comm_list_near()
{
    //Currently only consider global gridcut = -1.0, thus I have all child cell information
    // of the whole sim box.

    if (grid->cutoff >= 0.0) error->one(FLERR, "Macro Comm: cutoff >=0.0");
    if (!grid->exist_ghost) error->one(FLERR, "Macro Comm: Ghost cell not exist");

    // bb lo/hi = bounding box for my owned cells

    int i;
    double bblo[3], bbhi[3];
    double* lo, * hi;

    for (i = 0; i < 3; i++) {
        bblo[i] = BIG;
        bbhi[i] = -BIG;
    }
    auto& cells = grid->cells;
    auto& nlocal = grid->nlocal;
    for (int icell = 0; icell < nlocal; icell++) {
        if (cells[icell].nsplit <= 0) continue;
        lo = cells[icell].lo;
        hi = cells[icell].hi;
        for (i = 0; i < 3; i++) {
            bblo[i] = MIN(bblo[i], lo[i]);
            bbhi[i] = MAX(bbhi[i], hi[i]);
        }
    }

    // cut = max side length of all child cells in this proc

    double cut = 0.0;
    int maxChildLevel = 1000;
    for (int icell = 0; icell < nlocal; icell++) {
        if (cells[icell].nsplit > 0 && cells[icell].level < maxChildLevel) {
            maxChildLevel = cells[icell].level;
            for (i = 0; i < domain->dimension; ++i) {
                cut = MAX(cut, cells[icell].hi[i] - cells[icell].lo[i]);
            }
        }
    }


    // ebb lo/hi = bbox + cut
    // trim to simulation box in non-periodic dims
    // 
    // -----Warnning! condition below is not guaranteed ------------
    // if bblo/hi is at periodic boundary and cutoff is 0.0,
    //   add cell_epsilon to insure ghosts across periodic boundary acquired,
    //   else may be UNKNOWN to owned cell
    //--------------------------------------------------------------

    double* boxlo = domain->boxlo;
    double* boxhi = domain->boxhi;
    int* bflag = domain->bflag;

    double ebblo[3], ebbhi[3];
    cut *= 0.99999999;
    if(domain->dimension == 2) {
        ebblo[2] = -0.5; ebbhi[2] = 0.5;
    }
    for (i = 0; i < domain->dimension; i++) {
        ebblo[i] = bblo[i] - cut;
        ebbhi[i] = bbhi[i] + cut;
        if (bflag[2 * i] != PERIODIC) ebblo[i] = MAX(ebblo[i], boxlo[i]);
        if (bflag[2 * i] != PERIODIC) ebbhi[i] = MIN(ebbhi[i], boxhi[i]);
    }

    // box = ebbox split across periodic BC
    // 27 is max number of periodic images in 3d

    Grid::Box box[27];
    int nbox = grid->box_periodic(ebblo, ebbhi, box);

    // boxall = collection of boxes from all procs

    int me = comm->me;
    int nprocs = comm->nprocs;

    int nboxall;
    MPI_Allreduce(&nbox, &nboxall, 1, MPI_INT, MPI_SUM, world);

    int* recvcounts, * displs;
    memory->create(recvcounts, nprocs, "grid:recvcounts");
    memory->create(displs, nprocs, "grid:displs");

    int nsend = nbox * sizeof(Grid::Box);
    MPI_Allgather(&nsend, 1, MPI_INT, recvcounts, 1, MPI_INT, world);
    displs[0] = 0;
    for (i = 1; i < nprocs; i++) displs[i] = displs[i - 1] + recvcounts[i - 1];

    Grid::Box* boxall = new Grid::Box[nboxall];
    MPI_Allgatherv(box, nsend, MPI_CHAR, boxall, recvcounts, displs, MPI_CHAR, world);

    memory->destroy(recvcounts);
    memory->destroy(displs);

    // nlist = # of boxes that overlap with my bbox, skipping self boxes
    // list = indices into boxall of overlaps
    // overlap = true overlap or just touching

    int nlist = 0;
    int* list;
    memory->create(list, nboxall, "grid:list");

    for (i = 0; i < nboxall; i++) {
        if (boxall[i].proc == me) continue;
        if (grid->box_overlap(bblo, bbhi, boxall[i].lo, boxall[i].hi)) list[nlist++] = i;
    }

    // loop over my owned cells, not including sub cells
    // each may overlap with multiple boxes in list
    // on 1st pass, just tally memory to send copies of my cells
    // use lastproc to insure a cell only overlaps once per other proc
    // if oflag = 2 = my cell just touches box,
    // so flag grid cell as EMPTY ghost by setting nsurf = -1

    int j, oflag, lastproc, nsurf_hold;

    nsend = 0;
    memset(nsendeachproc, 0, sizeof(int) * nprocs);
    for (int icell = 0; icell < nlocal; icell++) {
        if (cells[icell].nsplit <= 0) continue;
        lo = cells[icell].lo;
        hi = cells[icell].hi;
        lastproc = -1;
        for (i = 0; i < nlist; i++) {
            j = list[i];
            oflag = grid->box_overlap(lo, hi, boxall[j].lo, boxall[j].hi);
            if (oflag != 1) continue;
            if (boxall[j].proc == lastproc) continue;
            lastproc = boxall[j].proc;
            nsendeachproc[lastproc] += 1;
            nsend++;
        }
    }

    // set gridCommMacro 
    nsendproc = 0;
    for (int i = 0; i < nprocs; i++) {
        if (nsendeachproc[i]) {
            proclist[nsendproc] = i;
            nsendproc++;
        }
    }
    for (int i = 0; i < nsendproc; i++) {
        sizelist[i] = nsendeachproc[proclist[i]] * sizeof(CommMacro);
    }
    ncellsendall = nsend;
    memory->destroy(sendcelllist);
    memory->create(sendcelllist, ncellsendall, "gridCommMacro : sendcelllist");
    memory->destroy(sbuf);
    memory->create(sbuf, ncellsendall * sizeof(CommMacro), "gridCommMacro : sbuf");


    int* sf = new int[nprocs]; //sendfirst
    sf[0] = 0;
    for (int i = 1; i < nprocs; ++i) {
        sf[i] = sf[i - 1] + nsendeachproc[i - 1];
    }

    memcpy(sendfirst, sf, sizeof(int)* nprocs);

    for (int icell = 0; icell < nlocal; icell++) {
        if (cells[icell].nsplit <= 0) continue;
        lo = cells[icell].lo;
        hi = cells[icell].hi;
        lastproc = -1;
        for (i = 0; i < nlist; i++) {
            j = list[i];
            oflag = grid->box_overlap(lo, hi, boxall[j].lo, boxall[j].hi);
            if (oflag != 1) continue;
            if (boxall[j].proc == lastproc) continue;
            lastproc = boxall[j].proc;
            sendcelllist[sf[lastproc]++] = icell;
        }
    }
    
    //if (me == 0) {
    //    fprintf(screen, "ncellsendall: %d\n", ncellsendall);
    //    fflush(screen);
    //    fprintf(screen, "sendcelllist: %d,%d,%d,%d,%d,%d\n",
    //        sendcelllist[0], sendcelllist[1], sendcelllist[2], 
    //        sendcelllist[3], sendcelllist[4], sendcelllist[5] );
    //    fflush(screen);
    //}


    for (int i = 0; i < ncellsendall; ++i) {
        if (me == 0)
        fprintf(screen, "cellId: %d\n", grid->cells[sendcelllist[i]].id);
        memcpy(sbuf + i * sizeof(CommMacro), &(grid->cells[sendcelllist[i]].id), sizeof(cellint));
    }

    //DEBUG
    for (int i = 0; i < nprocs - 1; ++i) {
        if (sf[i] != sendfirst[i + 1])
            error->one(FLERR, "sendcelllist set error!");
    }

    delete[] sf;
    // clean up

    memory->destroy(list);
    delete[] boxall;

    // perform irregular communication of list of neigh cells 
    // whose V&T needed to be transfer.

    nrecvproc = irregular->create_data_variable(nsendproc, proclist, sizelist,
        recvsize, 1); //must sort
    nrecvcell = recvsize / sizeof(CommMacro);
    memory->create(rbuf, recvsize, "gridCommMacro:rbuf");
    memset(rbuf, 0, recvsize);

    irregular->exchange_variable(sbuf, sizelist, rbuf);
    memory->destroy(recvicelllist);
    memory->create(recvicelllist, nrecvcell,"GridCommMacro:recvicellist");
    for (int i = 0; i < nrecvcell; ++i) {
        cellint id = 0;
        memcpy(&id, rbuf + i * sizeof(CommMacro), sizeof(cellint));
        if (grid->hash->find(id) != grid->hash->end()) {
            recvicelllist[i] = (*(grid->hash))[id];
        }
        else {
            error->one(FLERR, "GridCommMacro : no such owned or ghost cell");
        }
        //DEBUG
        if (me == 0)  fprintf(screen, "cellId: %d, local id: %d\n", id, recvicelllist[i]);           
    }

}

void GridCommMacro::runComm() 
{
    // pack macro, preparing for comm
    for (int i = 0; i < ncellsendall; ++i) {
        memcpy(sbuf + i * sizeof(CommMacro), 
            &(grid->cells[sendcelllist[i]].macro), sizeof(CommMacro));
    }

    irregular->exchange_variable(sbuf, sizelist, rbuf);


    for (int i = 0; i < nrecvcell; ++i) {
        memcpy(&(grid->cells[recvicelllist[i]].macro),
            rbuf + i * sizeof(CommMacro), sizeof(CommMacro));
    }
}

int GridCommMacro::interpolation(Particle::OnePart* ipart)
{
    double x[3];
    double* lo = domain->boxlo;
    double* hi = domain->boxhi;
    for (int i = 0; i < 3; ++i) {
        x[i] = ipart->x[i] + (random->uniform() - 0.5) *
            (grid->cells[ipart->icell].hi[i] - grid->cells[ipart->icell].lo[i]);;
    }
    if (x[0] < lo[0] || x[0] > hi[0] ||
        x[1] < lo[1] || x[1] > hi[1] || 
        x[2] < lo[2] || x[2] > hi[2] ) return ipart->icell;
    int id = grid->id_find_child(0, 0, domain->boxlo, domain->boxlo, x);
    if (id == -1) id = ipart->icell;
    return id;
}