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
// DEBUG
#include "update.h"

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
            for (i = 0; i < 3; ++i) {
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
    for (i = 0; i < 3; i++) {
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
    int sendsize = 0;
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
            //if (oflag == 2) {
            //    nsurf_hold = cells[icell].nsurf;
            //    cells[icell].nsurf = -1;
            //}
            sendsize += grid->pack_one(icell, NULL, 0, 0, 0, 0);
            //if (oflag == 2) cells[icell].nsurf = nsurf_hold;
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
    int sz = sizeof(CommMacro);
    for (int i = 0; i < ncellsendall; ++i) {
        memcpy(sbuf + i * sizeof(CommMacro), &grid->cells[sendcelllist[i]].id, sizeof(cellint));
    }

    int* sf = new int[nprocs]; //sendfirst
    for (int i = 1; i < nprocs; ++i) {
        sf[i] = sf[i - 1] + nsendeachproc[i - 1];
    }
    memcpy(sendfirst, sf, sizeof(int)* nprocs);
    // 2021年11月18日22:15:47 以下为需要更改的位置

    //char* sbuf;
    //memory->create(sbuf, sendsize, "grid:sbuf");
    //memset(sbuf, 0, sendsize);

    //int* proclist, * sizelist;
    //memory->create(proclist, nsend, "grid:proclist");
    //memory->create(sizelist, nsend, "grid:sizelist");

    // on 2nd pass over local cells, fill the send buf
    // use lastproc to insure a cell only overlaps once per other proc
    // if oflag = 2 = my cell just touches box,
    // so flag grid cell as EMPTY ghost by setting nsurf = -1

    //nsend = 0;
    //sendsize = 0;
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
            //if (oflag == 2) {
            //    nsurf_hold = cells[icell].nsurf;
            //    cells[icell].nsurf = -1;
            //}
            //sizelist[nsend] = grid->pack_one(icell, &sbuf[sendsize], 0, 0, 0, 1);
            //if (oflag == 2) cells[icell].nsurf = nsurf_hold;
            //proclist[nsend] = lastproc;
            //sendsize += sizelist[nsend];
            //nsend++;
        }
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

    // perform irregular communication of list of ghost cells

    //Irregular* irregular = new Irregular(sparta);
    //int recvsize;
    nrecvproc = irregular->create_data_variable(nsendproc, proclist, sizelist,
        recvsize, 1); //must sort
    nrecvcell = recvsize / sizeof(CommMacro);
    //char* rbuf;
    memory->create(rbuf, recvsize, "gridCommMacro:rbuf");
    memset(rbuf, 0, recvsize);

    irregular->exchange_variable(sbuf, sizelist, rbuf);
    //cellint* recvidlist;
    //memory->create(recvidlist, nrecvcell * sizeof(cellint), "gridCommMacro:recvidlist");
    for (int i = 0; i < ncellsendall; ++i) {
        cellint id = 0;
        memcpy(&id, rbuf + i * sizeof(CommMacro), sizeof(cellint));
        if (grid->hash->find(id) != grid->hash->end()) {
            recvicelllist[i] = (*grid->hash)[id];
        }
        else {
            error->one(FLERR, "GridCommMacro : no such owned or ghost cell");
        }
           
    }
    //delete irregular;

    // unpack received grid cells as ghost cells

    //int offset = 0;
    //for (i = 0; i < nrecvproc; i++)
    //    offset += grid->unpack_one(&rbuf[offset], 0, 0, 0);

    // more clean up

    //memory->destroy(proclist);
    //memory->destroy(sizelist);
    //memory->destroy(sbuf);
    //memory->destroy(rbuf);
    


    // create send buf and auxiliary irregular comm vectors
    //char* sbuf;
    //memory->create(sbuf, sendsize, "grid:sbuf");
    //memset(sbuf, 0, sendsize);

    //int* proclist, * sizelist;
    //memory->create(proclist, nsend, "grid:proclist");
    //memory->create(sizelist, nsend, "grid:sizelist");

    //// on 2nd pass over local cells, fill the send buf
    //// use lastproc to insure a cell only overlaps once per other proc
    //// if oflag = 2 = my cell just touches box,
    //// so flag grid cell as EMPTY ghost by setting nsurf = -1

    //nsend = 0;
    //sendsize = 0;
    //for (int icell = 0; icell < nlocal; icell++) {
    //    if (cells[icell].nsplit <= 0) continue;
    //    lo = cells[icell].lo;
    //    hi = cells[icell].hi;
    //    lastproc = -1;
    //    for (i = 0; i < nlist; i++) {
    //        j = list[i];
    //        oflag = grid->box_overlap(lo, hi, boxall[j].lo, boxall[j].hi);
    //        if (oflag != 1) continue;
    //        if (boxall[j].proc == lastproc) continue;
    //        lastproc = boxall[j].proc;
    //        sendcelllist[sf[lastproc]++] = icell;
    //        //if (oflag == 2) {
    //        //    nsurf_hold = cells[icell].nsurf;
    //        //    cells[icell].nsurf = -1;
    //        //}
    //        //sizelist[nsend] = grid->pack_one(icell, &sbuf[sendsize], 0, 0, 0, 1);
    //        //if (oflag == 2) cells[icell].nsurf = nsurf_hold;
    //        //proclist[nsend] = lastproc;
    //        //sendsize += sizelist[nsend];
    //        nsend++;
    //    }
    //}
    //for (int i = 0; i < nprocs - 1; ++i) {
    //    if (sf[i] != sendfirst[i + 1])
    //        error->one(FLERR, "sendcelllist set error!");
    //}
    //delete[] sf; 
    //// clean up

    //memory->destroy(list);
    //delete[] boxall;

    //// perform irregular communication of list of ghost cells

    //Irregular* irregular = new Irregular(sparta);
    //int recvsize;
    //int nrecv = irregular->create_data_variable(nsend, proclist, sizelist,
    //    recvsize, 1); //must sort

    //char* rbuf;
    //memory->create(rbuf, recvsize, "grid:rbuf");
    //memset(rbuf, 0, recvsize);

    //irregular->exchange_variable(sbuf, sizelist, rbuf);
    //delete irregular;

    //// unpack received grid cells as ghost cells

    //int offset = 0;
    //for (i = 0; i < nrecv; i++)
    //    offset += grid->unpack_one(&rbuf[offset], 0, 0, 0);

    //// more clean up

    //memory->destroy(proclist);
    //memory->destroy(sizelist);
    //memory->destroy(sbuf);
    //memory->destroy(rbuf);

}

//int Grid::unpack_one_comm_list(char* buf,
//    int ownflag, int partflag, int surfflag, int sortflag)
//{
//    char* ptr = buf;
//
//    // unpack child cell as owned or ghost
//
//    int icell;
//    if (ownflag) icell = nlocal;
//    else icell = nlocal + nghost;
//    grow_cells(1, ownflag);
//    if (ownflag) nlocal++;
//    else nghost++;
//
//    memcpy(&cells[icell], ptr, sizeof(ChildCell));
//    ptr += sizeof(ChildCell);
//    ptr = ROUNDUP(ptr);
//
//    if (ownflag) {
//        cells[icell].proc = me;
//        cells[icell].ilocal = icell;
//    }
//
//    // no surfs or any other info
//    // ditto for EMPTY ghost with nsurf < 0
//    // reset other fields for ghost cell (csurfs, nsplit, isplit)
//
//    if (!surfflag || cells[icell].nsurf < 0) {
//        cells[icell].csurfs = NULL;
//        cells[icell].nsplit = 1;
//        cells[icell].isplit = -1;
//        return ptr - buf;
//    }
//
//    // if nsurfs, unpack different info for explicit vs implicit
//    // explicit all: list of csurf indices
//    // implicit: add entire list of lines or triangles
//    // explicit distributed:
//    //    list of lines or triangles
//    //    check hash each to see if can skip b/c already have it
//    //    add new surfs to hash
//
//    int nsurf = cells[icell].nsurf;
//    if (nsurf) {
//        cells[icell].csurfs = csurfs->vget();
//
//        // explicit all surfs
//
//        if (!surf->implicit && !surf->distributed) {
//            memcpy(cells[icell].csurfs, ptr, nsurf * sizeof(surfint));
//            ptr += nsurf * sizeof(surfint);
//            ptr = ROUNDUP(ptr);
//
//            // implicit surfs
//
//        }
//        else if (surf->implicit) {
//            if (domain->dimension == 2) {
//                int sizesurf = sizeof(Surf::Line);
//                surfint* csurfs = cells[icell].csurfs;
//                for (int m = 0; m < nsurf; m++) {
//                    Surf::Line* line = (Surf::Line*)ptr;
//                    surf->add_line_copy(ownflag, line);
//                    if (ownflag) csurfs[m] = surf->nlocal - 1;
//                    else csurfs[m] = surf->nlocal + surf->nghost - 1;
//                    ptr += sizesurf;
//                    ptr = ROUNDUP(ptr);
//                }
//            }
//            else {
//                int sizesurf = sizeof(Surf::Tri);
//                surfint* csurfs = cells[icell].csurfs;
//                for (int m = 0; m < nsurf; m++) {
//                    Surf::Tri* tri = (Surf::Tri*)ptr;
//                    surf->add_tri_copy(ownflag, tri);
//                    if (ownflag) csurfs[m] = surf->nlocal - 1;
//                    else csurfs[m] = surf->nlocal + surf->nghost - 1;
//                    ptr += sizesurf;
//                    ptr = ROUNDUP(ptr);
//                }
//            }
//
//            // explicit distributed surfs
//
//        }
//        else {
//            Surf::MySurfHash* shash = surf->hash;
//
//            if (domain->dimension == 2) {
//                int sizesurf = sizeof(Surf::Line);
//                surfint* csurfs = cells[icell].csurfs;
//                for (int m = 0; m < nsurf; m++) {
//                    Surf::Line* line = (Surf::Line*)ptr;
//                    if (shash->find(line->id) == shash->end()) {
//                        surf->add_line_copy(ownflag, line);
//                        if (ownflag) csurfs[m] = surf->nlocal - 1;
//                        else csurfs[m] = surf->nlocal + surf->nghost - 1;
//                        (*shash)[line->id] = csurfs[m];
//                    }
//                    else csurfs[m] = (*shash)[line->id];
//                    ptr += sizesurf;
//                    ptr = ROUNDUP(ptr);
//                }
//            }
//            else {
//                int sizesurf = sizeof(Surf::Tri);
//                surfint* csurfs = cells[icell].csurfs;
//                for (int m = 0; m < nsurf; m++) {
//                    Surf::Tri* tri = (Surf::Tri*)ptr;
//                    if (shash->find(tri->id) == shash->end()) {
//                        surf->add_tri_copy(ownflag, tri);
//                        if (ownflag) csurfs[m] = surf->nlocal - 1;
//                        else csurfs[m] = surf->nlocal + surf->nghost - 1;
//                        (*shash)[tri->id] = csurfs[m];
//                    }
//                    else csurfs[m] = (*shash)[tri->id];
//                    ptr += sizesurf;
//                    ptr = ROUNDUP(ptr);
//                }
//            }
//        }
//
//        csurfs->vgot(nsurf);
//    }
//
//    if (ownflag) {
//        memcpy(&cinfo[icell], ptr, sizeof(ChildInfo));
//        ptr += sizeof(ChildInfo);
//        ptr = ROUNDUP(ptr);
//    }
//
//    // if split cell, unpack sinfo and sinfo.csplits and sinfo.csubs
//    // create Nsplit sub cells
//    // use sinfo.csubs to set cells.ilocal for new sub cells
//    // create new csub for new sub cell indices
//    // if ownflag, also unpack volumes from sub cells themselves
//
//    if (cells[icell].nsplit > 1) {
//        int isplit;
//        if (ownflag) isplit = nsplitlocal;
//        else isplit = nsplitlocal + nsplitghost;
//        cells[icell].isplit = isplit;
//        add_split_cell(ownflag);
//        memcpy(&sinfo[isplit], ptr, sizeof(SplitInfo));
//        ptr += sizeof(SplitInfo);
//        ptr = ROUNDUP(ptr);
//
//        sinfo[isplit].icell = icell;
//        int nsurf = cells[icell].nsurf;
//        sinfo[isplit].csplits = csplits->vget();
//        memcpy(sinfo[isplit].csplits, ptr, nsurf * sizeof(int));
//        csplits->vgot(nsurf);
//        ptr += nsurf * sizeof(int);
//        ptr = ROUNDUP(ptr);
//
//        int nsplit = cells[icell].nsplit;
//        sinfo[isplit].csubs = csubs->vget();
//        memcpy(sinfo[isplit].csubs, ptr, nsplit * sizeof(int));
//        csubs->vgot(nsplit);
//        ptr += nsplit * sizeof(int);
//        ptr = ROUNDUP(ptr);
//
//        double* dptr;
//        if (ownflag) {
//            dptr = (double*)ptr;
//            ptr += nsplit * sizeof(double);
//        }
//
//        int isub;
//        for (int i = 0; i < nsplit; i++) {
//            if (ownflag) isub = nlocal;
//            else isub = nlocal + nghost;
//            add_sub_cell(icell, ownflag);
//            cells[isub].ilocal = sinfo[isplit].csubs[i];
//            cells[isub].nsplit = -i;
//            if (ownflag) cinfo[isub].volume = dptr[i];
//            sinfo[isplit].csubs[i] = isub;
//        }
//
//    }
//    else {
//        if (ownflag) nunsplitlocal++;
//        else nunsplitghost++;
//    }
//
//    // unpack collision and fix info for new grid cell
//
//    if (ownflag) {
//        if (collide) {
//            ptr += collide->unpack_grid_one(icell, ptr);
//            ptr = ROUNDUP(ptr);
//        }
//        if (modify->n_pergrid) {
//            ptr += modify->unpack_grid_one(icell, ptr);
//            ptr = ROUNDUP(ptr);
//        }
//    }
//
//    // unpack particles, for unsplit cell or split cell
//
//    if (!partflag) return ptr - buf;
//
//    ptr += unpack_particles(ptr, icell, sortflag);
//
//    // unpack particles of sub cells
//
//    if (cells[icell].nsplit > 1) {
//        int isplit = cells[icell].isplit;
//        int nsplit = cells[icell].nsplit;
//        for (int i = 0; i < nsplit; i++) {
//            int m = sinfo[isplit].csubs[i];
//            ptr += unpack_particles(ptr, m, sortflag);
//        }
//    }
//
//    return ptr - buf;
//}