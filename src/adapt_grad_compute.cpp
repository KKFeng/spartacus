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

#include "string.h"
#include "adapt_grad_compute.h"
#include "update.h"
#include "domain.h"
#include "region.h"
#include "input.h"
#include "comm.h"
#include "grid.h"
#include "surf.h"
#include "modify.h"
#include "compute.h"
#include "fix.h"
#include "fix_ave_grid.h"
#include "error.h"
#include "memory.h"
#include "comm.h"
#include "grid_comm_macro.h"
#include "irregular.h"
#include "random_mars.h"
#include "random_park.h"

using namespace SPARTA_NS;

enum { XLO, XHI, YLO, YHI, ZLO, ZHI, INTERIOR };         // same as Domain
enum { NCHILD, NPARENT, NUNKNOWN, NPBCHILD, NPBPARENT, NPBUNKNOWN, NBOUND };  // Grid
enum{GRAD_MAX, GRAD_MIN, GRAD_NONE};
enum{CLEAR,VALUE,COMPUTE,FIX};
enum { UNKNOWN, OUTSIDE, INSIDE, OVERLAP };   // several files
//
//#define DELTA_RL 64 // how to grow region list
#define BIG 1.0e20
#define NVALUEMAX 5

/* ---------------------------------------------------------------------- */

AdaptGradCompute::AdaptGradCompute(SPARTA *sparta) : Pointers(sparta)
{
  me = comm->me;
  nprocs = comm->nprocs;
  mod = GRAD_NONE;
  int ncell = grid->nghost + grid->nlocal;
  q = new double[ncell];
  exist_q = new int[ncell] {};
}

/* ---------------------------------------------------------------------- */

AdaptGradCompute::~AdaptGradCompute()
{
    delete[] q;
    delete[] exist_q;
}

/* ---------------------------------------------------------------------- */

void AdaptGradCompute::command(int narg, char **arg)
{
  if (!grid->exist)
    error->all(FLERR,"Cannot compute grad when grid is not defined");

  if (narg < 1) error->all(FLERR,"Illegal adapt_grad_compute command");

  // process command-line args

  process_args(narg,arg);

  if (style == CLEAR) {
      grid->gradhashfilled = 0;
      grid->grad_l->clear();
      grid->grad_dt->clear();
      if (me == 0) {
          if (screen) fprintf(screen, "Gradient cleared succesfully !\n");
          if (logfile) fprintf(logfile, "Gradient cleared succesfully !\n");
      }
      return;
  }

  check_args();

  // perform adaptation

  if (me == 0) {
    if (screen) fprintf(screen,"Computing gradient for adapt ...\n");
    if (logfile) fprintf(logfile,"Computing gradient for adapt ...\n");
  }

  MPI_Barrier(world);
  double time1 = MPI_Wtime();

  if (style == VALUE) compute_grad_value();

  else error->all(FLERR, "wrong adapt_dt_weight style");

  //print_grad();

  MPI_Barrier(world);
  double time2 = MPI_Wtime();
  // stats

  double time_total = time2-time1;

  if (me == 0) {
    if (screen) {;
      fprintf(screen,"  CPU time = %g secs\n",time_total);
    }
    if (logfile) {
      fprintf(logfile,"  CPU time = %g secs\n",time_total);
    }
  }
}

/* ----------------------------------------------------------------------
   process command args for adapt_dt_weight
------------------------------------------------------------------------- */

void AdaptGradCompute::process_args(int narg, char **arg)
{
  if (narg == 1 && strcmp(arg[0], "clear") == 0) {
      style = CLEAR;
      return;
  }
  if (narg < 7) error->all(FLERR,"Illegal adapt_grad_compute command");

  int igroup = grid->find_group(arg[0]);
  if (igroup < 0) error->all(FLERR,"adapt_grad_compute group ID does not exist");
  groupbit = grid->bitmask[igroup];

  int iarg = 1;

  if (strcmp(arg[iarg], "value") == 0) 
  {
      if (iarg + 6 > narg) error->all(FLERR, "Illegal adapt_grad_compute command");
      style = VALUE;
      for (int i = 0; i < NVALUEMAX; ++i) {
          int tmp_iarg = iarg + 1 + i;
          if (strncmp(arg[tmp_iarg], "c_", 2) == 0) valuewhich_arr[i] = COMPUTE;
          else if (strncmp(arg[tmp_iarg], "f_", 2) == 0) valuewhich_arr[i] = FIX;
          else error->all(FLERR, "Illegal adapt_grad_compute command");

          int n = strlen(arg[tmp_iarg]);
          char* suffix = new char[n];
          strcpy(suffix, &arg[tmp_iarg][2]);

          char* ptr = strchr(suffix, '[');
          if (ptr) {
              if (suffix[strlen(suffix) - 1] != ']')
                  error->all(FLERR, "Illegal adapt_grad_compute command");
              valindex_arr[i] = atoi(ptr + 1);
              *ptr = '\0';
          }
          else valindex_arr[i] = 0;
          n = strlen(suffix) + 1;
          valueID_arr[i] = new char[n];
          strcpy(valueID_arr[i], suffix);
          delete[] suffix;
      }
      iarg += 6;

  }else error->all(FLERR,"Illegal adapt_grad_compute command");

  while (iarg < narg) {
      if (strcmp(arg[iarg], "mode") == 0) {
          if (iarg + 2 > narg) error->all(FLERR, "Illegal adapt_grad_compute command");
          if (strcmp(arg[iarg + 1], "max") == 0) mod = GRAD_MAX;
          else if (strcmp(arg[iarg + 1], "min") == 0) mod = GRAD_MIN;
          else if (strcmp(arg[iarg + 1], "none") == 0) mod = GRAD_NONE;
          else error->all(FLERR, "Illegal adapt_grad_compute command");
          iarg += 2;
      } else error->all(FLERR, "Illegal adapt_grad_compute command");
  }

}

/* ----------------------------------------------------------------------
   error check on value compute/fix for both adapt_dt_weight
------------------------------------------------------------------------- */

void AdaptGradCompute::check_args()
{
    // if fix ave/grid used require that:
    //   (1) fix adapt Nevery is multiple of fix ave Nfreq
    //   (2) fix ave/grid is defined before fix adapt (checked in fix adapt)

    if (style == VALUE) {
        for (int i = 0; i < NVALUEMAX; ++i) {
            if (valuewhich_arr[i] == COMPUTE) {
                icompute_arr[i] = modify->find_compute(valueID_arr[i]);
                if (icompute_arr[i] < 0)
                    error->all(FLERR, "Compute ID for adapt does not exist");
                if (modify->compute[icompute_arr[i]]->per_grid_flag == 0)
                    error->all(FLERR,
                        "Adapt compute does not calculate per-grid values");
                if (valindex_arr[i] == 0 && modify->compute[icompute_arr[i]]->size_per_grid_cols != 0)
                    error->all(FLERR, "Adapt compute does not calculate a per-grid vector");
                if (valindex_arr[i] && modify->compute[icompute_arr[i]]->size_per_grid_cols == 0)
                    error->all(FLERR, "Adapt compute does not calculate a per-grid array");
                if (valindex_arr[i] && valindex_arr[i] > modify->compute[icompute_arr[i]]->size_per_grid_cols)
                    error->all(FLERR, "Adapt compute array is accessed out-of-range");
            }
            else if (valuewhich_arr[i] == FIX) {
                ifix_arr[i] = modify->find_fix(valueID_arr[i]);
                if (ifix_arr[i] < 0)
                    error->all(FLERR, "Fix ID for adapt does not exist");
                if (modify->fix[ifix_arr[i]]->per_grid_flag == 0)
                    error->all(FLERR, "Adapt fix does not calculate per-grid values");
                if (valindex_arr[i] == 0 && modify->fix[ifix_arr[i]]->size_per_grid_cols != 0)
                    error->all(FLERR, "Adapt fix does not calculate a per-grid vector");
                if (valindex_arr[i] && modify->fix[ifix_arr[i]]->size_per_grid_cols == 0)
                    error->all(FLERR, "Adapt fix does not calculate a per-grid array");
                if (valindex_arr[i] && valindex_arr[i] > modify->fix[ifix_arr[i]]->size_per_grid_cols)
                    error->all(FLERR, "Adapt fix array is accessed out-of-range");
            }
        }
    }
}

void AdaptGradCompute::compute_grad_value() {
    Grid::MyGradHash* grad_l = grid->grad_l;
    Grid::MyGradHash* grad_dt = grid->grad_dt;
    struct AllGrad
    {
        cellint id;
        double dt, l;
    };
    if (mod == GRAD_NONE) {
        grad_l->clear();
        grad_dt->clear();
        grid->gradhashfilled = 0;
    }
    double value_arr[2]{};
    for (int i = 0; i < NVALUEMAX; ++i) {
        if (valuewhich_arr[i] == COMPUTE) {
            compute_arr[i] = modify->compute[icompute_arr[i]];
            compute_arr[i]->compute_per_grid();
            if (compute_arr[i]->post_process_grid_flag)
                compute_arr[i]->post_process_grid(valindex_arr[i], 1, NULL, NULL, NULL, 1);
        }
        else if (valuewhich_arr[i] == FIX) fix_arr[i] = modify->fix[ifix_arr[i]];
    }
    run_comm();
    Grid::ChildCell* cells = grid->cells;
    Grid::ChildInfo* cinfo = grid->cinfo;
    Grid::SplitInfo* sinfo = grid->sinfo;
    int nglocal = grid->nlocal;

    AllGrad* gradlist = new AllGrad[nglocal];
    int ngrad = 0, ngradmax = nglocal;

    for (int icell = 0; icell < nglocal; icell++) {
        if (!(cinfo[icell].mask & groupbit)) continue;
        if (cinfo[icell].type == INSIDE) continue;

        if (cells[icell].nsplit <= 1) {
            // unsplit cells or sub cells
            for (int i = 0; i < 2; ++i) {
                if (valuewhich_arr[i] == COMPUTE) value_arr[i] = value_compute(icell, i);
                else if (valuewhich_arr[i] == FIX) value_arr[i] = value_fix(icell, i);
            }
        }
        else continue;
        if (!exist_q[icell] || !q[icell]) {
            error->warning(FLERR, "!exist_q[icell] || !q[icell]");
            continue;
        }
        double grad = cal_grad(icell);
        double v_rms = sqrt(3 * update->boltz * value_arr[0] / value_arr[1]);
        double value_l = q[icell] / grad;
        double value_dt = q[icell] / grad / v_rms;
        if (isnan(value_l)|| isnan(value_dt)) {
            error->warning(FLERR, "value_l or value_dt is not a number, reset to BIG");
            value_l = value_dt = BIG;
        }
        cellint id = cells[icell].id;
        gradlist[ngrad].id = id;
        if ((mod == GRAD_NONE) || (grid->gradhashfilled == 0)
            || (grad_l->find(id) == grad_l->end())) {
            gradlist[ngrad].l = value_l;
            gradlist[ngrad].dt = value_dt;
        }
        else {
            if (mod == GRAD_MIN) {
                gradlist[ngrad].l = MIN(value_l, (*grad_l)[id]);
                gradlist[ngrad].dt = MIN(value_dt, (*grad_dt)[id]);
            }
            else {
                gradlist[ngrad].l = MAX(value_l, (*grad_l)[id]);
                gradlist[ngrad].dt = MAX(value_dt, (*grad_dt)[id]);
            }
        }
        ++ngrad;
    }
    int me = comm->me;
    int nprocs = comm->nprocs;
    int ngradall = 0;
    MPI_Allreduce(&ngrad, &ngradall, 1, MPI_INT, MPI_SUM, world);

    int* recvcounts, * displs;
    memory->create(recvcounts, nprocs, "grad:recvcounts");
    memory->create(displs, nprocs, "grad:displs");

    int nsend = ngrad * sizeof(AllGrad);
    MPI_Allgather(&nsend, 1, MPI_INT, recvcounts, 1, MPI_INT, world);
    displs[0] = 0;
    for (int i = 1; i < nprocs; i++) displs[i] = displs[i - 1] + recvcounts[i - 1];
    AllGrad* mygradlist = new AllGrad[ngrad];
    memcpy(mygradlist, gradlist, ngrad * sizeof(AllGrad));

    ngrad = ngradall;
    if (ngrad >= ngradmax) {
        ngradmax = ngrad;
        delete[] gradlist;
        gradlist = new AllGrad[ngradmax];
    }
    MPI_Allgatherv(mygradlist, nsend, MPI_CHAR, gradlist, recvcounts, displs, MPI_CHAR, world);

    delete[] mygradlist;
    memory->destroy(recvcounts);
    memory->destroy(displs);

    for (int i = 0; i < ngrad; ++i) {
        (*grad_l)[gradlist[i].id] = gradlist[i].l;
        (*grad_dt)[gradlist[i].id] = gradlist[i].dt;
    }
    delete[] gradlist;

    grid->gradhashfilled = 1;
}



/* ----------------------------------------------------------------------
   extract a value for icell,valindex from a compute
------------------------------------------------------------------------- */

double AdaptGradCompute::value_compute(int icell, int ico)
{
    double value;
    if (valindex_arr[ico] == 0 || compute_arr[ico]->post_process_grid_flag)
        value = compute_arr[ico]->vector_grid[icell];
    else value = compute_arr[ico]->array_grid[icell][valindex_arr[ico] - 1];
    return value;
}

/* ----------------------------------------------------------------------
   extract a value for icell,valindex from a fix
------------------------------------------------------------------------- */

double AdaptGradCompute::value_fix(int icell, int ifi)
{
    double value;
    if (valindex_arr[ifi] == 0) value = fix_arr[ifi]->vector_grid[icell];
    else value = fix_arr[ifi]->array_grid[icell][valindex_arr[ifi] - 1];
    return value;
}


void AdaptGradCompute::run_comm()
{
    GridCommMacro gcm(sparta);
    gcm.acquire_macro_comm_list_near();

    for (int icell = 0; icell < grid->nlocal; ++icell) {
        if (grid->cells[icell].nsplit <= 1) {
            // unsplit cells or sub cells
            double value_arr[3];
            for (int j = 2; j < 5; ++j) {
                if (valuewhich_arr[j] == COMPUTE) value_arr[j - 2] = value_compute(icell, j);
                else if (valuewhich_arr[j] == FIX) value_arr[j - 2] = value_fix(icell, j);
            }
            q[icell] = sqrt(value_arr[0] * value_arr[0] + value_arr[1] * value_arr[1]
                + value_arr[2] * value_arr[2]);
            exist_q[icell] = 1;
        }
    }
    // pack q, preparing for comm
    for (int i = 0; i < gcm.ncellsendall; ++i) {
        int icell = gcm.sendcelllist[i];
        if (icell < 0 || icell >= grid->ncell)
            error->all(FLERR, "cell index is out of own cell range");
        memcpy(gcm.sbuf + i * sizeof(CommMacro), &q[icell], sizeof(double));
        if (!exist_q[icell]) error->warning(FLERR, "using macro quantity of split cell");
    }    
    // communicating
    gcm.irregular->exchange_variable(gcm.sbuf, gcm.sizelist, gcm.rbuf);
    // unpack
    for (int i = 0; i < gcm.nrecvcell; ++i) {
        int icell = gcm.recvicelllist[i];
        if (icell < grid->nlocal || icell >= grid->nlocal + grid->nghost) {
            char str[128];
            sprintf(str, "cell index is out of ghost range. icell,nlocal,nghost = %d,%d,%d",
                icell, grid->nlocal , grid->nghost);
            error->warning(FLERR, str);

            continue;
        }
        exist_q[icell] = 1;
        memcpy(&q[icell], gcm.rbuf + i * sizeof(CommMacro), sizeof(double));
    }
}

// NOTE: currently cannot consider sub cell of ghost cell, if cell is split, ignore.
//       may fix it at some time
double AdaptGradCompute::cal_grad(int icell) {
    double grad[6]{};
    int exist_grad[6]{};
    int dim = domain->dimension;
    Grid::ChildInfo* cinfo = grid->cinfo;
    Grid::ChildCell* cells = grid->cells;
    Grid::ParentCell* pcells = grid->pcells;
    cellint *neigh = cells[icell].neigh;
    RanPark random(update->ranmaster->uniform());
    int nmask = cells[icell].nmask;
    double* lo = cells[icell].lo;
    double* hi = cells[icell].hi;
    int ic;
    int nsample = 10;
    for (int face = XLO; face <= ZHI; ++face) {
        if (face >= ZLO && dim != 3) break;
        int nflag = grid->neigh_decode(nmask, face);
        if (nflag == NCHILD || nflag == NPBCHILD) {
            ic = neigh[face];
            if (ic < 0 || cells[ic].nsplit > 1 || q[ic] == 0.0) continue;
            double dx = 0;
            for (int j = 0; j < 3; ++j) {
                if (dim != 3 && j == 2) break;
                double dxj = (cells[ic].hi[j] + cells[ic].lo[j] - hi[j] - lo[j]) / 2;
                dx += dxj * dxj;
            }
            dx = sqrt(dx);
            grad[face] = abs((q[ic] - q[icell]) / dx);
            exist_grad[face] = 1;
        } else if (nflag == NPARENT || nflag == NPBPARENT){
            double dx;
            if (face <= XHI) dx = hi[0] - lo[0];
            else if (face <= YHI) dx = hi[1] - lo[1];
            else  dx = hi[2] - lo[2];
            double sample_q = 0.0;
            int ngrad = 0;
            Grid::ParentCell* pcell = &pcells[neigh[face]];
            for (int i = 0; i < nsample; ++i) {
                double x[3];
                x[0] = pcell->lo[0] + random.uniform() * (pcell->hi[0] - pcell->lo[0]);
                x[1] = pcell->lo[1] + random.uniform() * (pcell->hi[1] - pcell->lo[1]);
                if (dim == 3)
                    x[2] = pcell->lo[2] + random.uniform() * (pcell->hi[2] - pcell->lo[2]);
                else x[2] = 0;
                ic = grid->id_find_child(pcell->id, cells[icell].level, pcell->lo, pcell->hi, x);
                if (ic < 0 || cells[ic].nsplit > 1 || q[ic] == 0.0) continue;
                sample_q += q[ic]; ++ngrad;
            }
            if (ngrad) {
                grad[face] = abs((sample_q/ngrad - q[icell]) / dx);
                exist_grad[face] = 1;
            }
        }
    }
    int useful = 0;
    double result = 0.0;
    int i = 0;
    while (i < 2 * dim) {
        if (exist_grad[i] && exist_grad[i + 1]) {
            result += (grad[i] + grad[i + 1])* (grad[i] + grad[i + 1]) / 4; ++useful;
        }
        else if (exist_grad[i] || exist_grad[i + 1]) {
            result += (grad[i] + grad[i + 1]) * (grad[i] + grad[i + 1]); ++useful;
        }
        i += 2;
    }
    if (useful == dim) return sqrt(result);
    else if (useful == 0) {
        error->warning(FLERR, "1 cell calulate grad failed, set to 0");
        return 0;
    } else {
        error->warning(FLERR, "1 cell calulate grad with not enough value");
        return sqrt(result / useful * dim);
    }
}

void AdaptGradCompute::print_grad() {
    Grid::MyGradHash* grad_l = grid->grad_l;
    Grid::MyGradHash* grad_dt = grid->grad_dt;
    Grid::ChildCell* cells = grid->cells;
    Grid::ChildInfo* cinfo = grid->cinfo;
    Grid::SplitInfo* sinfo = grid->sinfo;
    char str[128];
    sprintf(str, "grad.%d.dat", me);
    FILE* fp = fopen(str, "w");
    for (int icell = 0; icell < grid->nlocal; ++icell) {
        if (grid->cells[icell].nsplit > 1) continue;
        cellint id = cells[icell].id;
        double length = -1, time = -1;
        if ((grad_l->find(id) != grad_l->end())) {
            length = (*grad_l)[id];
            time = (*grad_dt)[id];
        }
        double x = (cells[icell].hi[0] + cells[icell].lo[0]) / 2;
        double y = (cells[icell].hi[1] + cells[icell].lo[1]) / 2;
        fprintf(fp, "%g %g %d %g %g %g %g\n", x, y, id, q[icell], cal_grad(icell), length, time);
    }
    fclose(fp);
}
