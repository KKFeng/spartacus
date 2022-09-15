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
#include "adapt_dt_weight.h"
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

using namespace SPARTA_NS;


enum{DT_MAX, DT_MIN, DT_NONE};
enum{SURF,NEAR_SURF,VALUE,COMPUTE,FIX,SAME};
enum { UNKNOWN, OUTSIDE, INSIDE, OVERLAP };   // several files

#define DELTA_RL 64 // how to grow region list
#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

AdaptDtWeight::AdaptDtWeight(SPARTA *sparta) : Pointers(sparta)
{
  me = comm->me;
  nprocs = comm->nprocs;
  mod = DT_NONE;
  nregion = maxregion = 0;
  regionlist = NULL;
}

/* ---------------------------------------------------------------------- */

AdaptDtWeight::~AdaptDtWeight()
{

}

/* ---------------------------------------------------------------------- */

void AdaptDtWeight::command(int narg, char **arg)
{
  if (!grid->exist)
    error->all(FLERR,"Cannot adapt dt_weight of grid when grid is not defined");

  grid->is_dt_weight = 1;

  if (narg < 1) error->all(FLERR,"Illegal adapt_grid command");

  // process command-line args

  process_args(narg,arg);
  check_args();

  // perform adaptation

  if (me == 0) {
    if (screen) fprintf(screen,"Adapting grid timestep ...\n");
    if (logfile) fprintf(logfile,"Adapting grid timestep ...\n");
  }

  MPI_Barrier(world);
  double time1 = MPI_Wtime();

  sparta->init();
  grid->remove_ghosts();

  if (style == SURF) set_weight_surf();
  else if (style == NEAR_SURF) set_weight_nearsurf();
  else if (style == VALUE) set_weight_value();
  else if (style == SAME) set_weight_same();
  else error->all(FLERR, "wrong adapt_dt_weight_style");
  
  grid->acquire_ghosts();
  grid->find_neighbors();

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

void AdaptDtWeight::process_args(int narg, char **arg)
{
  if (narg < 3) error->all(FLERR,"Illegal adapt command");

  int igroup = grid->find_group(arg[0]);
  if (igroup < 0) error->all(FLERR,"Adapt_dt_weight group ID does not exist");
  groupbit = grid->bitmask[igroup];

  int iarg = 1;

  if (strcmp(arg[iarg],"surf") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal adapt command");
      style = SURF;
      int igroup = surf->find_group(arg[iarg+1]);
      if (igroup < 0)
        error->all(FLERR,"Adapt command surface group does not exist");
      sgroupbit = surf->bitmask[igroup];
      surf_ndt = input->inumeric(FLERR,arg[iarg+2]);
      iarg += 3;

  } else if (strcmp(arg[iarg], "near_surf") == 0) {
      if (iarg + 4 > narg) error->all(FLERR, "Illegal adapt command");
      style = NEAR_SURF;
      int igroup = surf->find_group(arg[iarg + 1]);
      if (igroup < 0)
          error->all(FLERR, "Adapt command surface group does not exist");
      sgroupbit = surf->bitmask[igroup];
      surf_dist = input->numeric(FLERR, arg[iarg + 2]);
      surf_ndt = input->inumeric(FLERR, arg[iarg + 3]);
      iarg += 4;

  } else if (strcmp(arg[iarg],"value") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal adapt command");
      style = VALUE;
      if (strncmp(arg[iarg+1],"c_",2) == 0) valuewhich = COMPUTE;
      else if (strncmp(arg[iarg+1],"f_",2) == 0) valuewhich = FIX;
      else error->all(FLERR,"Illegal adapt command");

      int n = strlen(arg[iarg+1]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg+1][2]);

      char *ptr = strchr(suffix,'[');
      if (ptr) {
	if (suffix[strlen(suffix)-1] != ']')
	  error->all(FLERR,"Illegal adapt command");
	valindex = atoi(ptr+1);
	*ptr = '\0';
      } else valindex = 0;
      n = strlen(suffix) + 1;
      valueID = new char[n];
      strcpy(valueID,suffix);
      delete [] suffix;

      thresh = input->numeric(FLERR,arg[iarg+2]);
      max_dt = input->inumeric(FLERR,arg[iarg+3]);
      iarg += 4;

  } else if (strcmp(arg[iarg],"same") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal adapt command");
      style = SAME;
      same_dt = input->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;

  } else error->all(FLERR,"Illegal adapt command");

  while (iarg < narg) {
      if (strcmp(arg[iarg], "mode") == 0) {
          if (iarg + 2 > narg) error->all(FLERR, "Illegal adapt command");
          if (strcmp(arg[iarg + 1], "max") == 0) mod = DT_MAX;
          else if (strcmp(arg[iarg + 1], "min") == 0) mod = DT_MIN;
          else if (strcmp(arg[iarg + 1], "none") == 0) mod = DT_NONE;
          else error->all(FLERR, "Illegal adapt command");
          iarg += 2;
      }      
      else error->all(FLERR, "Illegal adapt command");
  }

}

/* ----------------------------------------------------------------------
   error check on value compute/fix for both adapt_dt_weight
------------------------------------------------------------------------- */

void AdaptDtWeight::check_args()
{
    // if fix ave/grid used require that:
    //   (1) fix adapt Nevery is multiple of fix ave Nfreq
    //   (2) fix ave/grid is defined before fix adapt (checked in fix adapt)

    if (style != VALUE) return;

    if (valuewhich == COMPUTE) {
        icompute = modify->find_compute(valueID);
        if (icompute < 0)
            error->all(FLERR, "Compute ID for adapt does not exist");
        if (modify->compute[icompute]->per_grid_flag == 0)
            error->all(FLERR,
                "Adapt compute does not calculate per-grid values");
        if (valindex == 0 && modify->compute[icompute]->size_per_grid_cols != 0)
            error->all(FLERR, "Adapt compute does not calculate a per-grid vector");
        if (valindex && modify->compute[icompute]->size_per_grid_cols == 0)
            error->all(FLERR, "Adapt compute does not calculate a per-grid array");
        if (valindex && valindex > modify->compute[icompute]->size_per_grid_cols)
            error->all(FLERR, "Adapt compute array is accessed out-of-range");

    }
    else if (valuewhich == FIX) {
        ifix = modify->find_fix(valueID);
        if (ifix < 0)
            error->all(FLERR, "Fix ID for adapt does not exist");
        if (modify->fix[ifix]->per_grid_flag == 0)
            error->all(FLERR, "Adapt fix does not calculate per-grid values");
        if (valindex == 0 && modify->fix[ifix]->size_per_grid_cols != 0)
            error->all(FLERR, "Adapt fix does not calculate a per-grid vector");
        if (valindex && modify->fix[ifix]->size_per_grid_cols == 0)
            error->all(FLERR, "Adapt fix does not calculate a per-grid array");
        if (valindex && valindex > modify->fix[ifix]->size_per_grid_cols)
            error->all(FLERR, "Adapt fix array is accessed out-of-range");
    }
}


void AdaptDtWeight::set_weight_surf() {
    Grid::ChildCell* cells = grid->cells;
    Grid::ChildInfo* cinfo = grid->cinfo;
    Surf::Line* lines = surf->lines;
    int dim = domain->dimension;
    Surf::Tri* tris = surf->tris;
    int nglocal = grid->nlocal;
    for (int icell = 0; icell < nglocal; icell++) {
        if (!(cinfo[icell].mask & groupbit)) continue;
        if (cinfo[icell].type != OVERLAP) continue;
        if (!cells[icell].nsurf) continue;
        int nsurf = cells[icell].nsurf;
        surfint* csurfs = cells[icell].csurfs;
        int j;
        for (j = 0; j < nsurf; j++) {
            int m = csurfs[j];
            if (dim == 2) {
                if ((lines[m].mask & sgroupbit)) break;
            }
            else {
                if ((tris[m].mask & sgroupbit)) break;
            }
        }
        if (j == nsurf) continue;
        if (mod==DT_MAX) cells[icell].dt_weight = MAX(surf_ndt, cells[icell].dt_weight);
        else if (mod==DT_MIN) cells[icell].dt_weight = MIN(surf_ndt, cells[icell].dt_weight);
        else cells[icell].dt_weight = surf_ndt;
    }
}

void AdaptDtWeight::set_weight_nearsurf() {
    Grid::ChildCell* cells = grid->cells;
    Grid::ChildInfo* cinfo = grid->cinfo;
    Surf::Line* lines = surf->lines;
    int dim = domain->dimension;
    Surf::Tri* tris = surf->tris;
    int nglocal = grid->nlocal;
    nregion = 0;
    for (int icell = 0; icell < nglocal; icell++) {
        if (cinfo[icell].type != OVERLAP) continue;
        if (!cells[icell].nsurf) continue;
        int nsurf = cells[icell].nsurf;
        surfint* csurfs = cells[icell].csurfs;
        int j;
        for (j = 0; j < nsurf; j++) {
            int m = csurfs[j];
            if (dim == 2) {
                if ((lines[m].mask & sgroupbit)) break;
            }
            else {
                if ((tris[m].mask & sgroupbit)) break;
            }
        }
        if (j == nsurf) continue;
        double x[3];
        for (int i = 0; i < 3; ++i) {
            x[i] = (cells[icell].lo[i] + cells[icell].hi[i])/2;
        }
        add_region(x, surf_dist);
    }

    gather_allregion();

    if (!nregion) {
        error->warning(FLERR, "adapt_dt_weight command: no specific group surf exist");
        return;
    }

    MyRegion maxregion = regionlist[0];
    double lo[3]{ BIG, BIG, BIG }, hi[3]{ -BIG, -BIG, -BIG };
    for (int i = 1; i < nregion; ++i) {
        for (int j = 0; j < 3; ++j) {
            lo[j] = MIN(lo[j], regionlist[i].x[j]);
            hi[j] = MAX(hi[j], regionlist[i].x[j]);
        }
    }
    if (dim == 3) {
        for (int j = 0; j < 3; ++j) maxregion.x[j] = (hi[j] + lo[j]) / 2;
        maxregion.radius = sqrt((hi[0] - lo[0]) * (hi[0] - lo[0])
            + (hi[1] - lo[1]) * (hi[1] - lo[1]) 
            + (hi[2] - lo[2]) * (hi[2] - lo[2])) / 2 + surf_dist;
    } else {
        maxregion.x[0] = (hi[0] + lo[0]) / 2;
        maxregion.x[1] = (hi[1] + lo[1]) / 2;
        maxregion.x[2] = 0.5;
        maxregion.radius = sqrt((hi[0] - lo[0]) * (hi[0] - lo[0])
            + (hi[1] - lo[1]) * (hi[1] - lo[1])) / 2 + surf_dist;
    }

    for (int icell = 0; icell < nglocal; icell++) {
        if (!(cinfo[icell].mask & groupbit)) continue;
        if (cinfo[icell].type == INSIDE) continue;
        double x[3];
        int i;
        for (i = 0; i < 3; ++i) x[i] = (cells[icell].lo[i] + cells[icell].hi[i]) / 2;
        if (!in_region(maxregion, x)) continue;
        for (i = 0; i < nregion; ++i) {
            if (in_region(regionlist[i], x)) break;
        }
        if (i == nregion) continue;
        if (mod == DT_MAX) cells[icell].dt_weight = MAX(surf_ndt, cells[icell].dt_weight);
        else if (mod == DT_MIN) cells[icell].dt_weight = MIN(surf_ndt, cells[icell].dt_weight);
        else cells[icell].dt_weight = surf_ndt;
    }
}


void AdaptDtWeight::set_weight_value() {
    int icell, nsplit, jcell;
    double value;
    int* csubs;

    // invoke compute each time refinement is done
    // grid could have changed from previous refinement or coarsening

    if (valuewhich == COMPUTE) {
        compute = modify->compute[icompute];
        compute->compute_per_grid();
        if (compute->post_process_grid_flag)
            compute->post_process_grid(valindex, 1, NULL, NULL, NULL, 1);
    }
    else if (valuewhich == FIX) fix = modify->fix[ifix];

    Grid::ChildCell* cells = grid->cells;
    Grid::ChildInfo* cinfo = grid->cinfo;
    Grid::SplitInfo* sinfo = grid->sinfo;
    int nglocal = grid->nlocal;
    for (int icell = 0; icell < nglocal; icell++) {
        if (!(cinfo[icell].mask & groupbit)) continue;
        if (cinfo[icell].type == INSIDE) continue;

        if (cells[icell].nsplit <= 1) {
            // unsplit cells or sub cells
            if (valuewhich == COMPUTE) value = value_compute(icell);
            else if (valuewhich == FIX) value = value_fix(icell);
        }
        else continue;
        cells[icell].dt_weight = (int)MIN(max_dt,MAX(1,value/thresh));
    }
}

void AdaptDtWeight::set_weight_same() {
    Grid::ChildCell* cells = grid->cells;
    Grid::ChildInfo* cinfo = grid->cinfo;
    int nglocal = grid->nlocal;
    for (int icell = 0; icell < nglocal; icell++) {
        if (!(cinfo[icell].mask & groupbit)) continue;
        if (cinfo[icell].type == INSIDE) continue;
        if (mod == DT_MAX) cells[icell].dt_weight = MAX(same_dt, cells[icell].dt_weight);
        else if (mod == DT_MIN) cells[icell].dt_weight = MIN(same_dt, cells[icell].dt_weight);
        else cells[icell].dt_weight = same_dt;
    }
}


/* ----------------------------------------------------------------------
   extract a value for icell,valindex from a compute
------------------------------------------------------------------------- */

double AdaptDtWeight::value_compute(int icell)
{
    double value;

    if (valindex == 0 || compute->post_process_grid_flag)
        value = compute->vector_grid[icell];
    else value = compute->array_grid[icell][valindex - 1];

    return value;
}

/* ----------------------------------------------------------------------
   extract a value for icell,valindex from a fix
------------------------------------------------------------------------- */

double AdaptDtWeight::value_fix(int icell)
{
    double value;

    if (valindex == 0) value = fix->vector_grid[icell];
    else value = fix->array_grid[icell][valindex - 1];

    return value;
}

void AdaptDtWeight::add_region(double* x, double radius) {
    if (nregion == maxregion) {
        maxregion += DELTA_RL;
        memory->grow(regionlist, maxregion, "adapt_dt_weight:regionlist");
    }
    memcpy(regionlist[nregion].x, x, 3 * sizeof(double));
    regionlist[nregion].radius = radius;
    ++nregion;
}

bool AdaptDtWeight::in_region(MyRegion& region, double* x) {
    if (domain->dimension == 3) {
        if ((region.x[0] - x[0]) * (region.x[0] - x[0])
            + (region.x[1] - x[1]) * (region.x[1] - x[1])
            + (region.x[2] - x[2]) * (region.x[2] - x[2]) <= region.radius * region.radius)
            return true;
    } else {
        if ((region.x[0] - x[0]) * (region.x[0] - x[0])
            + (region.x[1] - x[1]) * (region.x[1] - x[1]) <= region.radius * region.radius)
            return true;
    }
    return false;
}

void AdaptDtWeight::gather_allregion() {
    int me = comm->me;
    int nprocs = comm->nprocs; 
    int nregionall;
    MPI_Allreduce(&nregion, &nregionall, 1, MPI_INT, MPI_SUM, world);

    int* recvcounts, * displs;
    memory->create(recvcounts, nprocs, "grid:recvcounts");
    memory->create(displs, nprocs, "grid:displs");

    int nsend = nregion * sizeof(MyRegion);
    MPI_Allgather(&nsend, 1, MPI_INT, recvcounts, 1, MPI_INT, world);
    displs[0] = 0;
    for (int i = 1; i < nprocs; i++) displs[i] = displs[i - 1] + recvcounts[i - 1];
    MyRegion* myregionlist = new MyRegion[nregion];
    memcpy(myregionlist, regionlist, nregion * sizeof(MyRegion));
    nregion = nregionall;
    if (nregion >= maxregion) {
        maxregion = nregion;
        memory->grow(regionlist, maxregion, "adapt_dt_weight:regionlist");
    }
    MPI_Allgatherv(myregionlist, nsend, MPI_CHAR, regionlist, recvcounts, displs, MPI_CHAR, world);

    delete[] myregionlist;
    memory->destroy(recvcounts);
    memory->destroy(displs);
}
