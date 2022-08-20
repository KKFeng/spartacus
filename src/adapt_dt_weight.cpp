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

using namespace SPARTA_NS;


enum{SURF,VALUE,COMPUTE,FIX,SAME};
enum { UNKNOWN, OUTSIDE, INSIDE, OVERLAP };   // several files

#define INVOKED_PER_GRID 16
#define DELTA_LIST 1024
#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

AdaptDtWeight::AdaptDtWeight(SPARTA *sparta) : Pointers(sparta)
{
  me = comm->me;
  nprocs = comm->nprocs;
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

  if (style == SURF) set_weight_surf();
  else if (style == VALUE) set_weight_value();
  else if (style == SAME) set_weight_same();
  else error->all(FLERR, "wrong adapt_dt_weight_style");
  

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
   process command args for both adapt_grid and fix adapt
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
                if (!(lines[m].mask & sgroupbit)) continue;
            }
            else {
                if (!(tris[m].mask & sgroupbit)) continue;
            }
        }
        if (j == nsurf) continue;

        cells[icell].dt_weight = surf_ndt;
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
        //else {
        //    // split cells
        //    nsplit = cells[icell].nsplit;
        //    csubs = sinfo[cells[icell].isplit].csubs;
        //    value = -BIG;
        //    for (int j = 0; j < nsplit; j++) {
        //        jcell = csubs[j];
        //        if (valuewhich == COMPUTE)
        //            value = MAX(value, value_compute(jcell));
        //        else if (valuewhich == FIX)
        //            value = MAX(value, value_fix(jcell));
        //    }

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
        cells[icell].dt_weight = same_dt;
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