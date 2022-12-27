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

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "grid.h"
#include "update.h"
#include "particle.h"
#include "react.h"
#include "comm.h"
#include "random_park.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "collide_bgk.h"
#include "grid_comm_macro.h"
#include "surf_collide.h"
#include "output.h"
#include "mpi.h"

using namespace SPARTA_NS;
using namespace MathConst;

#define MAXLINE 1024
enum { USP, BGK, ESBGK, SBGK };
/* ---------------------------------------------------------------------- */

CollideBGK::CollideBGK(SPARTA* sparta, int narg, char** arg) :
    Collide(sparta, narg, arg)
{
    if (narg < 4) error->all(FLERR, "Illegal collide command");

    nmaxconserv = 0;
    conservMacro = NULL;

    // proc 0 reads file to extract params for current species
    // broadcasts params to all procs
    time_ave_coef = 0.99;
    nparams = particle->nspecies;
    resetWmax = 0.99;
    Pr = 0.666667;
    interpolate_flag = 1;
    if (nparams == 0)
        error->all(FLERR, "Cannot use collide command with no species defined");

    memory->create(params, nparams, "collide_bgk:params");
    if (strcmp(arg[2], "usp") == 0) {
        bgk_mod = USP;
    } else if (strcmp(arg[2], "bgk") == 0) {
        bgk_mod = BGK;
    }
    else if (strcmp(arg[2], "esbgk") == 0) {
        bgk_mod = ESBGK;
    }
    else if (strcmp(arg[2], "sbgk") == 0) {
        bgk_mod = SBGK;
    }
    else error->all(FLERR, "Illegal collide_bgk command: no such mod");
    if (narg > 4) {
        CollideBGKModify bgk_modify = CollideBGKModify(sparta);
        bgk_modify.command(narg - 4, arg + 4);
    }

    if (comm->me == 0) read_param_file(arg[3]);
    MPI_Bcast(params, nparams * sizeof(Params), MPI_BYTE, 0, world);

    count_try_relaxation = count_done_relaxation = count_fail_relaxation = 0;
    count_do_childcell = count_ignore_childcell = count_warning_ignore_childcell = 0;

    maxglocal = 0;
    resetWmax_flag = NULL;
    nplocalmax = 0;
    relax_flag = NULL;
}


/* ---------------------------------------------------------------------- */

CollideBGK::~CollideBGK()
{
    if (copymode) return;

    memory->destroy(params);
    //memory->destroy(prefactor);
}

/* ---------------------------------------------------------------------- */
void CollideBGK::reset_count() {
    count_try_relaxation = count_done_relaxation = count_fail_relaxation = 0;
    count_do_childcell = count_ignore_childcell = count_warning_ignore_childcell = 0;
}

/* ----------------------------------------------------------------------
* currently CollideBGK::init() will do nothing but call Collide::init();
   ---------------------------------------------------------------------- */

void CollideBGK::init()
{
    Collide::init();
}

/* ----------------------------------------------------------------------
* perform BGK-like collisions of all child cells I own, call perform_***bgk()
* to do per-particle job according to different bgk_mod
------------------------------------------------------------------------- */

void CollideBGK::collisions()
{
    // computing macro quantities for each model
    if (bgk_mod == USP) computeMacro<USP>();
    else if (bgk_mod == BGK) computeMacro<BGK>();
    else if (bgk_mod == SBGK) computeMacro<SBGK>();
    else if (bgk_mod == ESBGK) computeMacro<ESBGK>();

    if (nglocal > maxglocal) {
        maxglocal = ceil(nglocal * 1.2);
        memory->destroy(resetWmax_flag);
        memory->create(resetWmax_flag, maxglocal, "collideBGK:resetWmax_flag");
    }
    for (int icell = 0; icell < nglocal; icell++) resetWmax_flag[icell] = 1;
    reset_relaxflag();

    // loop over cells I own
    Grid::ChildCell* cells = grid->cells;
    Grid::ChildInfo* cinfo = grid->cinfo;
    Particle::OnePart* particles = particle->particles;
    int* next = particle->next;
    for (int icell = 0; icell < nglocal; icell++) {
        int np = cinfo[icell].count;
        if (!grid->cinfo[icell].macro.do_relaxation) continue;
        int ip = cinfo[icell].first;
        double volume = cinfo[icell].volume / cinfo[icell].weight * cells[icell].dt_weight;
        if (volume == 0.0) error->one(FLERR, "Collision cell volume is zero");

        // setup particle list for this cell

        if (np > npmax) {
            while (np > npmax) npmax += DELTAPART;
            memory->destroy(plist);
            memory->create(plist, npmax, "collide:plist");
        }

        Grid::ChildCell* cells = grid->cells;
        double bgk_attempt = attempt_collision(icell,0,cinfo[icell].macro.tao*2.0);
        int bgk_nattempt = static_cast<int> (bgk_attempt + (random->uniform()));

        int n = 0;
        while (ip >= 0) {
            plist[n++] = ip;
            ip = next[ip];
        }
        // Randomly swap particle lists, select the first bgk_nattempt part to relax
        if (bgk_nattempt < np / 2) {
            for (int i = 0; i < bgk_nattempt; i++) {
                int t = i + (np - i) * random->uniform();
                std::swap(plist[i], plist[t]);
            }
        }
        else {
            for (int i = np - 1; i > bgk_nattempt - 1; i--) {
                int t = i * random->uniform();
                std::swap(plist[i], plist[t]);
            }
        }
        for (int i = 0; i < bgk_nattempt; ++i) {
            relax_flag[plist[i]] = 1;
        }
    }

    // loop over all my part to improve cache hit ratio
    for (int i = 0; i < particle->nlocal; ++i) {
        if (!relax_flag[i]) continue;
        Particle::OnePart* ipart = &particles[i];
        int icell = ipart->icell;
        const CommMacro* interMacro = &grid->cells[icell].macro;
        if (interpolate_flag) {
            interMacro = grid->gridCommMacro->interpolation(ipart);
            if ((!interMacro) || (!(interMacro->Temp > 0))) {
                if (!interMacro)
                    error->warning(FLERR, "CollideBGK:interpolation failed!(!interMacro)");
                interMacro = &grid->cells[icell].macro;
            }
        }
        if (bgk_mod == USP) perform_uspbgk(ipart, icell, interMacro);
        else if (bgk_mod == BGK) perform_bgkbgk(ipart, icell, interMacro);
        else if (bgk_mod == SBGK) perform_sbgk(ipart, icell, interMacro);
        else if (bgk_mod == ESBGK) perform_esbgk(ipart, icell, interMacro);
    }

    for (int icell = 0; icell < nglocal; icell++) {
        if (resetWmax > 0.0 && resetWmax_flag[icell] &&
            (bgk_mod == USP || bgk_mod == SBGK))
            cinfo[icell].macro.Wmax *= resetWmax;
    }
    conservV();
    print_warning();
}

/* ----------------------------------------------------------------------
* Scale particles' new velocity to satisfy momentum & energy conservation
------------------------------------------------------------------------- */

void CollideBGK::conservV() {
    int nlocal = grid->nlocal;
    if (!(nmaxconserv >= 0)) error->one(FLERR,
        "CollideBGK::conservV(): !(nmaxconserv >= 0)");
    if (nlocal > nmaxconserv) {
        while (nlocal > nmaxconserv) nmaxconserv += DELTAPART;
        memory->destroy(conservMacro);
        memory->create(conservMacro, nmaxconserv, "collideBGK:postmacro");
    }
    for (int i = 0; i < nlocal; ++i) {
        NoCommMacro& nmacro = grid->cinfo[i].macro;
        nmacro.sum_vi[0] = nmacro.sum_vi[1] = nmacro.sum_vi[2] = 0.0;
        nmacro.sum_vij[0] = 0.0;
    }
    for (int ipart = 0; ipart < particle->nlocal; ++ipart) {
        Particle::OnePart& part = particle->particles[ipart];
        NoCommMacro& nmacro = grid->cinfo[part.icell].macro;
        for (int i = 0; i < 3; ++i) {
            nmacro.sum_vi[i] += part.v[i];
            nmacro.sum_vij[0] += part.v[i] * part.v[i];
        }
    }
    for (int icell = 0; icell < nlocal; ++icell) {
        conservMacro[icell].done_relaxation = grid->cinfo[icell].macro.do_relaxation;
        if (!conservMacro[icell].done_relaxation) continue;
        double theta = grid->cells[icell].macro.Temp / particle->species[0].mass * update->boltz;
        double np = grid->cinfo[icell].count;
        if (np <= 3) continue;
        NoCommMacro& nmacro = grid->cinfo[icell].macro;
        memcpy(conservMacro[icell].v_origin,
            grid->cells[icell].macro.v, sizeof(double) * 3);         
        memcpy(conservMacro[icell].v_post,
            nmacro.sum_vi, sizeof(double) * 3);
        for (int i = 0; i < 3; ++i) conservMacro[icell].v_post[i] /= np;
        double theta_post = (nmacro.sum_vij[0]
            - (nmacro.sum_vi[0] * nmacro.sum_vi[0] + nmacro.sum_vi[1] * nmacro.sum_vi[1]
                + nmacro.sum_vi[2] * nmacro.sum_vi[2]) / np) / np / 3;
        conservMacro[icell].coef = sqrt(theta / theta_post);
    }
    for (int ipart = 0; ipart < particle->nlocal; ++ipart) {
        Particle::OnePart& part = particle->particles[ipart];
        ConservMacro& cm = conservMacro[part.icell];
        if (!cm.done_relaxation) continue;
        for (int i = 0; i < 3; ++i) {
            part.v[i] = (part.v[i] - cm.v_post[i]) * cm.coef + cm.v_origin[i];
        }
    }
}

/* ----------------------------------------------------------------------
* perform per-part relaxation in differen mod: USP-BGK, original BGK, ES-BGK
* & SBGK, called by CollideBGK::collisions()
------------------------------------------------------------------------- */

void CollideBGK::perform_uspbgk(Particle::OnePart* ip, int icell, const CommMacro* interMacro)
{
    Grid::ChildInfo* cinfo = grid->cinfo;
    const double* sigma_ij = cinfo[icell].macro.sigma_ij;
    const double* q = cinfo[icell].macro.qi;
    double vn[3];
    int count_loop = 0;
    double theta = interMacro->Temp / particle->species[ip->ispecies].mass * update->boltz;
    while (true)
    {
        ++count_try_relaxation;
        ++count_loop;
        for (int i = 0; i < 3; i++) vn[i] = random->gaussian() * sqrt(theta);
        double C_2 = vn[0] * vn[0] + vn[1] * vn[1] + vn[2] * vn[2];
        double trace = C_2/3;
        double sigmacc =
            sigma_ij[0] * (vn[0] * vn[0] - trace)
            + sigma_ij[1] * (vn[1] * vn[1] - trace)
            + sigma_ij[2] * (vn[2] * vn[2] - trace)
            + sigma_ij[3] * vn[0] * vn[1] * 2
            + sigma_ij[4] * vn[0] * vn[2] * 2
            + sigma_ij[5] * vn[1] * vn[2] * 2;
        double qkck = (vn[0] * q[0] + vn[1] * q[1] + vn[2] * q[2]) *
            (C_2 / theta - 5);

        double W = 1.0 + cinfo[icell].macro.coef_A * sigmacc +
            cinfo[icell].macro.coef_B * qkck;
        if (W > cinfo[icell].macro.Wmax && W < 5) {
            cinfo[icell].macro.Wmax = W;
            resetWmax_flag[icell] = 0;
            break;
        }
        if (random->uniform() < W / cinfo[icell].macro.Wmax) break;

        if (count_loop > 100) {
            ++count_fail_relaxation;
            break;
        }
    }
    for (int i = 0; i < 3; i++) ip->v[i] = vn[i] + interMacro->v[i];
    ++count_done_relaxation;    
}

/* ---------------------------------------------------------------------- */

void CollideBGK::perform_bgkbgk(Particle::OnePart* ip, int , const CommMacro* interMacro)
{
    double theta = interMacro->Temp / particle->species[ip->ispecies].mass * update->boltz;
    for (int i = 0; i < 3; i++)
        ip->v[i] = random->gaussian() * sqrt(theta) + interMacro->v[i];
}

/* ---------------------------------------------------------------------- */

void CollideBGK::perform_esbgk(Particle::OnePart* ip, int icell, const CommMacro* interMacro)
{
    Grid::ChildInfo* cinfo = grid->cinfo;
    //(0, 1, 2, 3, 4, 5)
    //(00,11,22,01,02,12)
    const double* Sij = cinfo[icell].macro.sigma_ij;
    double vn[3];
    double theta = interMacro->Temp / particle->species[ip->ispecies].mass * update->boltz;
    for (int i = 0; i < 3; i++)
        vn[i] = random->gaussian() * sqrt(theta);
    ip->v[0] = vn[0]*Sij[0] + vn[1]*Sij[3] + vn[2]*Sij[4] +interMacro->v[0];
    ip->v[1] = vn[0]*Sij[3] + vn[1]*Sij[1] + vn[2]*Sij[5] +interMacro->v[1];
    ip->v[2] = vn[0]*Sij[4] + vn[1]*Sij[5] + vn[2]*Sij[2] +interMacro->v[2];
}

/* ---------------------------------------------------------------------- */

void CollideBGK::perform_sbgk(Particle::OnePart* ip, int icell, const CommMacro* interMacro)
{
    Grid::ChildInfo* cinfo = grid->cinfo;
    const double* q = cinfo[icell].macro.qi;
    double theta = interMacro->Temp / particle->species[ip->ispecies].mass * update->boltz;
    double vn[3];
    while (true)
    {
        for (int i = 0; i < 3; i++) vn[i] = random->gaussian() * sqrt(theta);
        double C_2 = vn[0] * vn[0] + vn[1] * vn[1] + vn[2] * vn[2];
        double qkck = (vn[0] * q[0] + vn[1] * q[1] + vn[2] * q[2]) *
            (C_2 / theta - 5);
        double W = 1.0 + cinfo[icell].macro.coef_B * qkck;
        if (W > cinfo[icell].macro.Wmax) {
            cinfo[icell].macro.Wmax = W;
            resetWmax_flag[icell] = 0;
            break;
        }
        if (random->uniform() < W / cinfo[icell].macro.Wmax) break;
    }
    for (int i = 0; i < 3; i++) ip->v[i] = vn[i] + interMacro->v[i];
}

/* ----------------------------------------------------------------------
   estimate a good value for vremax for a group pair in any grid cell
   called by Collide parent in init()

   NOTE: vremax is useless in BGK model, so always return 1.0
         
------------------------------------------------------------------------- */

double CollideBGK::vremax_init(int igroup, int jgroup)
{
    return 1.0;
}

/* ----------------------------------------------------------------------
* calculate number of part need relaxation this cell, 
* based on different bgk_mod
------------------------------------------------------------------------- */

double CollideBGK::attempt_collision(int icell, int, double tao)
{
    Grid::ChildInfo* cinfo = grid->cinfo;
    double fnum = update->fnum;
    double np = cinfo[icell].count;
    if (np < 4) {
        np += random->uniform() * 4;
        if (np < 4) return 0;
    }
    double bgk_nattempt;
    if (bgk_mod == ESBGK) bgk_nattempt = Pr * np * (1 - exp(-tao));
    else  bgk_nattempt = np * (1 - exp(-tao));
    return MIN(bgk_nattempt, (double)cinfo[icell].count);
}

/* ----------------------------------------------------------------------
  NOTE: perform_collision is replaced by perform_**bgk below, 
        should never be called
------------------------------------------------------------------------- */
int CollideBGK::perform_collision(Particle::OnePart*&,
    Particle::OnePart*&, Particle::OnePart*& )
{
    error->all(FLERR, "call perform_collision function of CollideBGK");
    return 0;
}

/* ----------------------------------------------------------------------
  compute macro quantities for all cells
  including velocity, temprature, shear stress and heat flux
  NOTE: Some computations may be omitted depending on BGK model
------------------------------------------------------------------------- */

template < int MOD > void CollideBGK::computeMacro() 
{
    for (int icell = 0; icell < nglocal; icell++)
    {
        NoCommMacro& nmacro = grid->cinfo[icell].macro;
        nmacro.sum_vi[0] = nmacro.sum_vi[1] = nmacro.sum_vi[2] = 0.0;
        nmacro.sum_vij[0] = nmacro.sum_vij[1] = nmacro.sum_vij[2] = 0.0;
        nmacro.sum_vij[3] = nmacro.sum_vij[4] = nmacro.sum_vij[5] = 0.0;
        nmacro.sum_C2vi[0] = nmacro.sum_C2vi[1] = nmacro.sum_C2vi[2] = 0.0;
    }

    // sum vi, vij viij for all child cells I own by iterating over all my part
    // Note: Currently only for single species !!!
    for (int ipart = 0; ipart < particle->nlocal; ++ipart) {
        Particle::OnePart& part = particle->particles[ipart];
        NoCommMacro& nmacro  = grid->cinfo[part.icell].macro;
        double* v = part.v;
        double C2 = 0.0;
        for (int i = 0; i < 3; ++i) {
            nmacro.sum_vi[i] += v[i];
            double vii = v[i] * v[i];
            nmacro.sum_vij[i] += vii;
            C2 += vii;
        }
        if (MOD == USP || MOD == ESBGK || MOD == SBGK) {
            nmacro.sum_vij[3] += v[0] * v[1];
            nmacro.sum_vij[4] += v[0] * v[2];
            nmacro.sum_vij[5] += v[1] * v[2];
        }
        if (MOD == USP || MOD == SBGK) {
            nmacro.sum_C2vi[0] += C2 * v[0];
            nmacro.sum_C2vi[1] += C2 * v[1];
            nmacro.sum_C2vi[2] += C2 * v[2];
        }  
    }
    
    for (int icell = 0; icell < nglocal; icell++)
    {
        NoCommMacro& mean_nmacro = grid->cinfo[icell].macro;
        CommMacro& cmacro = grid->cells[icell].macro;
        Grid::ChildCell& cell = grid->cells[icell];
        Grid::ChildInfo& cinfo = grid->cinfo[icell];
        Particle::OnePart* particles = particle->particles;
        int np = cinfo.count;
        if (np <= 3) {
            mean_nmacro.do_relaxation = 0;
            ++count_ignore_childcell;
            continue;
        }
        else
        {
            ++count_do_childcell;
            mean_nmacro.do_relaxation = 1;
        }
        // Currently assume all particles have same ispecies
        Particle::Species& 
            species = particle->species[particles[cinfo.first].ispecies];
        double mass = particle->species[particles[cinfo.first].ispecies].mass;
        Params& ps = params[particles[cinfo.first].ispecies];
        double pij[6]{}, qi[3]{};
        double* sum_vij = mean_nmacro.sum_vij;
        double* v = cmacro.v;
        for (int i = 0; i < 3; ++i) {
            v[i] = mean_nmacro.sum_vi[i] / np;
        }
        double V_2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
        double sum_C2 = sum_vij[0] + sum_vij[1] + sum_vij[2];
        cmacro.Temp = mass / update->boltz * (sum_C2 / np - V_2) / 3;
        if (!(cmacro.Temp > ps.T_ref * 0.01)) {
            // if particle is weighted, particles with same velocity maybe exist, thus
            // Temp ¡Ö 0 due to truncation error of floating point numbers
            ++count_warning_ignore_childcell;
            mean_nmacro.do_relaxation = 0;
            continue;
        }

        double nrho = cinfo.count * update->fnum * cinfo.weight / cell.dt_weight / cinfo.volume;
        mean_nmacro.tao = nrho * update->boltz * pow(ps.T_ref, ps.omega)
            * pow(cmacro.Temp, 1 - ps.omega) * update->dt / cell.dt_weight / ps.mu_ref / 2.0;
        double p = 0.0;
        if (MOD == USP|| MOD == SBGK) {
            double factor = mass * update->fnum * cinfo.weight / cell.dt_weight / cinfo.volume;
            for (int i = 0; i < 3; ++i) {
                pij[i] = factor * (sum_vij[i] - np * v[i] * v[i]);
            }
            pij[3] = factor * (sum_vij[3] - np * v[0] * v[1]);
            pij[4] = factor * (sum_vij[4] - np * v[0] * v[2]);
            pij[5] = factor * (sum_vij[5] - np * v[1] * v[2]);
            // time-average pij
            p = (pij[0] + pij[1] + pij[2]) / 3.0;
            for (int i = 0; i < 3; ++i) {
                mean_nmacro.sigma_ij[i] = mean_nmacro.sigma_ij[i] * time_ave_coef
                    + (pij[i] - p) * (1 - time_ave_coef) / (1 + mean_nmacro.tao);
            }
            for (int i = 3; i < 6; ++i) {
                mean_nmacro.sigma_ij[i] = mean_nmacro.sigma_ij[i] * time_ave_coef
                    + pij[i] * (1 - time_ave_coef) / (1 + mean_nmacro.tao);
            }
            double factor_q = factor / 
                (1 + Pr * mean_nmacro.tao);
            qi[0] = factor_q / 2 * (mean_nmacro.sum_C2vi[0]
                - v[0] * sum_C2 + 2 * np * V_2 * v[0]
                - 2 * (v[0] * sum_vij[0] + v[1] * sum_vij[3] + v[2] * sum_vij[4]));
            qi[1] = factor_q / 2 * (mean_nmacro.sum_C2vi[1]
                - v[1] * sum_C2 + 2 * np * V_2 * v[1]
                - 2 * (v[0] * sum_vij[3] + v[1] * sum_vij[1] + v[2] * sum_vij[5]));
            qi[2] = factor_q / 2 * (mean_nmacro.sum_C2vi[2]
                - v[2] * sum_C2 + 2 * np * V_2 * v[2]
                - 2 * (v[0] * sum_vij[4] + v[1] * sum_vij[5] + v[2] * sum_vij[2]));

            // time-average qi
            for (int i = 0; i < 3; ++i) {
                mean_nmacro.qi[i] = mean_nmacro.qi[i] * time_ave_coef
                + qi[i] * (1.0 - time_ave_coef) / (1 + Pr * mean_nmacro.tao);
            }
            // prefactor of weight in Acceptance-Rejection Method
            if (MOD == USP) {
                double p_theta = p * cmacro.Temp / mass * update->boltz;
                double tao_coth = mean_nmacro.tao * (1.0 + 2.0 / (exp(mean_nmacro.tao * 2.0) - 1.0));
                mean_nmacro.coef_A = (1.0 - tao_coth) / (2.0 * p_theta);
                mean_nmacro.coef_B = (1.0 - Pr * tao_coth) / (5.0 * p_theta);
            }
            else if (MOD == SBGK) {
                mean_nmacro.coef_B = (1.0 - Pr) / (5.0 * p * cmacro.Temp / mass * update->boltz);
            }
        }
        // NOTE: if MOD == ESBGK, sigma_ij is actually Sij in esbgk mod, no time-ave
        else if (MOD == ESBGK) {
            double* vi = mean_nmacro.sum_vi;
            double pf_Pr = (1.0 - Pr) / (Pr * 2.0);
            double pf_T = ((sum_vij[0] + sum_vij[1] + sum_vij[2]) -
                (vi[0] * vi[0] + vi[1] * vi[1] + vi[2] * vi[2]) / np) / 3.0;
            for (int i = 0; i < 3; ++i) {
                mean_nmacro.sigma_ij[i] = 1 + pf_Pr - pf_Pr / pf_T *
                    (sum_vij[i] - vi[i] * vi[i] / np);
            }
            mean_nmacro.sigma_ij[3] = - pf_Pr / pf_T *
                (sum_vij[3] - vi[0] * vi[1] / np);            
            mean_nmacro.sigma_ij[4] = - pf_Pr / pf_T *
                (sum_vij[4] - vi[0] * vi[2] / np);            
            mean_nmacro.sigma_ij[5] = - pf_Pr / pf_T *
                (sum_vij[5] - vi[1] * vi[2] / np);
        }
    }

    // run commMacro
    grid->gridCommMacro->runComm();

}

/* ----------------------------------------------------------------------
   read list of species defined in species file
   store info in filespecies and nfilespecies
   only invoked by proc 0
------------------------------------------------------------------------- */

void CollideBGK::read_param_file(char* fname)
{
    FILE* fp = fopen(fname, "r");
    if (fp == NULL) {
        char str[128];
        sprintf(str, "Cannot open BGK parameter file %s", fname);
        error->one(FLERR, str);
    }

    // set all species diameters to -1, so can detect if not read
    // set all cross-species parameters to -1 to catch no-reads, as
    // well as user-selected average

    for (int i = 0; i < nparams; i++) {
        params[i].mu_ref = -1.0;
    }

    // read file line by line
    // skip blank lines or comment lines starting with '#'
    // all other lines must have at least REQWORDS, which depends on VARIABLE flag

    int REQWORDS = 4;
    char** words = new char* [REQWORDS]; 
    char line[MAXLINE];
    int isp;

    while (fgets(line, MAXLINE, fp)) {
        int pre = strspn(line, " \t\n\r");
        if (pre == strlen(line) || line[pre] == '#') continue;

        int nwords = wordparse(REQWORDS, line, words);
        if (nwords < REQWORDS)
            error->one(FLERR, "Incorrect line format in BGK parameter file");

        isp = particle->find_species(words[0]);
        if (isp < 0) continue;

        else {
            if (nwords < REQWORDS) 
                error->one(FLERR, "Incorrect line format in BGK parameter file");
            params[isp].mu_ref  = atof(words[1]);
            params[isp].omega = atof(words[2]);
            params[isp].T_ref  = atof(words[3]);
        }
    }

    delete[] words;
    fclose(fp);

    // check that params were read for all species
    for (int i = 0; i < nparams; i++) {

        if (params[i].mu_ref < 0.0) {
            char str[128];
            sprintf(str, "Species %s did not appear in BGK parameter file",
                particle->species[i].id);
            error->one(FLERR, str);
        }
    }
}

/* ----------------------------------------------------------------------
   parse up to n=maxwords whitespace-delimited words in line
   store ptr to each word in words and count number of words
   same as CollideVSS::wordparse
------------------------------------------------------------------------- */

int CollideBGK::wordparse(int maxwords, char* line, char** words)
{
    int nwords = 1;
    char* word;

    words[0] = strtok(line, " \t\n");
    while ((word = strtok(NULL, " \t\n")) != NULL && nwords < maxwords) {
        words[nwords++] = word;
    }
    return nwords;
}

/* ---------------------------------------------------------------------- */

void CollideBGK::reset_relaxflag() {
    int nplocal = particle->nlocal;
    if (nplocal > nplocalmax) {
        nplocalmax = ceil(nplocal * 1.2);
        memory->destroy(relax_flag);
        memory->create(relax_flag, nplocalmax, "collide:relax_flag");
    }
    memset(relax_flag, 0, nplocalmax * sizeof(bool));
}

/* ---------------------------------------------------------------------- */

void CollideBGK::print_warning() {
    if (output->next_stats != update->ntimestep)return;
    bigint sum1, sum2, sum3, sum4;
    sum1 = sum2 = sum3 = sum4 = 0;
    MPI_Allreduce(&count_fail_relaxation, &sum1, 1, MPI_SPARTA_BIGINT, MPI_SUM, world);
    MPI_Allreduce(&count_done_relaxation, &sum2, 1, MPI_SPARTA_BIGINT, MPI_SUM, world);
    MPI_Allreduce(&count_warning_ignore_childcell, &sum3, 1, MPI_SPARTA_BIGINT, MPI_SUM, world);
    MPI_Allreduce(&count_do_childcell, &sum4, 1, MPI_SPARTA_BIGINT, MPI_SUM, world);
    if (comm->me == 0) {
        if (sum1) {
            char str[128];
            sprintf(str, "%d relaxation failed in total %d relaxation, percentage = %.4f",
                sum1, sum2, 100.0 * sum1 / sum2);
            error->warning(FLERR, str);
        }
        else if (sum3) {
            char str[128];
            sprintf(str, "%d cells ignored abnormally in total %d child cells, percentage = %.4f",
                sum3, sum3 + sum4, 100.0 * sum3 / (sum3 + sum4));
            error->warning(FLERR, str);
        }
    }
    reset_count();
}

/* ---------------------------------------------------------------------- */

CollideBGKModify::CollideBGKModify(SPARTA* sparta) : Pointers(sparta){}

/* ---------------------------------------------------------------------- */

CollideBGKModify::~CollideBGKModify(){}

/* ----------------------------------------------------------------------
* process collide_bgk_modify command, included in style_command.h
------------------------------------------------------------------------- */

void CollideBGKModify::command(int narg, char** arg)
{
    if (strcmp(collide->style, "bgk") != 0) {
        error->all(FLERR,
            "Using collide_bgk_modify command when collide.style != bgk");
    }
    if (narg == 0) error->all(FLERR, "Illegal collide_modify command");
    CollideBGK* collideBGK = dynamic_cast<CollideBGK*>(collide);
    if (!collideBGK) {
        error->all(FLERR, "CollideBGKModify: dynamic_cast fault");
    }
    int iarg = 0;
    while (iarg < narg) {
        if (strcmp(arg[iarg], "reset_wmax") == 0) {
            if (iarg + 2 > narg) error->all(FLERR, "Illegal collide_bgk_modify command");
            double reset = atof(arg[iarg + 1]);
            if (reset <= 0) {
                collideBGK->resetWmax = 0.0;
            }
            else if (reset >= 1.0) {
                error->all(FLERR, 
                    "Illegal collide_bgk_modify command: resetWmax > 1");
            }
            else {
                collideBGK->resetWmax = reset;
            }
            iarg += 2;
        }     
        else if (strcmp(arg[iarg], "pr_num") == 0) {
            if (iarg + 2 > narg) error->all(FLERR, "Illegal collide_bgk_modify command");
            collideBGK->Pr = atof(arg[iarg + 1]);
            if (collideBGK->Pr <= 0) 
                error->all(FLERR, "Illegal collide_bgk_modify Prantl number");
            iarg += 2;
        }
        else if (strcmp(arg[iarg], "time_ave") == 0) {
            if (iarg + 2 > narg) error->all(FLERR, "Illegal collide_bgk_modify command");
            collideBGK->time_ave_coef = atof(arg[iarg + 1]);
            if (collideBGK->time_ave_coef < 0 || collideBGK->time_ave_coef >=1)
                error->all(FLERR, "Illegal collide_bgk_modify time_ave_coef");
            iarg += 2;
        }else if (strcmp(arg[iarg], "interpolate") == 0) {
            if (iarg + 2 > narg) error->all(FLERR, "Illegal collide_bgk_modify command");
            if (strcmp(arg[iarg + 1], "yes") == 0) collideBGK->interpolate_flag = 1;
            else if (strcmp(arg[iarg + 1], "no") == 0) collideBGK->interpolate_flag = 0;
            else error->all(FLERR, "Illegal collide_bgk_modify command");
            iarg += 2;
        }
        else error->all(FLERR, "Illegal collide_bgk_modify command");

    }
}


