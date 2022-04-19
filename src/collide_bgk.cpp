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

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "collide_vss.h"
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

using namespace SPARTA_NS;
using namespace MathConst;

#define MAXLINE 1024
enum { USP, BGK, ESBGK, SBGK };
/* ---------------------------------------------------------------------- */

CollideBGK::CollideBGK(SPARTA* sparta, int narg, char** arg) :
    Collide(sparta, narg, arg)
{
    if (narg != 4) error->all(FLERR, "Illegal collide command");

    // proc 0 reads file to extract params for current species
    // broadcasts params to all procs
    time_ave_coef = 0.99;
    nparams = particle->nspecies;
    resetWmax = 0.9999;

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

    if (comm->me == 0) read_param_file(arg[3]);
    MPI_Bcast(params, nparams * sizeof(Params), MPI_BYTE, 0, world);

    // allocate per-species prefactor array
    //memory->create(prefactor, nparams, "collide_bgk:prefactor");
    //for (int i = 0; i < particle->nspecies; ++i) {
    //    Particle::Species& species = particle->species[i];
    //    prefactor[i] = update->boltz * pow(species.Tref, species.omega) / species.muref;
    //}
}

/* ---------------------------------------------------------------------- */

CollideBGK::~CollideBGK()
{
    if (copymode) return;

    memory->destroy(params);
    //memory->destroy(prefactor);
}

void CollideBGK::init()
{
    // initially read-in per-species params must match current species list

    if (nparams != particle->nspecies)
        error->all(FLERR, "BGK parameters do not match current species");

    Collide::init();
}

void CollideBGK::collisions()
{
    int i, j, k, m, n, ip, np;
    int nattempt, reactflag;
    double attempt, volume;
    Particle::OnePart* ipart, * jpart, * kpart;

    // loop over cells I own

    Grid::ChildInfo* cinfo = grid->cinfo;

    Particle::OnePart* particles = particle->particles;
    int* next = particle->next;

    for (int icell = 0; icell < nglocal; icell++) {
        np = cinfo[icell].count;
        if (np <= 3) continue;
        ip = cinfo[icell].first;
        volume = cinfo[icell].volume / cinfo[icell].weight;
        if (volume == 0.0) error->one(FLERR, "Collision cell volume is zero");

        // setup particle list for this cell

        if (np > npmax) {
            while (np > npmax) npmax += DELTAPART;
            memory->destroy(plist);
            memory->create(plist, npmax, "collide:plist");
        }

        if (bgk_mod == USP) computeMacro<USP>();
        else if (bgk_mod == BGK) computeMacro<BGK>();
        else if (bgk_mod == SBGK) computeMacro<SBGK>();
        else if (bgk_mod == ESBGK) computeMacro<ESBGK>();

        Grid::ChildCell* cells = grid->cells;
        double bgk_attempt = attempt_collision(icell,0,cinfo[icell].macro.tao);
        int bgk_nattempt = static_cast<int> (bgk_attempt + (random->uniform()));

        n = 0;
        while (ip >= 0) {
            plist[n++] = ip;
            ip = next[ip];
        }
        if (bgk_nattempt < np / 2) {
            for (i = 0; i < bgk_nattempt; i++) {
                int t = i + (np - i) * random->uniform();
                std::swap(plist[i], plist[t]);
            }
        }
        else {
            for (i = np - 1; i > bgk_nattempt - 1; i--) {
                int t = i * random->uniform();
                std::swap(plist[i], plist[t]);
            }
        }
        for (i = 0;i < bgk_nattempt; ++i) {        
            ipart = &particles[plist[i]];
            const CommMacro* interMacro = grid->gridCommMacro->interpolation(ipart);
            if ((!interMacro) || (!(interMacro->Temp > 0))) {
                if (!interMacro)
                    error->warning(FLERR, "CollideBGK:interpolation failed!(!interMacro)");
                interMacro = &grid->cells[icell].macro;
            }
            if (bgk_mod == USP) perform_uspbgk(ipart, icell, interMacro);
            else if (bgk_mod == BGK) perform_bgkbgk(ipart, icell, interMacro);
            else if (bgk_mod == SBGK) perform_sbgk(ipart, icell, interMacro);
            else if (bgk_mod == ESBGK) perform_esbgk(ipart, icell, interMacro);
        }   
    }

}
void CollideBGK::conservV() {
    int nlocal = grid->nlocal;
    if (nlocal > nmaxconserv) {
        while (nlocal > nmaxconserv) nmaxconserv += DELTAPART;
        memory->destroy(conservMacro);
        memory->create(conservMacro, nmaxconserv, "collideBGK:postmacro");
    }
    ConservMacro tmpmacro{ 0.0 };
    for (int i = 0; i < nlocal; ++i) {
        memcpy(conservMacro + i, &tmpmacro, sizeof(ConservMacro));
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
        NoCommMacro& nmacro = grid->cinfo[icell].macro;
        double np = grid->cinfo[icell].count;
        for (int i = 0; i < 3; ++i) {
            conservMacro[icell].v_origin[i] = grid->cells[icell].macro.v[i];
            conservMacro[icell].v_post[i] = nmacro.sum_vi[i]/np;
        }
        double theta_post = (nmacro.sum_vij[0]
            - (nmacro.sum_vi[0] * nmacro.sum_vi[0] + nmacro.sum_vi[1] * nmacro.sum_vi[1]
                + nmacro.sum_vi[2] * nmacro.sum_vi[2]) / np) / np / 3;
        conservMacro[icell].coef = sqrt(grid->cells[icell].macro.theta / theta_post);
    }
    for (int ipart = 0; ipart < particle->nlocal; ++ipart) {
        Particle::OnePart& part = particle->particles[ipart];
        ConservMacro& cm = conservMacro[part.icell];
        for (int i = 0; i < 3; ++i) {
            part.v[i] = (part.v[i] - cm.v_post[i]) * cm.coef + cm.v_origin[i];
        }
    }
}

void CollideBGK::perform_uspbgk(Particle::OnePart* ip, int icell, const CommMacro* interMacro)
{
    Grid::ChildInfo* cinfo = grid->cinfo;
    const double* pij = cinfo[icell].macro.sigma_ij;
    const double* sigma_ij = cinfo[icell].macro.sigma_ij;
    const double* q = cinfo[icell].macro.qi;
    double vn[3];
    while (true)
    {
        for (int i = 0; i < 3; i++) vn[i] = random->gaussian() * sqrt(interMacro->theta);
        double C_2 = vn[0] * vn[0] + vn[1] * vn[1] + vn[2] * vn[2];
        double sigmacc = sigma_ij[0] * (vn[0] * vn[0] - C_2 / 3)
            + sigma_ij[1] * (vn[1] * vn[1] - C_2 / 3)
            + sigma_ij[2] * (vn[2] * vn[2] - C_2 / 3)
            + 2 * sigma_ij[3] * vn[0] * vn[1]
            + 2 * sigma_ij[4] * vn[0] * vn[2]
            + 2 * sigma_ij[5] * vn[2] * vn[1];
        double qkck = (vn[0] * q[0] + vn[1] * q[1] + vn[2] * q[2]) *
            ((vn[0] * vn[0] + vn[1] * vn[1] + vn[2] * vn[2])
                / interMacro->theta - 5);

        double W = 1.0 + cinfo[icell].macro.coef_A * sigmacc +
            cinfo[icell].macro.coef_B * qkck;
        if (random->uniform() < W / cinfo[icell].macro.Wmax) {
            if (W > cinfo[icell].macro.Wmax) cinfo[icell].macro.Wmax = W;
            else if (resetWmax > 0.0) cinfo[icell].macro.Wmax *= resetWmax;
            break;
        }
    }
    for (int i = 0; i < 3; i++) ip->v[i] = vn[i] + interMacro->v[i];

    
}

/* ----------------------------------------------------------------------
   estimate a good value for vremax for a group pair in any grid cell
   called by Collide parent in init()

   NOTE: vremax is useless in BGK model
         
------------------------------------------------------------------------- */

double CollideBGK::vremax_init(int igroup, int jgroup)
{
    return 0.0;
}

/* ---------------------------------------------------------------------- */

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
    if (bgk_mod == ESBGK) bgk_nattempt = update->Pr * np * (1 - exp(-tao));
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
        if (MOD == USP || MOD == ESBGK) {
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
        if (np <= 3) continue;
        // Currently assume all particles have same ispecies
        Particle::Species& 
            species = particle->species[particles[cinfo.first].ispecies];
        double factor = species.mass * update->fnum * cinfo.weight / cinfo.volume;
        double pij[6]{}, qi[3]{};
        double* sum_vij = mean_nmacro.sum_vij;
        double* v = cmacro.v;
        for (int i = 0; i < 3; ++i) {
            v[i] = mean_nmacro.sum_vi[i] / np;
        }
        double V_2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
        double sum_C2 = sum_vij[0] + sum_vij[1] + sum_vij[2];
        cmacro.theta = (sum_C2 / np - V_2) / 3;
        cmacro.Temp = species.mass / update->boltz * cmacro.theta;
        double nrho = cinfo.count * update->fnum * cinfo.weight / cinfo.volume;
        mean_nmacro.tao = nrho * update->boltz * pow(species.Tref, species.omega)
            * pow(cmacro.Temp, 1 - species.omega) * update->dt / species.muref;
        double p = 0.0;
        if (MOD == USP || MOD == ESBGK) {
            for (int i = 0; i < 3; ++i) {
                pij[i] = factor * (sum_vij[i] - np * v[i] * v[i]);
            }
            pij[3] = factor * (sum_vij[3] - np * v[0] * v[1]);
            pij[4] = factor * (sum_vij[4] - np * v[0] * v[2]);
            pij[5] = factor * (sum_vij[5] - np * v[1] * v[2]);
            // time-average pij
            p = (pij[0] + pij[1] + pij[2]) / 3.0;
            for (int i = 0; i < 6; ++i) {
                mean_nmacro.sigma_ij[i] = mean_nmacro.sigma_ij[i] * time_ave_coef
                    + (pij[i] - p * (i < 3 ? 1.0 : 0.0)) * (1 - time_ave_coef);
            }
        }        
        if (MOD == USP || MOD == SBGK) {
            qi[0] = factor / 2 * (mean_nmacro.sum_C2vi[0]
                - v[0] * sum_C2 + 2 * np * V_2 * v[0]
                - 2 * (v[0] * sum_vij[0] + v[1] * sum_vij[3] + v[2] * sum_vij[4]));
            qi[1] = factor / 2 * (mean_nmacro.sum_C2vi[1]
                - v[1] * sum_C2 + 2 * np * V_2 * v[1]
                - 2 * (v[0] * sum_vij[3] + v[1] * sum_vij[1] + v[2] * sum_vij[5]));
            qi[2] = factor / 2 * (mean_nmacro.sum_C2vi[2]
                - v[2] * sum_C2 + 2 * np * V_2 * v[2]
                - 2 * (v[0] * sum_vij[4] + v[1] * sum_vij[5] + v[2] * sum_vij[2]));

            // time-average qi
            for (int i = 0; i < 3; ++i) {
                mean_nmacro.qi[i] = mean_nmacro.qi[i] * time_ave_coef
                    + qi[i] * (1.0 - time_ave_coef);
            }
        }
        if (MOD == USP) {
            double p_theta = p * cmacro.theta;
            double tao_coth = mean_nmacro.tao / 2 * (1 + 2 / (exp(mean_nmacro.tao) + 1));
            mean_nmacro.coef_A = (1 - tao_coth) / (2 * p_theta);
            mean_nmacro.coef_B = (1 - update->Pr * tao_coth) / (5 * p_theta);
        }
    }

    for (int isc = 0; isc < surf->nsc; ++isc) {
        double mass = particle->species[particle->particles[0].ispecies].mass;
        if (strcmp(surf->sc[isc]->style, "diffuse") != 0) continue;
        CommMacro* cmacro = surf->sc[isc]->returnComm();
        if (cmacro && cmacro->theta < 0 && cmacro->Temp>0) {

            cmacro->theta == cmacro->Temp / mass * update->boltz;
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
        params[i].diam = -1.0;
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
            if (nwords < REQWORDS + 1)  // one extra word in cross-species lines
                error->one(FLERR, "Incorrect line format in VSS parameter file");
            params[isp].diam  = atof(words[1]);
            params[isp].omega = atof(words[2]);
            params[isp].tref  = atof(words[3]);
            params[isp].alpha = atof(words[4]);
        }
    }

    delete[] words;
    fclose(fp);

    // check that params were read for all species
    for (int i = 0; i < nparams; i++) {

        if (params[i].diam < 0.0) {
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

