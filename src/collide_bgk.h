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

#ifdef COLLIDE_CLASS

CollideStyle(bgk,CollideBGK)

#else
#ifdef COMMAND_CLASS

CommandStyle(collide_bgk_modify, CollideBGKModify)

#else

#ifndef SPARTA_COLLIDE_BGK_H
#define SPARTA_COLLIDE_BGK_H

#include "collide.h"
#include "particle.h"

namespace SPARTA_NS {

class CollideBGK : public Collide {
 public:
  friend class CollideBGKModify;
  CollideBGK(class SPARTA *, int, char **);
  virtual ~CollideBGK();
  virtual void init();
  virtual void collisions();

  double vremax_init(int, int);
  virtual double attempt_collision(int, int, double);
  double attempt_collision(int, int, int, double) { return 0.0; };
  virtual int test_collision(int, int, int, Particle::OnePart*, Particle::OnePart*) { return 1; };
  virtual void setup_collision(Particle::OnePart*, Particle::OnePart*) { return; };
  virtual int perform_collision(Particle::OnePart *&, Particle::OnePart *&,
                        Particle::OnePart *&);
  void perform_uspbgk(Particle::OnePart*, int, const class CommMacro*);
  void perform_bgkbgk(Particle::OnePart*, int, const class CommMacro*);
  void perform_esbgk(Particle::OnePart*, int, const class CommMacro*);
  void perform_sbgk(Particle::OnePart*, int, const class CommMacro*);
  void conservV();
  double extract(int, int, const char*) { return 0.0; };

  struct Params {             // BGK model parameters
      double mu_ref;          // reference viscosity
      double omega;           // mu ~ T^omega
      double T_ref;           // reference temperature
  };
  struct ConservMacro
  {
      int done_relaxation;
      double coef;
      double v_origin[3], v_post[3];
  };

  // status
  bigint count_try_relaxation, count_done_relaxation, count_fail_relaxation;
  bigint count_do_childcell, count_ignore_childcell, count_warning_ignore_childcell;

 protected:
  int nmaxconserv;
  ConservMacro* conservMacro;
  
  Params *params;             // BGK params for each species
  int nparams;                // # of per-species params read in
  int maxglocal;
  int interpolate_flag;       // 1/0 = yes/no do interpolation, default = 1;
  int *resetWmax_flag;        // flag to decide whether resetWmax in current cell
  double resetWmax;           // coefficient to reduce Wmax, default = 0.9999,
                              // if resetWmax <= 0, don't do reset
  int bgk_mod;
  double Pr;                  // Prantl number
  double time_ave_coef;

  bool* relax_flag;
  int nplocalmax;

  template < int > void computeMacro();
  void read_param_file(char*);
  int wordparse(int, char*, char**);
  void reset_count();
  void print_warning();

  void reset_relaxflag();

};

class CollideBGKModify : protected Pointers {
public:
    CollideBGKModify(class SPARTA*);
    ~CollideBGKModify();
    void command(int, char**);
};

}

#endif
#endif
#endif
