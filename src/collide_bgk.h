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

#ifdef COLLIDE_CLASS

CollideStyle(bgk,CollideBGK)

#else

#ifndef SPARTA_COLLIDE_BGK_H
#define SPARTA_COLLIDE_BGK_H

#include "collide.h"
#include "particle.h"

namespace SPARTA_NS {

class CollideBGK : public Collide {
 public:
  CollideBGK(class SPARTA *, int, char **);
  virtual ~CollideBGK();
  virtual void init();
  virtual void collisions();

  double vremax_init(int, int);
  virtual double attempt_collision(int, int, double);
  double attempt_collision(int, int, int, double);
  virtual int test_collision(int, int, int, Particle::OnePart*, Particle::OnePart*) { return 1; };
  virtual void setup_collision(Particle::OnePart*, Particle::OnePart*) { return; };
  virtual int perform_collision(Particle::OnePart *&, Particle::OnePart *&,
                        Particle::OnePart *&);
  virtual void collisions_BGK();
  void perform_uspbgk(Particle::OnePart*, int, const class CommMacro*);
  void perform_bgkbgk(Particle::OnePart*, int, const class CommMacro*);
  void perform_esbgk(Particle::OnePart*, int, const class CommMacro*);
  void perform_sbgk(Particle::OnePart*, int, const class CommMacro*);
  void conservV();
  double extract(int, int, const char *);

  struct Params {             // BGK model parameters
      double diam;
      double omega;
      double tref;
      double alpha;
  };
  struct ConservMacro
  {
      double coef;
      double v_origin[3], v_post[3];
  };

 protected:
  int nmaxconserv;
  ConservMacro* conservMacro;
  //double* prefactor;          // static portion of collision attempt frequency
  Params *params;             // BGK params for each species
  double resetWmax;          // coefficient to reduce Wmax, default = 0.9999, only reset if >0
  int nparams;                // # of per-species params read in
  int bgk_mod;
  double time_ave_coef;
  template < int > void computeMacro();
  void read_param_file(char*);
  int wordparse(int, char*, char**);

};

}

#endif
#endif
