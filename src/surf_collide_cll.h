/* AD-CONVERTED: double->sfloat by tools/ad_convert.py (see sfloat.h) */
/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.github.io
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#ifdef SURF_COLLIDE_CLASS

SurfCollideStyle(cll,SurfCollideCLL)

#else

#ifndef SPARTA_SURF_COLLIDE_CLL_H
#define SPARTA_SURF_COLLIDE_CLL_H

#include "pointers.h"
#include "surf_collide.h"

namespace SPARTA_NS {

class SurfCollideCLL : public SurfCollide {
 public:
  SurfCollideCLL(class SPARTA *, int, char **);
  ~SurfCollideCLL();
  void init();
  Particle::OnePart *collide(Particle::OnePart *&, sfloat &,
                             int, sfloat *, int, int &);
  void wrapper(Particle::OnePart *, sfloat *, int *, sfloat*);
  void flags_and_coeffs(int *, sfloat *);

 private:
  sfloat acc_n,acc_t,acc_rot,acc_vib;   // surface accomodation coeffs
  sfloat vx,vy,vz;                      // translational velocity of surface
  sfloat wx,wy,wz;                      // angular velocity of surface
  sfloat px,py,pz;                      // point to rotate surface around
  sfloat eccen;                         // 1 if fully diffuse scattering
                                        // < 1 if partial diffuse scattering

  int tflag,rflag;           // flags for translation and rotation
  int trflag;                // 1 if either tflag or rflag is set
  int pflag;                 // 1 if partially energy accommodation
                             // with partial/fully diffuse scattering

  sfloat vstream[3];
  class RanKnuth *random;     // RNG for particle reflection

  void cll(Particle::OnePart *, sfloat *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

*/
