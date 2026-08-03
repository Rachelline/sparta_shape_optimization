/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer

   constraint.h: the Constraint interface. A general inequality constraint
   g_lo <= g(alpha) <= g_hi imposed on the optimizer's design vector,
   generic over Shape (uses only to_mesh()/jacobian(), so any shape
   family -- 2D or 3D -- gets constraints like MinSizeConstraint for
   free). chord/nseg (or nx/ntheta) are NOT passed here: they live inside
   each Shape's own constructor now, not threaded through every call.

   Constraints are pure geometry (no SPARTA/gas physics, no SPARTA run
   needed to evaluate one) -- contrast with Objective, which genuinely
   needs a full DSMC run.
------------------------------------------------------------------------- */

#ifndef SPARTA_CONSTRAINT_H
#define SPARTA_CONSTRAINT_H

#include "shape.h"

class Constraint {
 public:
  virtual ~Constraint() {}

  // g(alpha) -- pure geometry, no SPARTA run
  virtual double eval(const Shape &shape, const double *alpha) const = 0;

  // dg/d(alpha), length shape.ndesign() -- exact/analytic where possible
  // (see min_size_constraint.cpp), not FD/AD
  virtual void grad(const Shape &shape, const double *alpha,
                    double *grad) const = 0;

  virtual double lower() const = 0;
  virtual double upper() const = 0;

  // for config.txt / logging
  virtual const char *name() const = 0;
};

#endif
