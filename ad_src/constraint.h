/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer

   constraint.h: the Constraint interface. A general inequality constraint
   g_lo <= g(alpha) <= g_hi imposed on the optimizer's design vector,
   generic over Parametrization (uses only to_lines()/jacobian(), so any
   shape family gets constraints like MinAreaConstraint for free).

   Deliberately takes `chord`/`nseg` directly, not a whole RunConfig --
   constraints here are pure geometry (no SPARTA/gas physics, no SPARTA
   run needed to evaluate one), and the interface should make that
   cheapness visible rather than implying coupling that doesn't exist.
   Contrast with Objective, which genuinely needs a full DSMC run.
------------------------------------------------------------------------- */

#ifndef SPARTA_CONSTRAINT_H
#define SPARTA_CONSTRAINT_H

#include "parametrization.h"

class Constraint {
 public:
  virtual ~Constraint() {}

  // g(alpha) -- pure geometry, no SPARTA run
  virtual double eval(const Parametrization &shape, const double *alpha,
                      double chord, int nseg) const = 0;

  // dg/d(alpha), length shape.ndesign() -- exact/analytic where possible
  // (see min_area_constraint.cpp), not FD/AD
  virtual void grad(const Parametrization &shape, const double *alpha,
                    double chord, int nseg, double *grad) const = 0;

  virtual double lower() const = 0;
  virtual double upper() const = 0;

  // for config.txt / logging
  virtual const char *name() const = 0;
};

#endif
