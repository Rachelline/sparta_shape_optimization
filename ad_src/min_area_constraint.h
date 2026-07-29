/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer

   min_area_constraint: enforces |signed_area(alpha)| >= min_area, i.e. a
   floor on the 2D cross-sectional area of the body -- the concrete
   answer to "there needs to be a minimum volume that prevents the Bezier
   curve from looking like a needle." ("volume" is area in the current
   2D-only Parametrization; a 3D revolve/extrude shape family would want
   its own volume constraint, not this one, when one exists.)

   Generic over Parametrization (shoelace formula on whatever to_lines()
   returns, not hardcoded to Bezier) -- works for any future shape family
   too, unlike svg_shape which stays Bezier-specific on purpose.
------------------------------------------------------------------------- */

#ifndef SPARTA_MIN_AREA_CONSTRAINT_H
#define SPARTA_MIN_AREA_CONSTRAINT_H

#include "constraint.h"

class MinAreaConstraint : public Constraint {
 public:
  explicit MinAreaConstraint(double min_area) : min_area_(min_area) {}

  double eval(const Parametrization &shape, const double *alpha,
             double chord, int nseg) const override;

  // Exact analytic gradient -- no FD/AD. Shoelace area is linear in each
  // point (d(Area)/d(x_k) = 0.5*(y_{k+1}-y_{k-1}), d(Area)/d(y_k) =
  // 0.5*(x_{k-1}-x_{k+1})), chained through shape.jacobian() (d(pts)/d(alpha),
  // already exact) to get d(Area)/d(alpha).
  void grad(const Parametrization &shape, const double *alpha,
           double chord, int nseg, double *grad) const override;

  double lower() const override { return min_area_; }
  // IPOPT's own "no upper bound" sentinel (its nlp_upper_bound_inf default)
  double upper() const override { return 1e19; }

  const char *name() const override { return "min_area"; }

 private:
  double min_area_;
};

#endif
