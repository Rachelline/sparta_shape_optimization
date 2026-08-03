/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer

   min_size_constraint: enforces |measure(alpha)| >= min_size, i.e. a
   floor on the body's size -- the concrete answer to "there needs to be
   a minimum size that prevents the shape from collapsing to a needle
   chasing a smaller objective." measure() is the 2D cross-sectional
   area (shoelace) for a dim()==2 Shape, or the enclosed volume
   (divergence theorem) for a dim()==3 Shape. Replaces MinAreaConstraint,
   which only had a 2D-area notion of "size" (min_area_constraint.h
   itself conceded a 3D shape would need its own volume constraint --
   now it's the same class).

   Generic over Shape (uses only to_mesh()/jacobian()) -- works for any
   future shape family too, unlike svg_shape which stays Bezier-specific
   on purpose.
------------------------------------------------------------------------- */

#ifndef SPARTA_MIN_SIZE_CONSTRAINT_H
#define SPARTA_MIN_SIZE_CONSTRAINT_H

#include "constraint.h"

class MinSizeConstraint : public Constraint {
 public:
  explicit MinSizeConstraint(double min_size) : min_size_(min_size) {}

  // dim()==2: |shoelace area|. dim()==3: |enclosed volume|
  // (V = (1/6) sum_tri p0.(p1 x p2), valid for any closed, consistently
  // outward-wound triangle mesh regardless of choice of origin).
  double eval(const Shape &shape, const double *alpha) const override;

  // Exact analytic gradient -- no FD/AD, chained through shape.jacobian()
  // (d(pts)/d(alpha), already exact) in both branches. 2D: shoelace area
  // is linear in each point (d(Area)/d(x_k) = 0.5*(y_{k+1}-y_{k-1}),
  // d(Area)/d(y_k) = 0.5*(x_{k-1}-x_{k+1})). 3D: d(6V)/d(alpha) sums, over
  // every triangle (p0,p1,p2), (p1 x p2).dp0/dalpha + (p2 x p0).dp1/dalpha
  // + (p0 x p1).dp2/dalpha (cyclic permutations of the scalar triple
  // product's own linearity in each vertex).
  void grad(const Shape &shape, const double *alpha, double *grad) const override;

  double lower() const override { return min_size_; }
  // IPOPT's own "no upper bound" sentinel (its nlp_upper_bound_inf default)
  double upper() const override { return 1e19; }

  const char *name() const override { return "min_size"; }

 private:
  double min_size_;
};

#endif
