/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer

   bezier_shape.h: Shape adapter wrapping BezierGeom's symmetric-body
   free functions (bezier_geom.h) 1:1 -- replaces bezier_parametrization.h
   now that chord/nseg live in the constructor instead of being threaded
   through every call.

   ndesign() == 4: alpha = [x1, y1, x2, y2], the two interior control
   points of the upper-half cubic Bezier (see bezier_geom.h for the full
   geometry convention -- nose at origin, tail at (chord,0), lower half
   mirrored). nseg is segments PER HALF (BezierGeom's own convention):
   npoints() == 2*nseg unique points (no duplicate closing point --
   to_mesh() builds the explicit wraparound connectivity itself).
------------------------------------------------------------------------- */

#ifndef SPARTA_BEZIER_SHAPE_H
#define SPARTA_BEZIER_SHAPE_H

#include "shape.h"

class BezierShape : public Shape {
 public:
  BezierShape(double chord = 4.0, int nseg = 25) : chord_(chord), nseg_(nseg) {}

  int dim() const override { return 2; }
  int ndesign() const override { return 4; }
  int npoints() const { return 2 * nseg_; }

  // y in [0.05, 3.0]: load-bearing (y1,y2 > 0 is required for a
  // clockwise/valid body, see validate()). x bounds are a provisional
  // [0.2, 3.8], matching chord=4.0's worked example -- this interface
  // doesn't expose chord to callers of bounds(), so a real x-upper-bound
  // of chord-0.2 needs to be resolved by whoever constructs a BezierShape
  // with a non-default chord, not by this class.
  void bounds(double *lo, double *hi) const override;

  void to_mesh(const double *alpha, SurfMesh &m) const override;
  void jacobian(const double *alpha, std::vector<double> &jac) const override;
  bool validate(const SurfMesh &m, std::string *why) const override;

 private:
  double chord_;
  int nseg_;
};

#endif
