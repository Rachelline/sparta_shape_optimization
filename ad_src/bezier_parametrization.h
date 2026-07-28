/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer

   bezier_parametrization.h: thin Parametrization adapter wrapping
   BezierGeom's symmetric-body free functions (bezier_geom.h) 1:1.

   ndesign() == 4: alpha = [x1, y1, x2, y2], the two interior control
   points of the upper-half cubic Bezier (see bezier_geom.h for the full
   geometry convention -- nose at origin, tail at (chord,0), lower half
   mirrored). `nseg` throughout is segments PER HALF, matching
   BezierGeom's own convention: to_lines()/jacobian() produce 2*nseg
   total segments and 2*nseg+1 points (pts buffer length
   2*(2*nseg+1) doubles, norms buffer length 2*(2*nseg) doubles).
------------------------------------------------------------------------- */

#ifndef SPARTA_BEZIER_PARAMETRIZATION_H
#define SPARTA_BEZIER_PARAMETRIZATION_H

#include "parametrization.h"

class BezierParametrization : public Parametrization {
 public:
  int ndesign() const override { return 4; }

  int nsegments(int nseg) const override { return 2 * nseg; }

  void to_lines(const double *alpha, double chord, int nseg,
                double *pts, double *norms) const override;

  void jacobian(int nseg, double *jac) const override;

  bool validate(int nseg, const double *pts, std::string *why) const override;

  // y in [0.05, 3.0]: matches the reference optimizer's defaults, and
  // is load-bearing (y1,y2 > 0 is required for a clockwise/valid body,
  // see validate()). x bounds are a provisional [0.2, 3.8], matching
  // the reference's own worked example at chord=4.0 -- this interface
  // doesn't pass chord to bounds(), so a real x-upper-bound of
  // chord-0.2 needs to be resolved by whatever eventually constructs a
  // ShapeTNLP with a non-default chord, not by this class.
  void bounds(double *lo, double *hi) const override;
};

#endif
