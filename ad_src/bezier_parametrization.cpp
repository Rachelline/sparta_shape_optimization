/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   bezier_parametrization.cpp -- see bezier_parametrization.h
------------------------------------------------------------------------- */

#include "bezier_parametrization.h"
#include "bezier_geom.h"

void BezierParametrization::to_lines(const double *alpha, double chord,
                                     int nseg, double *pts,
                                     double *norms) const
{
  BezierGeom::symmetric_body_to_lines(alpha, chord, nseg, pts, norms);
}

void BezierParametrization::jacobian(int nseg, double *jac) const
{
  BezierGeom::symmetric_body_jacobian(nseg, jac);
}

bool BezierParametrization::validate(int nseg, const double *pts,
                                     std::string *why) const
{
  int ntot = 2 * nseg;
  if (BezierGeom::signed_area(ntot, pts) >= 0.0) {
    if (why) *why = "body not clockwise (need y1,y2 > 0 for a valid "
                    "symmetric body)";
    return false;
  }
  if (BezierGeom::min_segment_length(ntot, pts) <= 0.0) {
    if (why) *why = "degenerate (zero-length) segment";
    return false;
  }
  return true;
}

void BezierParametrization::bounds(double *lo, double *hi) const
{
  // order matches alpha = [x1, y1, x2, y2]
  lo[0] = 0.2;  hi[0] = 3.8;
  lo[1] = 0.05; hi[1] = 3.0;
  lo[2] = 0.2;  hi[2] = 3.8;
  lo[3] = 0.05; hi[3] = 3.0;
}
