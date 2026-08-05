/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   bezier_shape.cpp -- see bezier_shape.h
------------------------------------------------------------------------- */

#include "bezier_shape.h"
#include "bezier_geom.h"
#include "polygon2d.h"

#include <vector>

void BezierShape::bounds(double *lo, double *hi) const
{
  // order matches alpha = [x1, y1, x2, y2]
  lo[0] = 0.2;  hi[0] = 3.8;
  lo[1] = 0.05; hi[1] = 3.0;
  lo[2] = 0.2;  hi[2] = 3.8;
  lo[3] = 0.05; hi[3] = 3.0;
}

void BezierShape::to_mesh(const double *alpha, SurfMesh &m) const
{
  int npt = npoints();
  std::vector<double> pts_closed(2 * (npt + 1)), norms(2 * npt);
  BezierGeom::symmetric_body_to_lines(alpha, chord_, nseg_, pts_closed.data(),
                                      norms.data());

  m.dim = 2;
  m.pts.assign(pts_closed.begin(), pts_closed.begin() + 2 * npt);
  m.elems.resize(2 * npt);
  for (int i = 0; i < npt; i++) {
    m.elems[2 * i]     = i;
    m.elems[2 * i + 1] = (i + 1) % npt;
  }
}

void BezierShape::jacobian(const double *alpha, std::vector<double> &jac) const
{
  (void) alpha;   // alpha-independent: symmetric_body_jacobian is a linear map
  int npt = npoints();
  std::vector<double> jac_closed(2 * (npt + 1) * ndesign());
  BezierGeom::symmetric_body_jacobian(nseg_, jac_closed.data());
  jac.assign(jac_closed.begin(), jac_closed.begin() + 2 * npt * ndesign());
}

bool BezierShape::validate(const SurfMesh &m, std::string *why) const
{
  int npt = m.npoints();
  if (polygon_signed_area(m.pts.data(), npt) >= 0.0) {
    if (why) *why = "body not clockwise (need y1,y2 > 0 for a valid "
                    "symmetric body)";
    return false;
  }
  if (polygon_min_edge_length(m.pts.data(), npt) <= 0.0) {
    if (why) *why = "degenerate (zero-length) segment";
    return false;
  }
  return true;
}
