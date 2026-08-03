/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   polygon2d.cpp -- see polygon2d.h
------------------------------------------------------------------------- */

#include "polygon2d.h"

#include <cmath>

double polygon_signed_area(const double *pts, int npt)
{
  double area2 = 0.0;
  for (int i = 0; i < npt; i++) {
    int j = (i + 1) % npt;
    area2 += pts[2 * i] * pts[2 * j + 1] - pts[2 * j] * pts[2 * i + 1];
  }
  return 0.5 * area2;
}

void polygon_area_point_grad(const double *pts, int npt, double *dA_dpts)
{
  for (int k = 0; k < npt; k++) {
    int kp = (k + 1) % npt;
    int km = (k - 1 + npt) % npt;
    dA_dpts[2 * k]     = 0.5 * (pts[2 * kp + 1] - pts[2 * km + 1]);  // d/dx_k
    dA_dpts[2 * k + 1] = 0.5 * (pts[2 * km]     - pts[2 * kp]);      // d/dy_k
  }
}

double polygon_min_edge_length(const double *pts, int npt)
{
  double minlen = -1.0;
  for (int i = 0; i < npt; i++) {
    int j = (i + 1) % npt;
    double dx = pts[2 * j]     - pts[2 * i];
    double dy = pts[2 * j + 1] - pts[2 * i + 1];
    double len = std::sqrt(dx * dx + dy * dy);
    if (minlen < 0.0 || len < minlen) minlen = len;
  }
  return minlen;
}
