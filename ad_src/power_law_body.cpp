/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   power_law_body.cpp -- see power_law_body.h. Outward-normal convention
   verified analytically: for the parametric surface S(x,theta) =
   (x, r(x)cos(theta), r(x)sin(theta)), dS/dtheta x dS/dx points radially
   outward (away from the axis) for any r(x) >= 0, so every triangle below
   is built with edge order (+theta-ish, then +x-ish) from its first
   vertex, matching Surf::compute_tri_normal's (p2-p1)x(p3-p1) convention.
------------------------------------------------------------------------- */

#include "power_law_body.h"
#include <cmath>

namespace {
inline int idx_apex() { return 0; }
inline int idx_ring(int k, int j, int ntheta) { return 1 + (k - 1) * ntheta + j; }
inline int idx_base(int nx, int ntheta) { return 1 + nx * ntheta; }
}

void PowerLawBody::to_tris(const double *alpha, double *pts, int *tris) const
{
  double n = alpha[0];
  const double twopi = 2.0 * M_PI;

  pts[3 * idx_apex() + 0] = 0.0;
  pts[3 * idx_apex() + 1] = 0.0;
  pts[3 * idx_apex() + 2] = 0.0;

  for (int k = 1; k <= nx_; k++) {
    double xfrac = (double) k / (double) nx_;   // in (0,1]
    double x = L_ * xfrac;
    double r = Rmax_ * std::pow(xfrac, n);
    for (int j = 0; j < ntheta_; j++) {
      double theta = twopi * j / ntheta_;
      int p = idx_ring(k, j, ntheta_);
      pts[3 * p + 0] = x;
      pts[3 * p + 1] = r * std::cos(theta);
      pts[3 * p + 2] = r * std::sin(theta);
    }
  }

  int ib = idx_base(nx_, ntheta_);
  pts[3 * ib + 0] = L_;
  pts[3 * ib + 1] = 0.0;
  pts[3 * ib + 2] = 0.0;

  int t = 0;
  // apex fan
  for (int j = 0; j < ntheta_; j++) {
    int jp = (j + 1) % ntheta_;
    tris[3 * t + 0] = idx_apex();
    tris[3 * t + 1] = idx_ring(1, jp, ntheta_);
    tris[3 * t + 2] = idx_ring(1, j, ntheta_);
    t++;
  }
  // side quads (2 tris each)
  for (int k = 1; k < nx_; k++) {
    for (int j = 0; j < ntheta_; j++) {
      int jp = (j + 1) % ntheta_;
      tris[3 * t + 0] = idx_ring(k, j, ntheta_);
      tris[3 * t + 1] = idx_ring(k, jp, ntheta_);
      tris[3 * t + 2] = idx_ring(k + 1, j, ntheta_);
      t++;
      tris[3 * t + 0] = idx_ring(k, jp, ntheta_);
      tris[3 * t + 1] = idx_ring(k + 1, jp, ntheta_);
      tris[3 * t + 2] = idx_ring(k + 1, j, ntheta_);
      t++;
    }
  }
  // base fan
  for (int j = 0; j < ntheta_; j++) {
    int jp = (j + 1) % ntheta_;
    tris[3 * t + 0] = idx_ring(nx_, j, ntheta_);
    tris[3 * t + 1] = idx_ring(nx_, jp, ntheta_);
    tris[3 * t + 2] = idx_base(nx_, ntheta_);
    t++;
  }
}

void PowerLawBody::jacobian(const double *alpha, double *jac) const
{
  double n = alpha[0];
  const double twopi = 2.0 * M_PI;
  int npts = npoints();
  for (int i = 0; i < 3 * npts; i++) jac[i] = 0.0;   // apex, base fixed

  for (int k = 1; k <= nx_; k++) {
    double xfrac = (double) k / (double) nx_;
    double r = Rmax_ * std::pow(xfrac, n);
    double dr_dn = (xfrac > 0.0) ? r * std::log(xfrac) : 0.0;  // 0 at k==nx_ (xfrac==1)
    for (int j = 0; j < ntheta_; j++) {
      double theta = twopi * j / ntheta_;
      int p = idx_ring(k, j, ntheta_);
      jac[3 * p + 0] = 0.0;
      jac[3 * p + 1] = dr_dn * std::cos(theta);
      jac[3 * p + 2] = dr_dn * std::sin(theta);
    }
  }
}
