/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   min_area_constraint.cpp -- see min_area_constraint.h
------------------------------------------------------------------------- */

#include "min_area_constraint.h"

#include <cmath>
#include <vector>

namespace {

// Shoelace signed area over the first `npt` points of a closed loop
// (pts[npt] == pts[0], the duplicate closing point, is not touched).
double shoelace_signed_area(const double *pts, int npt)
{
  double area2 = 0.0;
  for (int i = 0; i < npt; i++) {
    int j = (i + 1) % npt;
    area2 += pts[2 * i] * pts[2 * j + 1] - pts[2 * j] * pts[2 * i + 1];
  }
  return 0.5 * area2;
}

// d(signed_area)/d(pts), length 2*npt: d/d(x_k) = 0.5*(y_{k+1}-y_{k-1}),
// d/d(y_k) = 0.5*(x_{k-1}-x_{k+1}), indices wrapping mod npt.
void shoelace_area_point_grad(const double *pts, int npt, double *dA_dpts)
{
  for (int k = 0; k < npt; k++) {
    int kp = (k + 1) % npt;
    int km = (k - 1 + npt) % npt;
    dA_dpts[2 * k]     = 0.5 * (pts[2 * kp + 1] - pts[2 * km + 1]);  // d/dx_k
    dA_dpts[2 * k + 1] = 0.5 * (pts[2 * km]     - pts[2 * kp]);      // d/dy_k
  }
}

}  // namespace

double MinAreaConstraint::eval(const Parametrization &shape, const double *alpha,
                               double chord, int nseg) const
{
  int npt = shape.nsegments(nseg);
  std::vector<double> pts(2 * (npt + 1)), norms(2 * npt);
  shape.to_lines(alpha, chord, nseg, pts.data(), norms.data());
  return std::fabs(shoelace_signed_area(pts.data(), npt));
}

void MinAreaConstraint::grad(const Parametrization &shape, const double *alpha,
                             double chord, int nseg, double *grad) const
{
  int npt = shape.nsegments(nseg);
  int ndesign = shape.ndesign();

  std::vector<double> pts(2 * (npt + 1)), norms(2 * npt);
  shape.to_lines(alpha, chord, nseg, pts.data(), norms.data());

  double signed_area = shoelace_signed_area(pts.data(), npt);
  double sign = (signed_area >= 0.0) ? 1.0 : -1.0;   // d|A|/dA

  std::vector<double> dA_dpts(2 * npt);
  shoelace_area_point_grad(pts.data(), npt, dA_dpts.data());

  std::vector<double> jac(2 * (npt + 1) * ndesign);
  shape.jacobian(nseg, jac.data());

  for (int c = 0; c < ndesign; c++) {
    double sum = 0.0;
    for (int r = 0; r < 2 * npt; r++) sum += dA_dpts[r] * jac[r * ndesign + c];
    grad[c] = sign * sum;
  }
}
