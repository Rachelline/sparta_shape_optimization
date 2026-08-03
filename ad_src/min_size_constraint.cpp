/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   min_size_constraint.cpp -- see min_size_constraint.h
------------------------------------------------------------------------- */

#include "min_size_constraint.h"
#include "polygon2d.h"

#include <cmath>
#include <vector>

namespace {

double cross_dot(const double a[3], const double b[3], const double c[3])
{
  // a . (b x c)
  return a[0] * (b[1] * c[2] - b[2] * c[1]) -
         a[1] * (b[0] * c[2] - b[2] * c[0]) +
         a[2] * (b[0] * c[1] - b[1] * c[0]);
}

void cross(const double a[3], const double b[3], double out[3])
{
  out[0] = a[1] * b[2] - a[2] * b[1];
  out[1] = a[2] * b[0] - a[0] * b[2];
  out[2] = a[0] * b[1] - a[1] * b[0];
}

double signed_volume6(const SurfMesh &m)
{
  double v6 = 0.0;
  int ntri = m.nelems();
  for (int t = 0; t < ntri; t++) {
    const double *p0 = &m.pts[3 * m.elems[3 * t]];
    const double *p1 = &m.pts[3 * m.elems[3 * t + 1]];
    const double *p2 = &m.pts[3 * m.elems[3 * t + 2]];
    v6 += cross_dot(p0, p1, p2);
  }
  return v6;
}

}  // namespace

double MinSizeConstraint::eval(const Shape &shape, const double *alpha) const
{
  SurfMesh m;
  shape.to_mesh(alpha, m);

  if (m.dim == 2) {
    return std::fabs(polygon_signed_area(m.pts.data(), m.npoints()));
  } else {
    return std::fabs(signed_volume6(m) / 6.0);
  }
}

void MinSizeConstraint::grad(const Shape &shape, const double *alpha,
                             double *grad) const
{
  SurfMesh m;
  shape.to_mesh(alpha, m);
  int ndesign = shape.ndesign();

  std::vector<double> jac;
  shape.jacobian(alpha, jac);

  if (m.dim == 2) {
    int npt = m.npoints();
    double signed_area = polygon_signed_area(m.pts.data(), npt);
    double sign = (signed_area >= 0.0) ? 1.0 : -1.0;   // d|A|/dA

    std::vector<double> dA_dpts(2 * npt);
    polygon_area_point_grad(m.pts.data(), npt, dA_dpts.data());

    for (int c = 0; c < ndesign; c++) {
      double sum = 0.0;
      for (int r = 0; r < 2 * npt; r++) sum += dA_dpts[r] * jac[r * ndesign + c];
      grad[c] = sign * sum;
    }
    return;
  }

  // 3D: d(6V)/d(alpha_k) = sum_tri (p1 x p2).dp0/dalpha_k
  //                               + (p2 x p0).dp1/dalpha_k
  //                               + (p0 x p1).dp2/dalpha_k
  double v6 = signed_volume6(m);
  double sign = (v6 >= 0.0) ? 1.0 : -1.0;   // d|V|/dV, V = v6/6

  for (int c = 0; c < ndesign; c++) grad[c] = 0.0;

  int ntri = m.nelems();
  for (int t = 0; t < ntri; t++) {
    int i0 = m.elems[3 * t], i1 = m.elems[3 * t + 1], i2 = m.elems[3 * t + 2];
    const double *p0 = &m.pts[3 * i0];
    const double *p1 = &m.pts[3 * i1];
    const double *p2 = &m.pts[3 * i2];

    double a[3], b[3], c3[3];
    cross(p1, p2, a);   // p1 x p2
    cross(p2, p0, b);   // p2 x p0
    cross(p0, p1, c3);  // p0 x p1

    for (int k = 0; k < ndesign; k++) {
      double contrib = 0.0;
      for (int d = 0; d < 3; d++) {
        contrib += a[d]  * jac[(3 * i0 + d) * ndesign + k];
        contrib += b[d]  * jac[(3 * i1 + d) * ndesign + k];
        contrib += c3[d] * jac[(3 * i2 + d) * ndesign + k];
      }
      grad[k] += contrib;
    }
  }

  for (int k = 0; k < ndesign; k++) grad[k] = sign * grad[k] / 6.0;
}
