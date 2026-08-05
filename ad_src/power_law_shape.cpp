/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   power_law_shape.cpp -- see power_law_shape.h
------------------------------------------------------------------------- */

#include "power_law_shape.h"

void PowerLawShape::to_mesh(const double *alpha, SurfMesh &m) const
{
  m.dim = 3;
  m.pts.resize(3 * body_.npoints());
  m.elems.resize(3 * body_.ntris());
  body_.to_tris(alpha, m.pts.data(), m.elems.data());
}

void PowerLawShape::jacobian(const double *alpha, std::vector<double> &jac) const
{
  jac.resize(3 * body_.npoints());   // ndesign() == 1
  body_.jacobian(alpha, jac.data());
}
