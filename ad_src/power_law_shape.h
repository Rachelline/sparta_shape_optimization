/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer

   power_law_shape.h: Shape adapter wrapping PowerLawBody (a 3D
   axisymmetric body of revolution, see power_law_body.h) -- the 3D
   counterpart of bezier_shape.h. Previously PowerLawBody could not
   implement Parametrization at all (its to_lines()-shaped 2D interface
   had no way to express a triangulated 3D mesh); Shape's explicit
   SurfMesh::elems connectivity removes that restriction, so the 3D body
   can now go through the same shared runner (shape_case.cpp), ShapeTNLP,
   and Constraint machinery as any 2D shape.

   ndesign() == 1 (the power-law exponent n). PowerLawBody's own
   jacobian(alpha, jac) already took alpha (dr/dn = r*ln(x/L) genuinely
   depends on n) -- it just had nowhere to plug in until this interface
   existed to carry alpha through jacobian() too.
------------------------------------------------------------------------- */

#ifndef SPARTA_POWER_LAW_SHAPE_H
#define SPARTA_POWER_LAW_SHAPE_H

#include "shape.h"
#include "power_law_body.h"

class PowerLawShape : public Shape {
 public:
  PowerLawShape(double L = 1.0, double Rmax = 0.3, int nx = 8, int ntheta = 12)
    : body_(L, Rmax, nx, ntheta) {}

  int dim() const override { return 3; }
  int ndesign() const override { return 1; }
  void bounds(double *lo, double *hi) const override { body_.bounds(lo, hi); }

  void to_mesh(const double *alpha, SurfMesh &m) const override;
  void jacobian(const double *alpha, std::vector<double> &jac) const override;

 private:
  PowerLawBody body_;
};

#endif
