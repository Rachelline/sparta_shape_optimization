/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer

   power_law_body: a 3D axisymmetric body of revolution,
       r(x) = Rmax * (x/L)^n,   x in [0, L]
   parametrized by a SINGLE design variable, the power-law exponent n.
   n < 1 gives a blunt, rapidly-flaring nose (r rises steeply near the
   apex); n = 1 is a cone; n > 1 gives a slender, pointy spike (r stays
   small for most of the length, then flares near the base). This is the
   classic power-law/spike forebody shape from hypersonic aerothermal
   trade studies -- exactly the "pointy vs blunt" knob heat-flux shape
   optimization is meant to explore.

   Deliberately NOT an implementation of parametrization.h: that
   interface's to_lines()/jacobian() are 2D-line-shaped (a closed
   polyline + per-segment normals for a `dimension 2` SPARTA surf file).
   A body of revolution needs a triangulated 3D mesh instead, so this is
   a parallel, minimal interface used by power_law_case.cpp/
   power_law_main.cpp only -- not plugged into ShapeTNLP/Objective, which
   are written against the 2D interface.

   Mesh: 1 apex point (x=0, r=0) + NX rings of NTHETA points each
   (x=L/NX .. L) + 1 base-center point (x=L, r=0), triangulated as an
   apex fan, NX-1 rings of side quads (2 tris each), and a base fan.
   Total points = NX*NTHETA + 2. Total triangles = 2*NX*NTHETA.

   Jacobian: exact/analytic. r(x,n) = Rmax*(x/L)^n, so
     dr/dn = r(x,n) * ln(x/L)          (x > 0; apex and base-center are
                                        fixed points, independent of n)
   d(point)/dn = (0, dr/dn*cos(theta), dr/dn*sin(theta)).
------------------------------------------------------------------------- */

#ifndef SPARTA_POWER_LAW_BODY_H
#define SPARTA_POWER_LAW_BODY_H

#include <vector>

class PowerLawBody {
 public:
  PowerLawBody(double L, double Rmax, int nx, int ntheta)
    : L_(L), Rmax_(Rmax), nx_(nx), ntheta_(ntheta) {}

  int ndesign() const { return 1; }
  int npoints() const { return nx_ * ntheta_ + 2; }
  int ntris()   const { return 2 * nx_ * ntheta_; }

  // alpha = {n}. pts: length 3*npoints() (x,y,z triples, body-local
  // coordinates: axis along +x from apex at x=0 to base at x=L).
  // tris: length 3*ntris() point INDICES (0-based), CCW winding when
  // viewed from outside the body (outward normal via (p1-p0)x(p2-p0)).
  void to_tris(const double *alpha, double *pts, int *tris) const;

  // jac: length 3*npoints()*ndesign(), row-major:
  // jac[3*i*ndesign()+3*k+c] = d(pts[3*i+c])/d(alpha[k]).  ndesign()==1
  // so this is just d(pts)/dn, length 3*npoints().
  void jacobian(const double *alpha, double *jac) const;

  void bounds(double *lo, double *hi) const { lo[0] = 0.3; hi[0] = 3.0; }

 private:
  double L_, Rmax_;
  int nx_, ntheta_;
};

#endif
