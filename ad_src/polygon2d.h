/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer

   polygon2d: general-purpose closed-polygon utilities over SurfMesh's
   "npt unique points, wraparound mod npt" convention (point i connects
   to point (i+1)%npt -- no duplicate closing point). Used by both
   BezierShape::validate() and MinSizeConstraint's 2D branch; previously
   duplicated as private helpers inside min_area_constraint.cpp and (in
   a slightly different, duplicate-closing-point convention) inside
   bezier_geom.h.
------------------------------------------------------------------------- */

#ifndef SPARTA_POLYGON2D_H
#define SPARTA_POLYGON2D_H

// Shoelace signed area. Positive => counter-clockwise (normals point
// INTO the body, wrong for flow around a body). Negative => clockwise
// (normals point outward into the gas, matching SPARTA bodies).
double polygon_signed_area(const double *pts, int npt);

// d(signed_area)/d(pts), length 2*npt: d/d(x_k) = 0.5*(y_{k+1}-y_{k-1}),
// d/d(y_k) = 0.5*(x_{k-1}-x_{k+1}), indices wrapping mod npt.
void polygon_area_point_grad(const double *pts, int npt, double *dA_dpts);

// Shortest edge length over the closed loop (wrapping mod npt).
double polygon_min_edge_length(const double *pts, int npt);

#endif
