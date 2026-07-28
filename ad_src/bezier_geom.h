/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer

   bezier_geom: pure geometry kernel mapping Bezier control points
   (the design vector alpha) to line segments + outward normals.

   Design rules (do not break -- this is the AD seam):
     - free functions, no stateful classes in the core API
     - no SPARTA headers, no Pointers, no MPI, no file I/O
     - inputs and outputs are plain scalar arrays owned by the caller
     - discretization (nseg) is an int and never part of the derivative

   Conventions (match SPARTA, see surf.h / Surf::compute_line_normal):
     - 2D, z = 0
     - segment i runs pts[i] -> pts[i+1]
     - outward normal = normalize( +z cross (p2 - p1) ) = (-dy, dx)/len,
       i.e. the normal points to the LEFT of the travel direction
       => traverse a closed body CLOCKWISE (negative signed area) so
          normals point away from the body, into the gas
          (verified against examples/circle/data.circle ordering)
------------------------------------------------------------------------- */

#ifndef SPARTA_BEZIER_GEOM_H
#define SPARTA_BEZIER_GEOM_H

namespace BezierGeom {

// -------------------------------------------------------------------
// core API (pure functions on arrays; templatable to AD types later)
// -------------------------------------------------------------------

// Evaluate the Bezier curve of degree nctrl-1 at parameter t in [0,1]
// via de Casteljau (numerically stable, AD-friendly: only +,*).
//   alpha : control points, length 2*nctrl, layout [x0,y0,x1,y1,...]
//   xy    : out, length 2, point on curve
void bezier_point(int nctrl, const double *alpha, double t, double *xy);

// Bernstein basis value B_{k,n}(t) with n = nctrl-1.
// The curve is linear in alpha:  x(t) = sum_k B_{k,n}(t) * alpha_x[k],
// so d(point_x)/d(alpha_x[k]) = B_{k,n}(t) and cross terms are zero.
// Used by tests and to assemble the analytic Jacobian.
double bernstein(int n, int k, double t);

// Tessellate the curve into nseg segments at t_i = i/nseg.
//   alpha : control points, length 2*nctrl
//   pts   : out, length 2*(nseg+1), shared segment endpoints
//   norms : out, length 2*nseg, unit outward normal per segment
// Requires nctrl >= 2, nseg >= 1, and no zero-length segments
// (degenerate segments leave the normal undefined; validate with
// min_segment_length()).
void bezier_to_lines(int nctrl, const double *alpha,
                     int nseg, double *pts, double *norms);

// Analytic Jacobian of the tessellation: d(pts)/d(alpha),
// row-major, shape [2*(nseg+1)] x [2*nctrl].
//   jac[r*(2*nctrl) + c] = d pts[r] / d alpha[c]
// Independent of alpha (Bezier is linear in its control points),
// exact and cheap. This is the "glue" factor for chaining with
// solver sensitivities:  d(drag)/d(alpha) = d(drag)/d(pts) * jac.
void bezier_to_lines_jacobian(int nctrl, int nseg, double *jac);

// -------------------------------------------------------------------
// symmetric body: the shape-optimization test case
// -------------------------------------------------------------------

// Closed 2D body, symmetric about the x-axis, built from two mirrored
// cubic Bezier halves in body coordinates (nose at origin):
//   upper: P0=(0,0), P1=alpha1, P2=alpha2, P3=(chord,0)
//   lower: mirror of upper (y -> -y)
// Design vector alpha = [x1,y1,x2,y2] (y1,y2 > 0 for a valid body).
//
//   alpha : length 4
//   chord : body length (fixed, not a design variable)
//   nseg  : segments PER HALF; total segments = 2*nseg
//   pts   : out, length 2*(2*nseg+1); traversal is nose -> tail along
//           the upper half, tail -> nose along the lower half
//           (clockwise => outward normals); last point == first point
//   norms : out, length 2*(2*nseg), unit outward normal per segment
void symmetric_body_to_lines(const double *alpha, double chord,
                             int nseg, double *pts, double *norms);

// Analytic Jacobian d(pts)/d(alpha), row-major [2*(2*nseg+1)] x [4].
// Alpha-independent (linear map). Nonzeros: with n=3, t_k = k/nseg,
//   upper pt k: d x/d x1 = B_{1,3}(t_k), d x/d x2 = B_{2,3}(t_k),
//               d y/d y1 = B_{1,3}(t_k), d y/d y2 = B_{2,3}(t_k)
//   lower pt  : same for x, negated for y
void symmetric_body_jacobian(int nseg, double *jac);

// -------------------------------------------------------------------
// validation helpers (used by the command layer and tests)
// -------------------------------------------------------------------

// Shortest segment length in a tessellation; caller checks it against
// a tolerance before handing geometry to SPARTA.
double min_segment_length(int nseg, const double *pts);

// 1 if the polyline is closed within tol (first point == last point),
// else 0. read_surf-style watertight bodies need a closed loop.
int is_closed(int nseg, const double *pts, double tol);

// Signed area of the polyline treated as a closed loop (shoelace).
// Positive => counter-clockwise traversal => normals point INTO the
// enclosed body (wrong for flow around a body). Negative => clockwise
// => normals point outward into the gas, matching SPARTA bodies.
double signed_area(int nseg, const double *pts);

// -------------------------------------------------------------------
// convenience owner for tests/drivers. NOT part of the AD seam --
// the command layer may use raw arrays directly.
// -------------------------------------------------------------------

struct Tessellation {
  int nseg;
  double *pts;     // length 2*(nseg+1)
  double *norms;   // length 2*nseg

  Tessellation(int nctrl, const double *alpha, int nseg_);
  ~Tessellation();

  // non-copyable (owns raw buffers)
  Tessellation(const Tessellation &) = delete;
  Tessellation &operator=(const Tessellation &) = delete;
};

}  // namespace BezierGeom

#endif
