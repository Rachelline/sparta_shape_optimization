/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer

   bezier_geom implementation. See bezier_geom.h for the design rules:
   pure array-to-array functions, no SPARTA dependencies, no I/O.

   AD note: every arithmetic path below uses only +,-,*,/ and sqrt on
   the scalar type, so converting to template <typename T> later is
   mechanical (move definitions header-side or add explicit
   instantiations for double and the AD type).
------------------------------------------------------------------------- */

#include "bezier_geom.h"

#include <cmath>

namespace BezierGeom {

/* ----------------------------------------------------------------------
   de Casteljau evaluation of the curve at parameter t
   work arrays sized on the stack for small nctrl, heap otherwise
------------------------------------------------------------------------- */

void bezier_point(int nctrl, const double *alpha, double t, double *xy)
{
  if (nctrl < 1) {
    xy[0] = xy[1] = 0.0;
    return;
  }

  const int NSTACK = 16;
  double bufx[NSTACK],bufy[NSTACK];
  double *wx = bufx, *wy = bufy;
  if (nctrl > NSTACK) {
    wx = new double[nctrl];
    wy = new double[nctrl];
  }

  for (int i = 0; i < nctrl; i++) {
    wx[i] = alpha[2*i];
    wy[i] = alpha[2*i+1];
  }

  for (int level = 1; level < nctrl; level++)
    for (int i = 0; i < nctrl-level; i++) {
      wx[i] = (1.0-t)*wx[i] + t*wx[i+1];
      wy[i] = (1.0-t)*wy[i] + t*wy[i+1];
    }

  xy[0] = wx[0];
  xy[1] = wy[0];

  if (nctrl > NSTACK) {
    delete [] wx;
    delete [] wy;
  }
}

/* ----------------------------------------------------------------------
   Bernstein basis B_{k,n}(t), computed by the same triangular
   recurrence as de Casteljau (stable for all t in [0,1]):
   seed unit vector e_k, run the recurrence, read off the result
------------------------------------------------------------------------- */

double bernstein(int n, int k, double t)
{
  if (k < 0 || k > n) return 0.0;

  const int NSTACK = 16;
  double buf[NSTACK];
  double *w = buf;
  if (n+1 > NSTACK) w = new double[n+1];

  for (int i = 0; i <= n; i++) w[i] = 0.0;
  w[k] = 1.0;

  for (int level = 1; level <= n; level++)
    for (int i = 0; i <= n-level; i++)
      w[i] = (1.0-t)*w[i] + t*w[i+1];

  double result = w[0];
  if (n+1 > NSTACK) delete [] w;
  return result;
}

/* ----------------------------------------------------------------------
   tessellate curve into nseg segments, endpoints at t_i = i/nseg
   normal convention matches Surf::compute_line_normal:
     norm = normalize( +z cross (p2-p1) ) = (-dy,dx)/len
------------------------------------------------------------------------- */

void bezier_to_lines(int nctrl, const double *alpha,
                     int nseg, double *pts, double *norms)
{
  // endpoints, shared between adjacent segments

  for (int i = 0; i <= nseg; i++) {
    double t = (double) i / (double) nseg;
    bezier_point(nctrl,alpha,t,&pts[2*i]);
  }

  // per-segment unit outward normals

  for (int i = 0; i < nseg; i++) {
    double dx = pts[2*(i+1)]   - pts[2*i];
    double dy = pts[2*(i+1)+1] - pts[2*i+1];
    double len = sqrt(dx*dx + dy*dy);
    norms[2*i]   = -dy/len;
    norms[2*i+1] =  dx/len;
  }
}

/* ----------------------------------------------------------------------
   analytic Jacobian d(pts)/d(alpha), row-major [2*(nseg+1)] x [2*nctrl]
   pts_x(t_i) = sum_k B_{k,n}(t_i) alpha_x[k]  (linear in alpha),
   so jac is alpha-independent:
     d pts[2i]   / d alpha[2k]   = B_{k,n}(t_i)
     d pts[2i+1] / d alpha[2k+1] = B_{k,n}(t_i)
   and zero for x/y cross terms
------------------------------------------------------------------------- */

void bezier_to_lines_jacobian(int nctrl, int nseg, double *jac)
{
  int ncol = 2*nctrl;
  int nrow = 2*(nseg+1);
  int n = nctrl-1;

  for (int r = 0; r < nrow*ncol; r++) jac[r] = 0.0;

  for (int i = 0; i <= nseg; i++) {
    double t = (double) i / (double) nseg;
    for (int k = 0; k < nctrl; k++) {
      double b = bernstein(n,k,t);
      jac[(2*i)  *ncol + 2*k]   = b;   // x wrt x-coord of ctrl pt k
      jac[(2*i+1)*ncol + 2*k+1] = b;   // y wrt y-coord of ctrl pt k
    }
  }
}

/* ----------------------------------------------------------------------
   symmetric body: two mirrored cubic halves, closed, clockwise
   traversal: upper nose->tail (t = 0..1), lower tail->nose (mirrored)
------------------------------------------------------------------------- */

void symmetric_body_to_lines(const double *alpha, double chord,
                             int nseg, double *pts, double *norms)
{
  // upper-half control points in body coordinates

  double ctrl[8];
  ctrl[0] = 0.0;      ctrl[1] = 0.0;
  ctrl[2] = alpha[0]; ctrl[3] = alpha[1];
  ctrl[4] = alpha[2]; ctrl[5] = alpha[3];
  ctrl[6] = chord;    ctrl[7] = 0.0;

  double *up = new double[2*(nseg+1)];
  double *un = new double[2*nseg];
  bezier_to_lines(4,ctrl,nseg,up,un);

  // assemble closed loop: 2*nseg segments, 2*nseg+1 points,
  // last point == first point == nose (0,0)

  for (int i = 0; i <= nseg; i++) {           // upper: nose -> tail
    pts[2*i]   = up[2*i];
    pts[2*i+1] = up[2*i+1];
  }
  for (int j = 1; j <= nseg; j++) {           // lower: tail -> nose
    int i = nseg + j;
    int k = nseg - j;
    pts[2*i]   =  up[2*k];
    pts[2*i+1] = -up[2*k+1];
  }

  delete [] up;
  delete [] un;

  // per-segment unit outward normals (left of clockwise travel)

  int ntot = 2*nseg;
  for (int i = 0; i < ntot; i++) {
    double dx = pts[2*(i+1)]   - pts[2*i];
    double dy = pts[2*(i+1)+1] - pts[2*i+1];
    double len = sqrt(dx*dx + dy*dy);
    norms[2*i]   = -dy/len;
    norms[2*i+1] =  dx/len;
  }
}

/* ----------------------------------------------------------------------
   analytic Jacobian d(pts)/d(alpha), row-major [2*(2*nseg+1)] x [4]
   columns: 0=x1 1=y1 2=x2 3=y2
------------------------------------------------------------------------- */

void symmetric_body_jacobian(int nseg, double *jac)
{
  const int ncol = 4;
  int npt = 2*nseg+1;

  for (int r = 0; r < 2*npt*ncol; r++) jac[r] = 0.0;

  for (int k = 0; k <= nseg; k++) {
    double t = (double) k / (double) nseg;
    double b1 = bernstein(3,1,t);
    double b2 = bernstein(3,2,t);

    // upper point k
    jac[(2*k)  *ncol + 0] = b1;   // dx/dx1
    jac[(2*k)  *ncol + 2] = b2;   // dx/dx2
    jac[(2*k+1)*ncol + 1] = b1;   // dy/dy1
    jac[(2*k+1)*ncol + 3] = b2;   // dy/dy2

    // mirrored lower point at loop index nseg + (nseg-k), skip
    // duplicates of tail (k=nseg) and nose (k=0 -> last==first row)
    int i = 2*nseg - k;
    if (i > nseg && i < 2*nseg+1) {
      jac[(2*i)  *ncol + 0] =  b1;
      jac[(2*i)  *ncol + 2] =  b2;
      jac[(2*i+1)*ncol + 1] = -b1;
      jac[(2*i+1)*ncol + 3] = -b2;
    }
  }

  // closing point (loop index 2*nseg) coincides with the nose (0,0),
  // which is fixed -> its rows stay zero

  return;
}

/* ----------------------------------------------------------------------
   validation helpers
------------------------------------------------------------------------- */

double min_segment_length(int nseg, const double *pts)
{
  double minlen = -1.0;
  for (int i = 0; i < nseg; i++) {
    double dx = pts[2*(i+1)]   - pts[2*i];
    double dy = pts[2*(i+1)+1] - pts[2*i+1];
    double len = sqrt(dx*dx + dy*dy);
    if (minlen < 0.0 || len < minlen) minlen = len;
  }
  return minlen;
}

int is_closed(int nseg, const double *pts, double tol)
{
  double dx = pts[2*nseg]   - pts[0];
  double dy = pts[2*nseg+1] - pts[1];
  return (sqrt(dx*dx + dy*dy) <= tol) ? 1 : 0;
}

double signed_area(int nseg, const double *pts)
{
  // shoelace over the loop, implicitly closing last -> first

  double area2 = 0.0;
  for (int i = 0; i < nseg; i++) {
    double x1 = pts[2*i],     y1 = pts[2*i+1];
    double x2 = pts[2*(i+1)], y2 = pts[2*(i+1)+1];
    area2 += x1*y2 - x2*y1;
  }
  area2 += pts[2*nseg]*pts[1] - pts[0]*pts[2*nseg+1];
  return 0.5*area2;
}

/* ----------------------------------------------------------------------
   convenience owner
------------------------------------------------------------------------- */

Tessellation::Tessellation(int nctrl, const double *alpha, int nseg_)
  : nseg(nseg_)
{
  pts   = new double[2*(nseg+1)];
  norms = new double[2*nseg];
  bezier_to_lines(nctrl,alpha,nseg,pts,norms);
}

Tessellation::~Tessellation()
{
  delete [] pts;
  delete [] norms;
}

}  // namespace BezierGeom
