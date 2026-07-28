/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer

   parametrization.h: the Parametrization interface. Maps a design
   vector alpha -> a closed, clockwise polyline (points + outward
   normals) in body coordinates, ready for shape_case.cpp to translate
   into the simulation box and hand to SPARTA's read_surf.

   This is the ONE interface a new shape family implements. Everything
   downstream (shape_case.cpp, and eventually ShapeTNLP) is written
   against this interface and never against a concrete shape -- adding a
   shape family is a new small class, not a change anywhere else.

   Design rule (matches bezier_geom.h's own "AD seam" rule): to_lines()/
   jacobian() implementations should stay pure array-in/array-out -- no
   SPARTA headers, no I/O -- so a future AD engine can seed jacobian()'s
   output directly into Surf::Line points without this interface, or any
   implementation of it, needing to change.
------------------------------------------------------------------------- */

#ifndef SPARTA_PARAMETRIZATION_H
#define SPARTA_PARAMETRIZATION_H

#include <string>

class Parametrization {
 public:
  virtual ~Parametrization() {}

  // size of the design vector alpha
  virtual int ndesign() const = 0;

  // total segment count to_lines()/jacobian() will produce for a given
  // requested discretization `nseg` (shape-family-defined -- e.g. for
  // BezierParametrization, nseg is segments PER HALF and this returns
  // 2*nseg). The caller derives every buffer size from this:
  //   npoints = nsegments(nseg) + 1        (closed loop, last == first)
  //   pts     length 2*npoints doubles
  //   norms   length 2*nsegments(nseg) doubles
  //   jac     length 2*npoints*ndesign() doubles, row-major
  virtual int nsegments(int nseg) const = 0;

  // alpha -> closed clockwise polyline, in BODY coordinates (not yet
  // translated into the simulation box -- shape_case.cpp does that).
  //   alpha : design vector, length ndesign()
  //   chord : a fixed (non-design) length scale; shape families that
  //           don't need one may ignore it
  //   nseg  : requested discretization; total segment/point counts and
  //           required buffer sizes are shape-family-defined (see each
  //           implementation's header for its exact convention)
  //   pts   : out, caller-allocated, closed loop (last point == first)
  //   norms : out, caller-allocated, one unit outward normal per segment
  virtual void to_lines(const double *alpha, double chord, int nseg,
                        double *pts, double *norms) const = 0;

  // d(pts)/d(alpha), row-major: jac[r*ndesign() + c] = d(pts[r])/d(alpha[c]).
  // Exact at whatever alpha to_lines() was (or would be) called with;
  // the interface doesn't assume alpha-independence, only that shape
  // families for which the Jacobian genuinely doesn't depend on alpha
  // (e.g. any Bezier-linear map) are free to ignore that and precompute.
  virtual void jacobian(int nseg, double *jac) const = 0;

  // Shape-family-specific validity checks (self-intersection rules,
  // parameter ordering, anything beyond generic polyline sanity --
  // closed/clockwise/min-segment-length -- which shape_case.cpp checks
  // itself via BezierGeom::signed_area/min_segment_length on every
  // shape's output, since those are general-purpose polyline utilities
  // that happen to live in bezier_geom.h today). Default: no extra
  // checks.
  virtual bool validate(int nseg, const double *pts, std::string *why) const {
    (void) nseg; (void) pts; (void) why;
    return true;
  }

  // Box constraints an optimizer should apply to alpha, length
  // ndesign() each. Not consumed anywhere yet (ShapeTNLP is a later
  // milestone) -- part of the interface now so implementations declare
  // sane defaults from day one.
  virtual void bounds(double *lo, double *hi) const = 0;
};

#endif
