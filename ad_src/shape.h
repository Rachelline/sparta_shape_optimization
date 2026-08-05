/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer

   shape.h: the Shape interface. Maps a design vector alpha -> a closed
   surface mesh (2D polyline or 3D triangulation) in body coordinates,
   ready for shape_case.cpp to translate into the simulation box and
   hand to SPARTA's read_surf. Replaces parametrization.h -- unifies
   what used to be two parallel, unconnected interfaces (Parametrization
   for 2D Bezier bodies wired into ShapeTNLP/Constraint, and PowerLawBody
   as a standalone 3D-only interface used by nothing else).

   Four differences from the old Parametrization interface:
     - jacobian() now takes alpha. The old signature omitted it, silently
       assuming alpha-independence -- true for Bezier (a linear map), but
       false for a power-law body (dr/dn = r*ln(x/L) genuinely depends on
       n). That was a latent correctness trap for any alpha-dependent
       shape family; this interface doesn't let a new implementation
       make that mistake.
     - Explicit connectivity (SurfMesh::elems) replaces the old "closed
       loop, last point == first, the writer silently drops the
       duplicate" convention -- and generalizes past loops to triangle
       meshes.
     - chord/nseg (or nx/ntheta) move into each Shape's own constructor,
       since they're per-shape-family discretization choices, not
       something a generic caller should thread through every call.
     - The unused per-segment `norms` output is dropped: SPARTA derives
       surface normals itself from element winding (verified: nothing
       downstream ever read the old norms buffer).

   Design rule (matches bezier_geom.h's own "AD seam" rule): to_mesh()/
   jacobian() implementations should stay pure array-in/array-out -- no
   SPARTA headers, no I/O -- so a future AD engine can seed jacobian()'s
   output directly without this interface, or any implementation of it,
   needing to change.
------------------------------------------------------------------------- */

#ifndef SPARTA_SHAPE_H
#define SPARTA_SHAPE_H

#include <string>
#include <vector>

// dim*npoints() coordinates + connectivity, in BODY coordinates (not yet
// translated into the simulation box -- shape_case.cpp does that).
// elems holds `dim` point indices per element (2 = line segment
// endpoints, 3 = triangle vertices), 0-based into pts.
struct SurfMesh {
  int dim = 2;
  std::vector<double> pts;
  std::vector<int> elems;

  int npoints() const { return dim > 0 ? (int) pts.size() / dim : 0; }
  int nelems()  const { return dim > 0 ? (int) elems.size() / dim : 0; }
};

class Shape {
 public:
  virtual ~Shape() {}

  virtual int dim() const = 0;        // 2 or 3
  virtual int ndesign() const = 0;    // size of the design vector alpha

  // Box constraints an optimizer should apply to alpha, length
  // ndesign() each. Consumed by ShapeTNLP::get_bounds_info().
  virtual void bounds(double *lo, double *hi) const = 0;

  // alpha -> closed surface mesh, in body coordinates.
  virtual void to_mesh(const double *alpha, SurfMesh &m) const = 0;

  // d(mesh.pts)/d(alpha), row-major: jac[(dim*i+d)*ndesign()+k] =
  // d(pts[dim*i+d])/d(alpha[k]). Length dim*npoints()*ndesign(). Exact
  // at the alpha to_mesh() was (or would be) called with; shape families
  // whose Jacobian genuinely doesn't depend on alpha (any linear map,
  // e.g. Bezier) are free to ignore the argument and precompute.
  virtual void jacobian(const double *alpha, std::vector<double> &jac) const = 0;

  // Shape-family-specific validity checks (self-intersection rules,
  // parameter ordering, anything beyond what shape_case.cpp already
  // checks generically -- the simulation-box bounds). Default: none.
  virtual bool validate(const SurfMesh &m, std::string *why) const {
    (void) m; (void) why;
    return true;
  }
};

#endif
