/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer

   svg_shape: write a self-contained SVG overlaying the initial and the
   optimized symmetric-Bezier body outlines. No external dependencies --
   opens in any browser. Geometry comes from BezierGeom.

   Ported near-verbatim from the reference ad_src's svg_shape.h/.cpp
   (drag_init/drag_final params renamed to value_init/value_final, plus a
   new `objective_name` param for the header text, since the annotated
   scalar is now either DragObjective's or HeatFluxObjective's value).
   Stays Bezier-specific -- the API already takes alpha[4] directly; a
   second shape family would need its own renderer, deferred until one
   exists.
------------------------------------------------------------------------- */

#ifndef SVG_SHAPE_H
#define SVG_SHAPE_H

#include <string>

// Draw both bodies into one SVG at `path`:
//   - initial shape: faint dashed grey outline + control points
//   - optimized shape: bold colored fill + outline + control points
// alpha_init / alpha_final are length-4 design vectors [x1,y1,x2,y2].
// objective_name (e.g. "drag" or "heatflux") labels the header/legend.
// Returns 0 on success, nonzero if the file could not be written.
int write_shapes_svg(const std::string &path,
                     const double alpha_init[4], double value_init,
                     const double alpha_final[4], double value_final,
                     double chord, int nseg,
                     const std::string &objective_name);

#endif
