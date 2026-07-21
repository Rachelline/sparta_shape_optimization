/* AD-CONVERTED: double->sfloat by tools/ad_convert.py (see sfloat.h) */
/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.github.io
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#ifndef SPARTA_GEOMETRY_H
#define SPARTA_GEOMETRY_H

#include "sfloat.h"

namespace Geometry {
  int line_quad_intersect(sfloat *, sfloat *, sfloat *,
                          sfloat *, sfloat *);
  int quad_line_intersect_point(sfloat *, sfloat *, sfloat *,
                                sfloat *, sfloat *, sfloat *);
  int line_touch_quad_face(sfloat *, sfloat *, int, sfloat *, sfloat *);

  int tri_hex_intersect(sfloat *, sfloat *, sfloat *, sfloat *,
                        sfloat *, sfloat *);
  int hex_tri_intersect_point(sfloat *, sfloat *, sfloat *, sfloat *,
                              sfloat *, sfloat *, sfloat *);
  int tri_touch_hex_face(sfloat *, sfloat *, sfloat *, int, sfloat *, sfloat *);
  int tri_on_hex_face(sfloat *, sfloat *, sfloat *, sfloat *, sfloat *);
  int edge_on_hex_face(sfloat *, sfloat *, sfloat *, sfloat *);

  bool line_line_intersect(sfloat *, sfloat *,
                           sfloat *, sfloat *, sfloat *,
                           sfloat *, sfloat &param, int &, int=0);

  bool axi_line_intersect(sfloat, sfloat *, sfloat *, int, sfloat *, sfloat *,
                          sfloat *, sfloat *, sfloat *, int,
                          sfloat *, sfloat *, sfloat &, int &);
  bool axi_horizontal_line(sfloat, sfloat *, sfloat *, sfloat,
                           int &, sfloat &, sfloat &);

  bool line_tri_intersect(sfloat *, sfloat *,
                          sfloat *, sfloat *, sfloat *, sfloat *,
                          sfloat *, sfloat &param, int &);
  bool line_tri_intersect_noeps(sfloat *, sfloat *,
                                sfloat *, sfloat *, sfloat *, sfloat *,
                                sfloat *, sfloat &param, int &);
  int whichside(sfloat *, sfloat *, sfloat, sfloat, sfloat);
  int point_on_hex(sfloat *, sfloat *, sfloat *);
  int point_in_hex(sfloat *, sfloat *, sfloat *);
  int point_in_tri(sfloat *, sfloat *, sfloat *, sfloat *, sfloat *);

  sfloat distsq_point_line(sfloat *, sfloat *, sfloat *);
  sfloat distsq_point_tri(sfloat *, sfloat *, sfloat *, sfloat *, sfloat *);

  sfloat dist_line_quad(sfloat *, sfloat *, sfloat *, sfloat *);
  sfloat dist_tri_hex(sfloat *, sfloat *, sfloat *, sfloat *,
                      sfloat *, sfloat *);

  sfloat line_fraction(sfloat *, sfloat *, sfloat *);
  sfloat tri_fraction(sfloat *, sfloat *, sfloat *, sfloat *);
  sfloat poly_area(int, sfloat *, sfloat *);
}

#endif
