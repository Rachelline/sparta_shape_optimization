/* AD-CONVERTED: double->sfloat by tools/ad_convert.py (see sfloat.h) */
/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.github.io
   Steve Plimpton, sjplimp@gmail.com
   Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#ifdef COMMAND_CLASS

CommandStyle(create_isurf,CreateISurf)

#else

#ifndef SPARTA_CREATE_ISURF_H
#define SPARTA_CREATE_ISURF_H

#include "stdio.h"
#include "pointers.h"
#include "surf.h"

namespace SPARTA_NS {

class CreateISurf : protected Pointers {
 public:
  CreateISurf(class SPARTA *);
  virtual ~CreateISurf();
  virtual void command(int, char **);

 protected:
  int me,nprocs;
  int dim;
  int nglocal;

  // for generating implicit surfaces

  int ggroup;               // group id for grid cells
  int groupbit;
  int ncorner;              // number of corners
  int nmulti;               // number of adjacent neighbors
  int nedge;                // number of cell edges
  sfloat thresh;            // lower threshold for corner values
  sfloat corner[3];         // corners of grid group
  sfloat xyzsize[3];        // size of lowest level cell (must be uniform grid)
  int nxyz[3], Nxyz;        // dimensions of grid
  sfloat **cvalues;         // array of corner point values
  sfloat ***mulvalues;      // array of multi values
  sfloat **tmp_cvalues;     // temporary array of corner point values
  sfloat ***tmp_mulvalues;  // temporary array of multi values
  sfloat **mvalues;         // minimum intersection value
  int **svalues;            // marks corners as in or out
  sfloat ***ivalues;        // point of intersection between corner points

  // buffer between corner point and intersection

  sfloat surfbuffer;

  sfloat **icvalues;        // corner values for Fix Ablate
  int *tvalues;             // vector of per grid cell surf types

  Surf::Line *llines;       // local copy of Surf lines
  Surf::Tri *ltris;         // local copy of Surf tris
  sfloat **cuvalues;        // local copy of custom per-surf data

  int ctype;                // flag for how corners in unknown cells are set
  sfloat mind;              // minimum cell length
  sfloat cin, cout;         // in and out corner values
  sfloat cbufmin, cbufmax;  // corner value buffer
  class FixAblate *ablate;  // ablate fix

  // for communicating

  int **ixyz;             // ix,iy,iz indices (1 to Nxyz) of my cells
                          // in 2d/3d ablate grid (iz = 1 for 2d)

  // various arrays to pass to other processors

  int **sghost;
  sfloat **cghost;
  sfloat ***inghost;
  sfloat ***ighost;       // ditto for my ghost cells communicated to me
  int maxgrid;            // max size of per-cell vectors/arrays
  int maxghost;           // max size of cdelta_ghost

  int *proclist;
  cellint *locallist;
  int *numsend;
  int maxsend;

  sfloat *sbuf;
  int maxsbuf;

  union ubuf {
    double d;
    int64_t i;
    uint64_t u;
    ubuf(double arg) : d(arg) {}
    ubuf(int64_t arg) : i(arg) {}
    ubuf(int arg) : i(arg) {}
    ubuf(uint64_t arg) : u(arg) {}
    ubuf(uint32_t arg) : u(arg) {}
  };

  void process_args(int, char **);

  // functions to set corner/multi values

  void set_corners();
  void set_multi();

  // send/recv values between neighborind cells (similar to fix_ablate)

  void sync(int);
  void sync_voxels();
  void comm_neigh_corners(int);
  void grow_send();
  int walk_to_neigh(int, int, int, int);

  // sets corner values whose cells have surfaces

  void surface_edge2d();
  void surface_edge3d();

  // marks corners which have no surfaces

  void set_inout();

  // find remaining corners

  int find_side_2d();
  int find_side_3d();
  void set_cvalues();
  void set_cvalues_inout();
  void set_cvalues_voxel();
  void set_cvalues_ave();
  void set_cvalues_multi();

  // detects intersection between surfaces and cell edges

  int corner_hit2d(sfloat*, sfloat*, Surf::Line*, sfloat&, int&);
  int corner_hit3d(sfloat*, sfloat*, Surf::Tri*, sfloat&, int&);

  // remove old surfaces

  void remove_old();

  // misc functions

  sfloat param2cval(sfloat, sfloat);
  sfloat interpolate(sfloat, sfloat);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Could not find surf_modify surf-ID

Self-explanatory.

E: Could not find surf_modify sc-ID

Self-explanatory.

E: %d surface elements not assigned to a collision model

All surface elements must be assigned to a surface collision model via
the surf_modify command before a simulation is perforemd.

E: Reuse of surf_collide ID

A surface collision model ID cannot be used more than once.

E: Invalid surf_collide style

Self-explanatory.

*/
