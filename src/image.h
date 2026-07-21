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

#ifndef SPARTA_IMAGE_H
#define SPARTA_IMAGE_H

#include "math.h"
#include "stdio.h"
#include "pointers.h"

namespace SPARTA_NS {

class Image : protected Pointers {
 public:
  int width,height;             // size of image
  sfloat theta,phi;             // view image from theta,phi
  sfloat xctr,yctr,zctr;        // center of image in user coords
  sfloat up[3];                 // up direction in image
  sfloat zoom;                  // zoom factor
  sfloat persp;                 // perspective factor
  sfloat shiny;                 // shininess of objects
  int ssao;                     // SSAO on or off
  int seed;                     // RN seed for SSAO
  sfloat ssaoint;               // strength of shading from 0 to 1
  sfloat *boxcolor;             // color to draw box outline with
  sfloat *gridcolor;            // color to draw grid lines with
  int background[3];            // RGB values of background

  Image(class SPARTA *, int);
  ~Image();
  void buffers();
  void clear();
  void merge();
  void write_JPG(FILE *);
  void write_PNG(FILE *);
  void write_PPM(FILE *);
  void view_params(sfloat, sfloat, sfloat, sfloat, sfloat, sfloat);

  void draw_line(sfloat *, sfloat *, sfloat *, sfloat);
  void draw_box(sfloat (*)[3], sfloat *, sfloat);
  void draw_box2d(sfloat (*)[3], sfloat *, sfloat);
  void draw_axes(sfloat (*)[3], sfloat);
  void draw_sphere(sfloat *, sfloat *, sfloat);
  void draw_brick(sfloat *, sfloat *, sfloat *);
  void draw_cylinder(sfloat *, sfloat *, sfloat *, sfloat, int);
  void draw_triangle(sfloat *, sfloat *, sfloat *, sfloat *);

  int map_dynamic(int);
  int map_reset(int, int, char **);
  int map_minmax(int, sfloat, sfloat);
  sfloat *map_value2color(int, sfloat);

  int addcolor(char *, sfloat, sfloat, sfloat);
  sfloat *element2color(char *);
  sfloat element2diam(char *);
  sfloat *color2rgb(const char *, int index=0);
  int default_colors();

 private:
  int me,nprocs;
  int npixels;

  class ColorMap **maps;
  int nmap;

  sfloat *depthBuffer,*surfaceBuffer;
  sfloat *depthcopy,*surfacecopy;
  char *imageBuffer,*rgbcopy,*writeBuffer;

  // MPI_Gatherv

  int *recvcounts,*displs;

  // constant view params

  sfloat FOV;
  sfloat ambientColor[3];

  sfloat keyLightTheta;
  sfloat keyLightPhi;
  sfloat keyLightColor[3];

  sfloat fillLightTheta;
  sfloat fillLightPhi;
  sfloat fillLightColor[3];

  sfloat backLightTheta;
  sfloat backLightPhi;
  sfloat backLightColor[3];

  sfloat specularHardness;
  sfloat specularIntensity;

  sfloat SSAORadius;
  int SSAOSamples;
  sfloat SSAOJitter;

  // dynamic view params

  sfloat zdist;
  sfloat tanPerPixel;
  sfloat camDir[3],camUp[3],camRight[4],camPos[3];
  sfloat keyLightDir[3],fillLightDir[3],backLightDir[3];
  sfloat keyHalfDir[3];

  // color values

  int ncolors;
  char **username;
  sfloat **userrgb;

  // color maps

  struct MapEntry {
    int single,lo,hi;              // NUMERIC or MINVALUE or MAXVALUE
    sfloat svalue,lvalue,hvalue;   // actual value
    sfloat *color;                 // RGB values
  };

  MapEntry *mentry;
  int nentry;
  sfloat interpolate[3];

  // SSAO RNG

  class RanKnuth *random;

  // internal methods

  void draw_pixel(int, int, sfloat, sfloat *, sfloat*);
  void compute_SSAO();

  // inline functions

  inline sfloat saturate(sfloat v) {
    if (v < 0.0) return 0.0;
    else if (v > 1.0) return 1.0;
    else return v;
  }

  inline sfloat distance(sfloat* a, sfloat* b) {
    return sqrt((a[0] - b[0]) * (a[0] - b[0]) +
                (a[1] - b[1]) * (a[1] - b[1]) +
                (a[2] - b[2]) * (a[2] - b[2]));
  }
};

// ColorMap class

class ColorMap : protected Pointers {
 public:
  int dynamic;                     // 0/1 if lo/hi bounds are static/dynamic

  ColorMap(class SPARTA *, class Image*);
  ~ColorMap();
  int reset(int, char **);
  int minmax(sfloat, sfloat);
  sfloat *value2color(sfloat);

 private:
  class Image *image;              // caller with color2rgb() method
  int mstyle,mrange;               // 2-letter style/range of color map
  int mlo,mhi;                     // bounds = NUMERIC or MINVALUE or MAXVALUE
  sfloat mlovalue,mhivalue;        // user bounds if NUMERIC
  sfloat locurrent,hicurrent;      // current bounds for this snapshot
  sfloat mbinsize,mbinsizeinv;     // bin size for sequential color map
  sfloat interpolate[3];           // local storage for returned RGB color

  struct MapEntry {
    int single,lo,hi;              // NUMERIC or MINVALUE or MAXVALUE
    sfloat svalue,lvalue,hvalue;   // actual value
    sfloat *color;                 // RGB values
  };

  MapEntry *mentry;
  int nentry;
};

}

#endif

/* ERROR/WARNING messages:

E: Invalid image up vector

Up vector cannot be (0,0,0).

*/
