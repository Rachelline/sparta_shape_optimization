/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   shape_case.cpp -- see shape_case.h
------------------------------------------------------------------------- */

#include "shape_case.h"
#include "sparta_util.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <unistd.h>

namespace {

// SurfMesh -> SPARTA surf-data file. dim==2 writes a Lines section (each
// element 2 point indices); dim==3 writes a Triangles section (each
// element 3 point indices). Point/element counts and the element kind
// keyword are the only things that differ between the two.
void write_surf_tmp(const char *path, const SurfMesh &m)
{
  FILE *fp = std::fopen(path, "w");
  if (!fp) die("cannot open surf tmp file for writing");

  int npt = m.npoints(), nel = m.nelems();
  const char *elem_kind = (m.dim == 2) ? "lines" : "triangles";
  std::fprintf(fp, "shape_case generated surf\n\n");
  std::fprintf(fp, "%d points\n%d %s\n\nPoints\n\n", npt, nel, elem_kind);
  for (int i = 0; i < npt; i++) {
    std::fprintf(fp, "%d", i + 1);
    for (int d = 0; d < m.dim; d++) std::fprintf(fp, " %.15g", m.pts[m.dim * i + d]);
    std::fprintf(fp, "\n");
  }
  std::fprintf(fp, "\n%s\n\n", (m.dim == 2) ? "Lines" : "Triangles");
  for (int i = 0; i < nel; i++) {
    std::fprintf(fp, "%d", i + 1);
    for (int d = 0; d < m.dim; d++) std::fprintf(fp, " %d", m.elems[m.dim * i + d] + 1);
    std::fprintf(fp, "\n");
  }
  std::fclose(fp);
}

#ifdef SPARTA_AD
// Side file read_surf.cpp's SPARTA_AD_SEED_JACFILE hook (src/read_surf.cpp,
// read_points()) consumes: one line per point, `dim` doubles per design
// column (dx/dalpha_k dy/dalpha_k [dz/dalpha_k], k=0..ndesign-1). jac is
// shape.jacobian()'s own row-major output, [dim*npoints] x [ndesign].
void write_jac_tmp(const char *path, const double *jac, int npt, int dim,
                   int ndesign)
{
  FILE *fp = std::fopen(path, "w");
  if (!fp) die("cannot open jacobian tmp file for writing");
  for (int i = 0; i < npt; i++) {
    for (int c = 0; c < ndesign; c++) {
      for (int d = 0; d < dim; d++)
        std::fprintf(fp, "%.15g%s", jac[(dim * i + d) * ndesign + c],
                    (d + 1 < dim) ? " " : "");
      std::fprintf(fp, "%s", (c + 1 < ndesign) ? "  " : "\n");
    }
  }
  std::fclose(fp);
}
#endif

}  // namespace

double evaluate(const Shape &shape, const Objective &obj,
                const double *alpha, int seed, const RunConfig &c,
                double *grad)
{
  if (access(c.species_file, R_OK) != 0)
    die("cannot read species_file (SPARTA fails silently with -screen "
        "none if this is wrong)");
  if (c.collisions && access(c.vss_file, R_OK) != 0)
    die("cannot read vss_file (SPARTA fails silently with -screen none "
        "if this is wrong)");

  int ndesign = shape.ndesign();
  int dim = shape.dim();

  SurfMesh m;
  shape.to_mesh(alpha, m);

  std::string why;
  if (!shape.validate(m, &why)) die(why);

  int npt = m.npoints();
  for (int i = 0; i < npt; i++) {
    for (int d = 0; d < dim; d++) {
      m.pts[dim * i + d] += c.origin[d];
      if (m.pts[dim * i + d] <= c.boxlo[d] || m.pts[dim * i + d] >= c.boxhi[d])
        die("body outside box");
    }
  }

  char surfpath[256];
  std::snprintf(surfpath, sizeof(surfpath), "tmp_surf_%d_%d.data",
               (int) getpid(), seed);
  write_surf_tmp(surfpath, m);

#ifdef SPARTA_AD
  char jacpath[256];
  std::snprintf(jacpath, sizeof(jacpath), "tmp_jac_%d_%d.data",
               (int) getpid(), seed);
  {
    std::vector<double> jac;
    shape.jacobian(alpha, jac);
    write_jac_tmp(jacpath, jac.data(), npt, dim, ndesign);
  }
#endif

  void *spa = open_sparta(c.verbose);

  char line[512];

  std::snprintf(line, sizeof(line), "seed %d", seed);
  cmd(spa, line);
  std::snprintf(line, sizeof(line), "dimension %d", dim);
  cmd(spa, line);
  cmd(spa, "global gridcut 0.0 comm/sort yes");
  std::snprintf(line, sizeof(line), "boundary %s", c.boundary);
  cmd(spa, line);
  std::snprintf(line, sizeof(line), "create_box %.15g %.15g %.15g %.15g %.15g %.15g",
               c.boxlo[0], c.boxhi[0], c.boxlo[1], c.boxhi[1], c.boxlo[2], c.boxhi[2]);
  cmd(spa, line);
  std::snprintf(line, sizeof(line), "create_grid %d %d %d",
               c.grid[0], c.grid[1], c.grid[2]);
  cmd(spa, line);
  cmd(spa, "balance_grid rcb cell");
  std::snprintf(line, sizeof(line), "global nrho %.15g fnum %.15g", c.nrho, c.fnum);
  cmd(spa, line);
  std::snprintf(line, sizeof(line), "species %s %s", c.species_file, c.species_names);
  cmd(spa, line);
  std::snprintf(line, sizeof(line), "mixture gas %s vstream %.15g 0 0 temp %.15g",
               c.species_names, c.vstream, c.tinf);
  cmd(spa, line);
#ifdef SPARTA_AD
  // Consumed by src/read_surf.cpp's SPARTA_AD_SEED_JACFILE hook, inside
  // the read_surf command issued right below -- must be set before that
  // command runs, not after.
  setenv("SPARTA_AD_SEED_JACFILE", jacpath, 1);

  // Consumed by src/compute_surf.cpp's ComputeSurf::init() (read once per
  // `compute surf` command). Explicitly set/unset (not just "set when
  // true") so a process that runs evaluate() with different RunConfigs
  // never sees a stale value from an earlier call.
  if (c.score_correction) setenv("SPARTA_AD_SCORE_CORRECTION", "1", 1);
  else unsetenv("SPARTA_AD_SCORE_CORRECTION");
#endif

  std::snprintf(line, sizeof(line), "read_surf %s", surfpath);
  cmd(spa, line);

  if (c.specular) cmd(spa, "surf_collide 1 specular");
  else {
    std::snprintf(line, sizeof(line), "surf_collide 1 diffuse %.15g %.15g",
                 c.wall_temp, c.wall_accom);
    cmd(spa, line);
  }
  cmd(spa, "surf_modify all collide 1");

  if (c.collisions) {
    std::snprintf(line, sizeof(line), "collide vss gas %s", c.vss_file);
    cmd(spa, line);
  }

  cmd(spa, "fix in emit/face gas xlo twopass");
  std::snprintf(line, sizeof(line), "timestep %.15g", c.tstep);
  cmd(spa, line);
  if (c.verbose) {
    std::snprintf(line, sizeof(line), "stats %d", c.stats_every);
    cmd(spa, line);
  }

  std::snprintf(line, sizeof(line), "run %d", c.nsettle);
  cmd(spa, line);

  // fix ave/surf's Nfreq window must land on an ABSOLUTE timestep that's
  // a multiple of navg (Nevery=1, Nfreq=Nrepeat=navg below) -- not just
  // navg steps after whenever the fix happened to be created. Without
  // this reset, the window only lands correctly when nsettle is itself
  // a multiple of navg (true by coincidence whenever nsettle==navg, the
  // RunConfig default -- but a silent footgun for any other nsettle/navg
  // combination, where it fails with "Fix used in compute reduce not
  // computed at compatible time"). Resetting to 0 right before the
  // averaging fix is created decouples nsettle and navg completely.
  cmd(spa, "reset_timestep 0");

  obj.setup(spa, c.navg);
  std::snprintf(line, sizeof(line), "run %d", c.navg);
  cmd(spa, line);

  double value = obj.extract(spa, ndesign, grad);

  close_sparta(spa);
  unlink(surfpath);
#ifdef SPARTA_AD
  unlink(jacpath);
#endif

  return value;
}

double evaluate_avg(const Shape &shape, const Objective &obj,
                    const double *alpha, const int *seeds, int nseeds,
                    const RunConfig &c, double *grad)
{
  int ndesign = shape.ndesign();
  std::vector<double> g;
  if (grad) {
    g.assign(ndesign, 0.0);
    for (int j = 0; j < ndesign; j++) grad[j] = 0.0;
  }

  double sum = 0.0;
  for (int k = 0; k < nseeds; k++) {
    sum += evaluate(shape, obj, alpha, seeds[k], c, grad ? g.data() : nullptr);
    if (grad) for (int j = 0; j < ndesign; j++) grad[j] += g[j];
  }
  if (grad) for (int j = 0; j < ndesign; j++) grad[j] /= nseeds;

  return sum / nseeds;
}

double grad_fd(const Shape &shape, const Objective &obj,
              const double *alpha, const int *seeds, int nseeds,
              const RunConfig &c, double h, double *grad)
{
  int ndesign = shape.ndesign();
  std::vector<double> a(alpha, alpha + ndesign);

  for (int j = 0; j < ndesign; j++) {
    double orig = a[j];

    a[j] = orig + h;
    double plus = evaluate_avg(shape, obj, a.data(), seeds, nseeds, c);

    a[j] = orig - h;
    double minus = evaluate_avg(shape, obj, a.data(), seeds, nseeds, c);

    a[j] = orig;
    grad[j] = (plus - minus) / (2.0 * h);
  }

  return evaluate_avg(shape, obj, alpha, seeds, nseeds, c);
}
