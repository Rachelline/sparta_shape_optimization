/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   shape_case.cpp -- see shape_case.h
------------------------------------------------------------------------- */

#include "shape_case.h"
#include "library.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <unistd.h>

namespace {

void die(const char *msg)
{
  std::fprintf(stderr, "shape_case: %s\n", msg);
  std::exit(1);
}

void cmd(void *spa, const char *str)
{
  char buf[1024];
  std::snprintf(buf, sizeof(buf), "%s", str);
  sparta_command(spa, buf);
}

// SPARTA surf-data format: npt points, npt lines (closed loop,
// line i connects point i -> point (i+1)%npt). Generic over any
// closed polyline -- not shape-family-specific.
void write_surf_tmp(const char *path, const double *pts, int npt)
{
  FILE *fp = std::fopen(path, "w");
  if (!fp) die("cannot open surf tmp file for writing");
  std::fprintf(fp, "shape_case generated surf\n\n");
  std::fprintf(fp, "%d points\n%d lines\n\nPoints\n\n", npt, npt);
  for (int i = 0; i < npt; i++)
    std::fprintf(fp, "%d %.15g %.15g\n", i + 1, pts[2 * i], pts[2 * i + 1]);
  std::fprintf(fp, "\nLines\n\n");
  for (int i = 0; i < npt; i++)
    std::fprintf(fp, "%d %d %d\n", i + 1, i + 1, (i + 1) % npt + 1);
  std::fclose(fp);
}

#ifdef SPARTA_AD
// Side file read_surf.cpp's SPARTA_AD_SEED_JACFILE hook (src/read_surf.cpp,
// read_points()) consumes: one line per point (same order write_surf_tmp
// used), 2*ndesign doubles per line (dx/dalpha_j dy/dalpha_j, j=0..ndesign-1).
// jac is shape.jacobian()'s own row-major output, [2*(npt+1)] x [ndesign];
// only the first npt points' rows are written, matching write_surf_tmp's
// own "drop the duplicate closing point" convention.
void write_jac_tmp(const char *path, const double *jac, int npt, int ndesign)
{
  FILE *fp = std::fopen(path, "w");
  if (!fp) die("cannot open jacobian tmp file for writing");
  for (int i = 0; i < npt; i++) {
    for (int c = 0; c < ndesign; c++) {
      std::fprintf(fp, "%.15g %.15g%s",
                   jac[(2 * i) * ndesign + c], jac[(2 * i + 1) * ndesign + c],
                   c + 1 < ndesign ? "  " : "\n");
    }
  }
  std::fclose(fp);
}
#endif

}  // namespace

double evaluate(const Parametrization &shape, const Objective &obj,
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
  int npt = shape.nsegments(c.nseg);      // unique loop points == line count
  int npt_closed = npt + 1;               // to_lines()'s own convention

  double *pts = new double[2 * npt_closed];
  double *norms = new double[2 * npt];
  shape.to_lines(alpha, c.chord, c.nseg, pts, norms);

  std::string why;
  if (!shape.validate(c.nseg, pts, &why)) die(why.c_str());

  for (int i = 0; i < npt; i++) {
    pts[2 * i]     += c.origin[0];
    pts[2 * i + 1] += c.origin[1];
    if (pts[2 * i] <= 0.0 || pts[2 * i] >= c.boxhi[0] ||
        pts[2 * i + 1] <= 0.0 || pts[2 * i + 1] >= c.boxhi[1])
      die("body outside box");
  }

  char surfpath[256];
  std::snprintf(surfpath, sizeof(surfpath), "tmp_surf_%d_%d.data",
               (int) getpid(), seed);
  write_surf_tmp(surfpath, pts, npt);
  delete[] pts;
  delete[] norms;

#ifdef SPARTA_AD
  char jacpath[256];
  std::snprintf(jacpath, sizeof(jacpath), "tmp_jac_%d_%d.data",
               (int) getpid(), seed);
  {
    double *jac = new double[2 * npt_closed * ndesign];
    shape.jacobian(c.nseg, jac);
    write_jac_tmp(jacpath, jac, npt, ndesign);
    delete[] jac;
  }
#endif

  void *spa;
  char a0[] = "sparta", a1[] = "-log", a2[] = "none";
  char a3[] = "-screen";
  char *argv_verbose[] = {a0, a1, a2};
  char *argv_quiet[]   = {a0, a1, a2, a3, a2};
  if (c.verbose) sparta_open_no_mpi(3, argv_verbose, &spa);
  else           sparta_open_no_mpi(5, argv_quiet, &spa);

  char line[512];

  std::snprintf(line, sizeof(line), "seed %d", seed);
  cmd(spa, line);
  cmd(spa, "dimension 2");
  cmd(spa, "global gridcut 0.0 comm/sort yes");
  cmd(spa, "boundary o r p");
  std::snprintf(line, sizeof(line), "create_box 0 %.15g 0 %.15g -0.5 0.5",
               c.boxhi[0], c.boxhi[1]);
  cmd(spa, line);
  std::snprintf(line, sizeof(line), "create_grid %d %d 1", c.grid[0], c.grid[1]);
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
  if (c.verbose) cmd(spa, "stats 100");

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

  sparta_close(spa);
  unlink(surfpath);
#ifdef SPARTA_AD
  unlink(jacpath);
#endif

  return value;
}

double evaluate_avg(const Parametrization &shape, const Objective &obj,
                    const double *alpha, const int *seeds, int nseeds,
                    const RunConfig &c, double *grad)
{
  int ndesign = shape.ndesign();
  double *g = grad ? new double[ndesign] : 0;
  if (grad) for (int j = 0; j < ndesign; j++) grad[j] = 0.0;

  double sum = 0.0;
  for (int k = 0; k < nseeds; k++) {
    sum += evaluate(shape, obj, alpha, seeds[k], c, g);
    if (grad) for (int j = 0; j < ndesign; j++) grad[j] += g[j];
  }
  if (grad) for (int j = 0; j < ndesign; j++) grad[j] /= nseeds;
  delete[] g;

  return sum / nseeds;
}

double grad_fd(const Parametrization &shape, const Objective &obj,
              const double *alpha, const int *seeds, int nseeds,
              const RunConfig &c, double h, double *grad)
{
  int ndesign = shape.ndesign();
  double *a = new double[ndesign];
  for (int j = 0; j < ndesign; j++) a[j] = alpha[j];

  for (int j = 0; j < ndesign; j++) {
    double orig = a[j];

    a[j] = orig + h;
    double plus = evaluate_avg(shape, obj, a, seeds, nseeds, c);

    a[j] = orig - h;
    double minus = evaluate_avg(shape, obj, a, seeds, nseeds, c);

    a[j] = orig;
    grad[j] = (plus - minus) / (2.0 * h);
  }

  delete[] a;
  return evaluate_avg(shape, obj, alpha, seeds, nseeds, c);
}
