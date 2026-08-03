/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer

   power_law_main.cpp: a 3D heat-flux shape-optimization demo on a
   PowerLawShape (see power_law_shape.h) -- a single design variable n
   controlling nose bluntness (n<1 blunt/flaring, n=1 cone, n>1 pointy
   spike). Diffuse wall (wall_accom=1) so HeatFluxObjective (etot) is
   nondegenerate, same convention as the 2D Bezier heatflux runs.

   Drives the same shared runner as shape_main.cpp/opt_main.cpp
   (shape_case.h's evaluate_avg()/grad_fd()) instead of a hand-rolled
   duplicate -- previously this file carried its own ~100-line copy of
   evaluate() because PowerLawBody had no way to plug into the (2D-only)
   Parametrization interface; Shape's explicit SurfMesh connectivity
   removed that restriction.

   Deliberately still a standalone driver, not built on ShapeTNLP/IPOPT:
   with ndesign()==1 a full NLP solver is unnecessary machinery. Instead
   a plain gradient-descent loop (AD gradient every step, one SPARTA run
   per step -- the actual point of AD over FD here) with a fixed
   learning rate and box-clamped n. See power_law_body.h::bounds().

   Output (output/pl_<date>_<experiment>/), mirroring opt_main.cpp's own
   convention:
     config.txt        -- run parameters
     trajectory.csv     -- iter,n,value
     result.txt          -- summary
     shape_before.txt    -- initial geometry (points + triangle indices)
     shape_after.txt     -- final geometry, same format
   Rendering shape_before/after.txt into an image is a separate,
   SPARTA-free post-process (tools/ad_verify-style), not done here.
------------------------------------------------------------------------- */

#include "power_law_shape.h"
#include "heat_flux_objective.h"
#include "run_output.h"
#include "shape_case.h"
#include "sparta_util.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <fstream>
#include <string>
#include <unistd.h>
#include <vector>

namespace {

void write_geom_txt(const char *path, const double *pts, int npt,
                    const int *tris, int ntri, double n)
{
  std::ofstream f(path);
  f << "# power-law body, n=" << n << "\n";
  f << "# npoints=" << npt << " ntriangles=" << ntri << "\n";
  f << "POINTS\n";
  for (int i = 0; i < npt; i++)
    f << pts[3*i] << " " << pts[3*i+1] << " " << pts[3*i+2] << "\n";
  f << "TRIANGLES\n";
  for (int i = 0; i < ntri; i++)
    f << tris[3*i] << " " << tris[3*i+1] << " " << tris[3*i+2] << "\n";
}

struct Config {
  double L = 1.0, Rmax = 0.3;
  int nx = 8, ntheta = 12;
  double n0 = 1.4;
  double lr = 0.15;
  int max_iter = 6;
  std::string experiment = "run";
  RunConfig run;
};

int main_impl(int argc, char **argv)
{
  Config cfg;

  // Deck defaults matching this demo's own physical setup (distinct from
  // shape_case.h's 2D-Bezier-oriented RunConfig defaults).
  cfg.run.nsettle = 1500;
  cfg.run.navg = 3000;
  cfg.run.vstream = 1000.0;
  cfg.run.tinf = 300.0;
  cfg.run.wall_temp = 300.0;
  cfg.run.wall_accom = 1.0;
  cfg.run.nrho = 1.0e20;
  cfg.run.fnum = 5.0e15;
  cfg.run.tstep = 1e-6;
  cfg.run.collisions = 0;   // free-molecular: 3D collisional DSMC is far
                            // more expensive than 2D and not needed to
                            // demonstrate the shape gradient
  cfg.run.boundary = "o p p";
  cfg.run.stats_every = 200;
  // apex at (x0,0,0); y,z already centered on the axis at 0, so only x
  // is translated -- box: x in [0,4], y,z in [-1,1]
  cfg.run.origin[0] = 1.0; cfg.run.origin[1] = 0.0; cfg.run.origin[2] = 0.0;
  cfg.run.boxlo[0]  = 0.0; cfg.run.boxhi[0]  = 4.0;
  cfg.run.boxlo[1]  = -1.0; cfg.run.boxhi[1] = 1.0;
  cfg.run.boxlo[2]  = -1.0; cfg.run.boxhi[2] = 1.0;
  cfg.run.grid[0] = 20; cfg.run.grid[1] = 10; cfg.run.grid[2] = 10;

  int seed = 12345;

  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "--n0") && i+1 < argc) cfg.n0 = atof(argv[++i]);
    else if (!strcmp(argv[i], "--lr") && i+1 < argc) cfg.lr = atof(argv[++i]);
    else if (!strcmp(argv[i], "--max-iter") && i+1 < argc) cfg.max_iter = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--nsettle") && i+1 < argc) cfg.run.nsettle = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--navg") && i+1 < argc) cfg.run.navg = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--experiment") && i+1 < argc) cfg.experiment = argv[++i];
    else if (!strcmp(argv[i], "--verbose")) cfg.run.verbose = 1;
    else if (!strcmp(argv[i], "--help") || !strcmp(argv[i], "-h")) {
      std::printf("usage: power_law_opt [--n0 N] [--lr LR] [--max-iter K]\n"
                 "  [--nsettle N] [--navg N] [--experiment NAME] [--verbose]\n"
                 "3D power-law-body heatflux shape optimization (1 design var: n).\n"
                 "Run from a dir with N.species. Output in output/pl_<ad|fd>_<date>_<experiment>/.\n");
      return 0;
    } else { std::fprintf(stderr, "unknown arg: %s\n", argv[i]); return 1; }
  }

  if (access(cfg.run.species_file, R_OK) != 0) die("cannot read N.species from cwd");

  PowerLawShape shape(cfg.L, cfg.Rmax, cfg.nx, cfg.ntheta);
  HeatFluxObjective obj;

#ifdef SPARTA_AD
  std::string dir = make_output_dir("pl_ad", cfg.experiment);
#else
  std::string dir = make_output_dir("pl_fd", cfg.experiment);
#endif
  std::printf("power-law 3D heatflux optimization (%s gradient)\n",
#ifdef SPARTA_AD
             "AD"
#else
             "central finite-difference"
#endif
             );
  std::printf("  n0 = %.4g   output dir = %s\n", cfg.n0, dir.c_str());

  double lo, hi;
  shape.bounds(&lo, &hi);

  // save initial geometry
  {
    SurfMesh m;
    double alpha[1] = {cfg.n0};
    shape.to_mesh(alpha, m);
    write_geom_txt((dir + "/shape_before.txt").c_str(), m.pts.data(),
                  m.npoints(), m.elems.data(), m.nelems(), cfg.n0);
  }

  std::ofstream traj(dir + "/trajectory.csv");
  traj << "iter,n,value,grad\n";

  const double fd_h = 0.05;   // stock-build-only central FD step
  int seeds[1] = {seed};

  double n = cfg.n0;
  double value = 0.0, grad = 0.0;
  for (int it = 0; it <= cfg.max_iter; it++) {
    double alpha[1] = {n};
#ifdef SPARTA_AD
    value = evaluate_avg(shape, obj, alpha, seeds, 1, cfg.run, &grad);
#else
    value = grad_fd(shape, obj, alpha, seeds, 1, cfg.run, fd_h, &grad);
#endif
    traj << it << "," << n << "," << value << "," << grad << "\n";
    traj.flush();
    std::printf("  iter %d: n=%.5f  value=%.6e  grad=%.6e\n", it, n, value, grad);
    if (it == cfg.max_iter) break;
    // Normalized (sign-based) step: the raw gradient's magnitude reflects
    // the objective's own O(1e4) scale, not a usable step size in n-space
    // (O(1)); --lr is a fixed step size per iteration, direction only
    // from the gradient's sign.
    double step = (grad > 0.0) ? -cfg.lr : (grad < 0.0 ? cfg.lr : 0.0);
    double n_new = n + step;
    if (n_new < lo) n_new = lo;
    if (n_new > hi) n_new = hi;
    n = n_new;
  }

  // save final geometry
  {
    SurfMesh m;
    double alpha[1] = {n};
    shape.to_mesh(alpha, m);
    write_geom_txt((dir + "/shape_after.txt").c_str(), m.pts.data(),
                  m.npoints(), m.elems.data(), m.nelems(), n);
  }

  {
    std::ofstream cf(dir + "/config.txt");
    cf << "experiment   = " << cfg.experiment << "\n";
    cf << "objective    = heatflux (etot, diffuse wall, wall_accom="
       << cfg.run.wall_accom << ")\n";
#ifdef SPARTA_AD
    cf << "gradient     = AD (forward-mode, uncorrected -- see "
          "docs/AD_GRADIENTS.md; SPARTA_AD_SCORE_CORRECTION "
          "not set by this driver)\n";
#else
    cf << "gradient     = finite-difference (h=" << fd_h << ")\n";
#endif
    cf << "L            = " << cfg.L << "\n";
    cf << "Rmax         = " << cfg.Rmax << "\n";
    cf << "nx, ntheta   = " << cfg.nx << ", " << cfg.ntheta << "\n";
    cf << "n0           = " << cfg.n0 << "\n";
    cf << "n bounds     = [" << lo << ", " << hi << "]\n";
    cf << "lr           = " << cfg.lr << "\n";
    cf << "max_iter     = " << cfg.max_iter << "\n";
    cf << "nsettle      = " << cfg.run.nsettle << "\n";
    cf << "navg         = " << cfg.run.navg << "\n";
    cf << "vstream      = " << cfg.run.vstream << "\n";
    cf << "wall_temp    = " << cfg.run.wall_temp << "\n";
    cf << "wall_accom   = " << cfg.run.wall_accom << "\n";
    cf << "collisions   = 0 (free-molecular)\n";
    cf << "seed         = " << seed << "\n";
  }
  {
    std::ofstream rf(dir + "/result.txt");
    rf << "experiment    : " << cfg.experiment << "\n";
    rf << "objective     : heatflux\n";
    rf << "shape         : power_law_body (ndesign=1)\n";
#ifdef SPARTA_AD
    rf << "gradient      : AD (forward-mode, uncorrected)\n";
#else
    rf << "gradient      : finite-difference (h=" << fd_h << ")\n";
#endif
    rf << "n start       : " << cfg.n0 << "\n";
    rf << "n final       : " << n << "\n";
    rf << "value final   : " << value << "\n";
    rf << "grad final    : " << grad << "\n";
    rf << "results in    : " << dir << "/\n";
  }

  std::printf("done: n %.5f -> %.5f, value=%.6e\n  results in %s/\n",
             cfg.n0, n, value, dir.c_str());
  return 0;
}

}  // namespace

int main(int argc, char **argv)
{
  try {
    return main_impl(argc, argv);
  } catch (const std::exception &e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
}
