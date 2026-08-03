/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer

   power_law_main.cpp: a 3D heat-flux shape-optimization demo on a
   PowerLawBody (see power_law_body.h) -- a single design variable n
   controlling nose bluntness (n<1 blunt/flaring, n=1 cone, n>1 pointy
   spike). Diffuse wall (wall_accom=1) so HeatFluxObjective (etot) is
   nondegenerate, same convention as the 2D Bezier heatflux runs.

   Deliberately a standalone driver, not built on ShapeTNLP/IPOPT: with
   ndesign()==1 a full NLP solver is unnecessary machinery. Instead a
   plain gradient-descent loop (AD gradient every step, one SPARTA run
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

#include "power_law_body.h"
#include "heat_flux_objective.h"
#include "run_output.h"
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

void write_surf3d_tmp(const char *path, const double *pts, int npt,
                      const int *tris, int ntri)
{
  FILE *fp = std::fopen(path, "w");
  if (!fp) die("cannot open surf tmp file for writing");
  std::fprintf(fp, "power_law_main generated 3D surf\n\n");
  std::fprintf(fp, "%d points\n%d triangles\n\nPoints\n\n", npt, ntri);
  for (int i = 0; i < npt; i++)
    std::fprintf(fp, "%d %.15g %.15g %.15g\n", i + 1,
                pts[3 * i], pts[3 * i + 1], pts[3 * i + 2]);
  std::fprintf(fp, "\nTriangles\n\n");
  for (int i = 0; i < ntri; i++)
    std::fprintf(fp, "%d %d %d %d\n", i + 1,
                tris[3 * i] + 1, tris[3 * i + 1] + 1, tris[3 * i + 2] + 1);
  std::fclose(fp);
}

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

#ifdef SPARTA_AD
void write_jac3d_tmp(const char *path, const double *jac, int npt)
{
  FILE *fp = std::fopen(path, "w");
  if (!fp) die("cannot open jacobian tmp file for writing");
  for (int i = 0; i < npt; i++)
    std::fprintf(fp, "%.15g %.15g %.15g\n",
                jac[3 * i], jac[3 * i + 1], jac[3 * i + 2]);
  std::fclose(fp);
}
#endif

struct Config {
  double L = 1.0, Rmax = 0.3;
  int nx = 8, ntheta = 12;
  double n0 = 1.4;
  double lr = 0.15;
  int max_iter = 6;
  int nsettle = 1500, navg = 3000;
  double vstream = 1000.0, tinf = 300.0, wall_temp = 300.0, wall_accom = 1.0;
  double nrho = 1.0e20, fnum = 5.0e15, tstep = 1e-6;
  double x0 = 1.0;                          // apex position in box
  double boxhi[3] = {4.0, 1.0, 1.0};         // box: x in [0,boxhi0], y,z in [-boxhi1,boxhi1]
  int grid[3] = {20, 10, 10};
  const char *species_file = "N.species";
  const char *species_names = "N";
  std::string experiment = "run";
  int verbose = 0;
};

// One SPARTA run: build geometry at n, run, extract heatflux + dHeatflux/dn.
double evaluate(const PowerLawBody &body, const HeatFluxObjective &obj,
                double n, int seed, const Config &c, double *grad)
{
  int npt = body.npoints();
  int ntri = body.ntris();
  std::vector<double> pts(3 * npt);
  std::vector<int> tris(3 * ntri);
  double alpha[1] = {n};
  body.to_tris(alpha, pts.data(), tris.data());

  // translate into the box: apex at (x0, ybox_mid, zbox_mid) = (x0,0,0)
  // since y,z already centered on the axis at 0.
  for (int i = 0; i < npt; i++) {
    pts[3*i + 0] += c.x0;
    if (pts[3*i+0] <= 0.0 || pts[3*i+0] >= c.boxhi[0] ||
        std::fabs(pts[3*i+1]) >= c.boxhi[1] || std::fabs(pts[3*i+2]) >= c.boxhi[2])
      die("body outside box");
  }

  char surfpath[256];
  std::snprintf(surfpath, sizeof(surfpath), "tmp_pl_surf_%d_%d.data",
               (int) getpid(), seed);
  write_surf3d_tmp(surfpath, pts.data(), npt, tris.data(), ntri);

#ifdef SPARTA_AD
  char jacpath[256];
  std::snprintf(jacpath, sizeof(jacpath), "tmp_pl_jac_%d_%d.data",
               (int) getpid(), seed);
  {
    std::vector<double> jac(3 * npt);
    body.jacobian(alpha, jac.data());
    write_jac3d_tmp(jacpath, jac.data(), npt);
  }
#endif

  void *spa = open_sparta(c.verbose);

  char line[512];
  std::snprintf(line, sizeof(line), "seed %d", seed);
  cmd(spa, line);
  cmd(spa, "dimension 3");
  cmd(spa, "global gridcut 0.0 comm/sort yes");
  cmd(spa, "boundary o p p");
  std::snprintf(line, sizeof(line), "create_box 0 %.15g %.15g %.15g %.15g %.15g",
               c.boxhi[0], -c.boxhi[1], c.boxhi[1], -c.boxhi[2], c.boxhi[2]);
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
  setenv("SPARTA_AD_SEED_JACFILE", jacpath, 1);
#endif
  std::snprintf(line, sizeof(line), "read_surf %s", surfpath);
  cmd(spa, line);

  std::snprintf(line, sizeof(line), "surf_collide 1 diffuse %.15g %.15g",
               c.wall_temp, c.wall_accom);
  cmd(spa, line);
  cmd(spa, "surf_modify all collide 1");

  // free-molecular (no `collide vss`): 3D collisional DSMC is far more
  // expensive than 2D and not needed to demonstrate the shape gradient.
  cmd(spa, "fix in emit/face gas xlo twopass");
  std::snprintf(line, sizeof(line), "timestep %.15g", c.tstep);
  cmd(spa, line);
  if (c.verbose) cmd(spa, "stats 200");

  std::snprintf(line, sizeof(line), "run %d", c.nsettle);
  cmd(spa, line);
  cmd(spa, "reset_timestep 0");

  obj.setup(spa, c.navg);
  std::snprintf(line, sizeof(line), "run %d", c.navg);
  cmd(spa, line);

  double value = obj.extract(spa, 1, grad);

  close_sparta(spa);
  unlink(surfpath);
#ifdef SPARTA_AD
  unlink(jacpath);
#endif
  return value;
}

int main_impl(int argc, char **argv)
{
  Config cfg;
  int seed = 12345;

  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "--n0") && i+1 < argc) cfg.n0 = atof(argv[++i]);
    else if (!strcmp(argv[i], "--lr") && i+1 < argc) cfg.lr = atof(argv[++i]);
    else if (!strcmp(argv[i], "--max-iter") && i+1 < argc) cfg.max_iter = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--nsettle") && i+1 < argc) cfg.nsettle = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--navg") && i+1 < argc) cfg.navg = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--experiment") && i+1 < argc) cfg.experiment = argv[++i];
    else if (!strcmp(argv[i], "--verbose")) cfg.verbose = 1;
    else if (!strcmp(argv[i], "--help") || !strcmp(argv[i], "-h")) {
      std::printf("usage: power_law_opt [--n0 N] [--lr LR] [--max-iter K]\n"
                 "  [--nsettle N] [--navg N] [--experiment NAME] [--verbose]\n"
                 "3D power-law-body heatflux shape optimization (1 design var: n).\n"
                 "Run from a dir with N.species. Output in output/pl_<ad|fd>_<date>_<experiment>/.\n");
      return 0;
    } else { std::fprintf(stderr, "unknown arg: %s\n", argv[i]); return 1; }
  }

  if (access(cfg.species_file, R_OK) != 0) die("cannot read N.species from cwd");

  PowerLawBody body(cfg.L, cfg.Rmax, cfg.nx, cfg.ntheta);
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
  body.bounds(&lo, &hi);

  // save initial geometry
  {
    int npt = body.npoints(), ntri = body.ntris();
    std::vector<double> pts(3*npt);
    std::vector<int> tris(3*ntri);
    double alpha[1] = {cfg.n0};
    body.to_tris(alpha, pts.data(), tris.data());
    write_geom_txt((dir + "/shape_before.txt").c_str(), pts.data(), npt,
                  tris.data(), ntri, cfg.n0);
  }

  std::ofstream traj(dir + "/trajectory.csv");
  traj << "iter,n,value,grad\n";

  const double fd_h = 0.05;   // stock-build-only central FD step

  double n = cfg.n0;
  double value = 0.0, grad = 0.0;
  for (int it = 0; it <= cfg.max_iter; it++) {
#ifdef SPARTA_AD
    value = evaluate(body, obj, n, seed, cfg, &grad);
#else
    value = evaluate(body, obj, n, seed, cfg, nullptr);
    double vp = evaluate(body, obj, n + fd_h, seed, cfg, nullptr);
    double vm = evaluate(body, obj, n - fd_h, seed, cfg, nullptr);
    grad = (vp - vm) / (2.0 * fd_h);
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
    int npt = body.npoints(), ntri = body.ntris();
    std::vector<double> pts(3*npt);
    std::vector<int> tris(3*ntri);
    double alpha[1] = {n};
    body.to_tris(alpha, pts.data(), tris.data());
    write_geom_txt((dir + "/shape_after.txt").c_str(), pts.data(), npt,
                  tris.data(), ntri, n);
  }

  {
    std::ofstream cf(dir + "/config.txt");
    cf << "experiment   = " << cfg.experiment << "\n";
    cf << "objective    = heatflux (etot, diffuse wall, wall_accom="
       << cfg.wall_accom << ")\n";
#ifdef SPARTA_AD
    cf << "gradient     = AD (forward-mode, uncorrected -- see "
          "docs/ad_phase_c_investigation/FINDINGS.md; SPARTA_AD_SCORE_CORRECTION "
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
    cf << "nsettle      = " << cfg.nsettle << "\n";
    cf << "navg         = " << cfg.navg << "\n";
    cf << "vstream      = " << cfg.vstream << "\n";
    cf << "wall_temp    = " << cfg.wall_temp << "\n";
    cf << "wall_accom   = " << cfg.wall_accom << "\n";
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
