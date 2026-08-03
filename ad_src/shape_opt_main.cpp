/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer

   shape_opt_main: IPOPT-driven shape-optimization driver, dimension-
   agnostic. --shape bezier|powerlaw selects a 2D or 3D Shape;
   --objective drag|heatflux selects the QOI. Gradient source depends on
   which binary this is: shape_opt (stock) always uses grad_fd();
   shape_opt_ad (built with -DSPARTA_AD) gets the whole gradient from one
   SPARTA run instead -- see shape_tnlp.h for both, including the
   important caveat that AD gradients are used exactly as the solver
   returns them, with a known low bias, not corrected here (see
   docs/AD_GRADIENTS.md and the --score-correction flag below).

   Absorbs what used to be power_law_main.cpp's own standalone gradient-
   descent driver: with a 3D PowerLawShape now going through the same
   Shape/ShapeTNLP/Constraint machinery as the 2D Bezier case, a second
   driver isn't needed. Tradeoff accepted explicitly: 3D optimization now
   requires IPOPT, where the old standalone driver had a built-in
   gradient-descent fallback; in exchange it gets IPOPT's objective
   scaling, the principled fix for the old driver's hand-rolled
   "normalized step" hack (raw gradient O(1e4) vs design variable O(1)).

   Bounds come from Shape::bounds(), not CLI flags -- the reference's
   --x-lo/--x-hi/--y-lo/--y-hi are Bezier-specific naming baked into
   what's meant to be a generic driver.

   Optional general constraints: --min-area A imposes |measure(alpha)|
   >= A via MinSizeConstraint, so the optimizer can't collapse the body
   to a needle (2D: area) or a sliver (3D: volume) chasing a smaller
   objective. Uses an exact analytic gradient chained through
   Shape::jacobian(), not FD/AD -- see min_size_constraint.cpp.

   Results are written to a dated output folder:
     output/<fd|ad>_<YYYY-MM-DD>_<experiment>[_N]/
       config.txt         resolved parameters
       ipopt.log          IPOPT's own iteration log
       trajectory.csv     per-iteration alpha / value / infeasibility
       result.txt         start & final alpha, value, status, size measure
       shape_initial.txt  initial mesh (points + line/triangle indices)
       shape_final.txt    final mesh, same format
       shapes.svg         initial (faint) vs optimized (bold) outline
                           -- only written when the shape is 2D

   IMPORTANT for --objective heatflux: total heat flux needs a
   thermally-accommodating wall (--wall-accom > 0) and enough averaging
   to see past its noise floor (energy is a higher/quadratic velocity
   moment than momentum, so it's inherently noisier than drag at the
   same sample size). This driver does not auto-inject those for the
   bezier shape; pass them explicitly, e.g. --wall-accom 1.0 --navg 4000
   (--shape powerlaw defaults wall_accom to 1.0 already, matching the
   old standalone driver's own default).

   Run from ad_src/ (where N.species/N.vss live):
     ./build/shape_opt --alpha 1.3,1.0,2.7,0.8 --experiment baseline
     ./build/shape_opt --shape powerlaw --objective heatflux \
       --alpha 1.4 --experiment pl_baseline
------------------------------------------------------------------------- */

#include "bezier_shape.h"
#include "cli.h"
#include "drag_objective.h"
#include "heat_flux_objective.h"
#include "min_size_constraint.h"
#include "power_law_shape.h"
#include "run_output.h"
#include "shape_case.h"
#include "shape_tnlp.h"
#include "svg_shape.h"

#include "IpIpoptApplication.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <exception>
#include <fstream>
#include <string>
#include <unistd.h>
#include <vector>

using namespace Ipopt;

static const char *USAGE =
  "usage: shape_opt [--input FILE] [options]\n"
  "  --objective drag|heatflux   (default: drag)\n"
  "  --shape bezier|powerlaw     (default: bezier)\n"
  "  --alpha x1,y1,x2,y2         bezier: 4 values. powerlaw: 1 (n).\n"
  "                              (required unless set in --input)\n"
  "  --seeds s1,s2,...           (default: 12345,67890,13579 -- more\n"
  "                               than shape_main's single-seed default,\n"
  "                               since a real optimization loop needs a\n"
  "                               less noisy gradient than a one-off eval)\n"
  "  --max-iter N  --tol T  --acceptable-tol T  --acceptable-iter N\n"
  "  --h H                       FD step (default: 0.05)\n"
  "  --experiment NAME\n"
  "  --nsettle N  --navg N  --vstream V\n"
  "  --chord L  --nseg N          bezier-only discretization\n"
  "  --pl-L L  --pl-rmax R  --pl-nx N  --pl-ntheta N   powerlaw-only mesh\n"
  "  --min-area A                 constrain |measure(alpha)| >= A (area\n"
  "                                for a 2D shape, volume for 3D). Default: none.\n"
  "  --wall-temp T  --wall-accom A   (diffuse wall; see note above for\n"
  "                                   --objective heatflux)\n"
  "  --specular  --nocoll  --verbose\n"
  "  --score-correction          AD build only: enable the flux-measure\n"
  "                               score-function correction in\n"
  "                               compute_surf.cpp (see docs/AD_GRADIENTS.md).\n"
  "                               No effect on stock builds. Off by default.\n"
  "\n"
  "  Minimizes the chosen objective with IPOPT, subject to the shape's\n"
  "  own box bounds (Shape::bounds()). Run from a dir with\n"
  "  N.species/N.vss. Output written to output/<fd|ad>_<date>_<experiment>/.\n";

// ---- config ------------------------------------------------------------

struct Config {
  std::string objective_name = "drag";
  std::string shape_name = "bezier";
  std::vector<double> alpha;
  std::vector<int> seeds;
  // objective is a noisy Monte-Carlo estimate, so tight tolerances are
  // unreachable: IPOPT would grind with vanishing line-search steps
  // hunting a precision the physics can't deliver. Defaults are chosen
  // for the stochastic objective; acc_tol/acc_iter are the fallback that
  // accepts a noise-level optimum. (Same rationale, same numbers, as the
  // reference's DragCase.)
  double h = 0.05, tol = 1e-4, acc_tol = 1e-3;
  int max_iter = 40, acc_iter = 5;
  std::string experiment = "run";
  double chord = 4.0;    // BezierShape ctor arg
  int nseg = 25;
  double pl_L = 1.0, pl_Rmax = 0.3;   // PowerLawShape ctor args
  int pl_nx = 8, pl_ntheta = 12;
  RunConfig run;
  bool alpha_set = false;
  bool min_area_set = false;
  double min_area = 0.0;
};

// A --shape powerlaw run needs a different deck "profile" (box, grid,
// gas scales, boundary conditions -- see shape_case.h's RunConfig
// comment) than the 2D-Bezier-oriented defaults Config::run starts with.
// Applied BEFORE the real argument parsing below (which can still
// override any individual field), matching what the old standalone
// power_law_main.cpp hardcoded.
static void apply_powerlaw_profile(Config &cfg)
{
  cfg.objective_name = "heatflux";
  cfg.run.nsettle = 1500;
  cfg.run.navg = 3000;
  cfg.run.vstream = 1000.0;
  cfg.run.tinf = 300.0;
  cfg.run.wall_temp = 300.0;
  cfg.run.wall_accom = 1.0;
  cfg.run.nrho = 1.0e20;
  cfg.run.fnum = 5.0e15;
  cfg.run.tstep = 1e-6;
  cfg.run.collisions = 0;    // free-molecular
  cfg.run.boundary = "o p p";
  cfg.run.stats_every = 200;
  cfg.run.origin[0] = 1.0; cfg.run.origin[1] = 0.0; cfg.run.origin[2] = 0.0;
  cfg.run.boxlo[0]  = 0.0; cfg.run.boxhi[0]  = 4.0;
  cfg.run.boxlo[1]  = -1.0; cfg.run.boxhi[1] = 1.0;
  cfg.run.boxlo[2]  = -1.0; cfg.run.boxhi[2] = 1.0;
  cfg.run.grid[0] = 20; cfg.run.grid[1] = 10; cfg.run.grid[2] = 10;
}

// Looks only for the shape choice (--shape / --input's "shape=" key,
// same precedence as the real parse below) so apply_powerlaw_profile()
// can run before any other field is parsed -- otherwise it would
// clobber whatever the user explicitly set.
static std::string prescan_shape_name(int narg, char **arg)
{
  std::string shape_name = "bezier";
  for (int i = 1; i < narg; i++) {
    if (!strcmp(arg[i], "--input") && i + 1 < narg) {
      std::ifstream f(arg[i + 1]);
      std::string line;
      while (std::getline(f, line)) {
        size_t hash = line.find('#');
        if (hash != std::string::npos) line = line.substr(0, hash);
        size_t eq = line.find('=');
        if (eq == std::string::npos) continue;
        std::string key = trim(line.substr(0, eq));
        std::string val = trim(line.substr(eq + 1));
        if (key == "shape") shape_name = val;
      }
    }
  }
  for (int i = 1; i < narg; i++)
    if (!strcmp(arg[i], "--shape") && i + 1 < narg) shape_name = arg[i + 1];
  return shape_name;
}

// ---- input-file parsing (drag.in-style, '#' starts a comment) -----------

static bool parse_input_file(const char *path, Config &cfg)
{
  std::ifstream f(path);
  if (!f) { fprintf(stderr, "ERROR: cannot open input file '%s'\n", path); return false; }
  std::string line;
  while (std::getline(f, line)) {
    size_t hash = line.find('#');
    if (hash != std::string::npos) line = line.substr(0, hash);
    size_t eq = line.find('=');
    if (eq == std::string::npos) continue;
    std::string key = trim(line.substr(0, eq));
    std::string val = trim(line.substr(eq + 1));
    if (key.empty() || val.empty()) continue;

    if (key == "objective") { cfg.objective_name = val;
    } else if (key == "shape") { cfg.shape_name = val;
    } else if (key == "alpha") {
      cfg.alpha = parse_doubles(val);
      if (!cfg.alpha.empty()) cfg.alpha_set = true;
    } else if (key == "seeds") { cfg.seeds = parse_ints(val);
    } else if (key == "max_iter") { cfg.max_iter = atoi(val.c_str());
    } else if (key == "tol") { cfg.tol = atof(val.c_str());
    } else if (key == "acceptable_tol") { cfg.acc_tol = atof(val.c_str());
    } else if (key == "acceptable_iter") { cfg.acc_iter = atoi(val.c_str());
    } else if (key == "h") { cfg.h = atof(val.c_str());
    } else if (key == "experiment") { cfg.experiment = val;
    } else if (key == "min_area") { cfg.min_area = atof(val.c_str()); cfg.min_area_set = true;
    } else if (key == "nseg") { cfg.nseg = atoi(val.c_str());
    } else if (key == "chord") { cfg.chord = atof(val.c_str());
    } else if (key == "pl_L") { cfg.pl_L = atof(val.c_str());
    } else if (key == "pl_rmax") { cfg.pl_Rmax = atof(val.c_str());
    } else if (key == "pl_nx") { cfg.pl_nx = atoi(val.c_str());
    } else if (key == "pl_ntheta") { cfg.pl_ntheta = atoi(val.c_str());
    } else if (key == "nsettle") { cfg.run.nsettle = atoi(val.c_str());
    } else if (key == "navg") { cfg.run.navg = atoi(val.c_str());
    } else if (key == "vstream") { cfg.run.vstream = atof(val.c_str());
    } else if (key == "wall_temp") { cfg.run.wall_temp = atof(val.c_str());
    } else if (key == "wall_accom") { cfg.run.wall_accom = atof(val.c_str());
    } else if (key == "specular") { cfg.run.specular = atoi(val.c_str());
    } else if (key == "collisions") { cfg.run.collisions = atoi(val.c_str());
    } else if (key == "verbose") { cfg.run.verbose = atoi(val.c_str());
    } else if (key == "score_correction") { cfg.run.score_correction = atoi(val.c_str()) != 0;
    }
    // unknown keys ignored
  }
  return true;
}

// ---- generic mesh dump (both dims) --------------------------------------

static void write_mesh_txt(const std::string &path, const Shape &shape,
                           const double *alpha)
{
  SurfMesh m;
  shape.to_mesh(alpha, m);
  std::ofstream f(path);
  f << "# dim=" << m.dim << " npoints=" << m.npoints()
    << " nelems=" << m.nelems() << "\n";
  f << "POINTS\n";
  for (int i = 0; i < m.npoints(); i++) {
    for (int d = 0; d < m.dim; d++)
      f << m.pts[m.dim * i + d] << (d + 1 < m.dim ? " " : "\n");
  }
  f << ((m.dim == 2) ? "LINES\n" : "TRIANGLES\n");
  for (int i = 0; i < m.nelems(); i++) {
    for (int d = 0; d < m.dim; d++)
      f << m.elems[m.dim * i + d] << (d + 1 < m.dim ? " " : "\n");
  }
}

static const char *status_str(ApplicationReturnStatus s)
{
  switch (s) {
    case Solve_Succeeded:               return "Solve_Succeeded";
    case Solved_To_Acceptable_Level:    return "Solved_To_Acceptable_Level";
    case Infeasible_Problem_Detected:   return "Infeasible_Problem_Detected";
    case Search_Direction_Becomes_Too_Small:
                                        return "Search_Direction_Too_Small";
    case Diverging_Iterates:            return "Diverging_Iterates";
    case User_Requested_Stop:           return "User_Requested_Stop";
    case Feasible_Point_Found:          return "Feasible_Point_Found";
    case Maximum_Iterations_Exceeded:   return "Maximum_Iterations_Exceeded";
    case Restoration_Failed:            return "Restoration_Failed";
    case Error_In_Step_Computation:     return "Error_In_Step_Computation";
    case Maximum_CpuTime_Exceeded:      return "Maximum_CpuTime_Exceeded";
    case Invalid_Option:                return "Invalid_Option";
    case Invalid_Number_Detected:       return "Invalid_Number_Detected";
    default:                            return "Error/Other";
  }
}

// ---- main -----------------------------------------------------------------

int main_impl(int narg, char **arg)
{
  Config cfg;

  if (prescan_shape_name(narg, arg) == "powerlaw") apply_powerlaw_profile(cfg);

  // pass 1: --input FILE (so CLI flags below can still override it)
  for (int i = 1; i < narg; i++) {
    if (!strcmp(arg[i], "--input") && i + 1 < narg) {
      if (!parse_input_file(arg[i + 1], cfg)) return 1;
    } else if (!strcmp(arg[i], "--help") || !strcmp(arg[i], "-h")) {
      printf("%s", USAGE); return 0;
    }
  }

  // pass 2: CLI overrides
  for (int i = 1; i < narg; i++) {
    if (!strcmp(arg[i], "--input")) { i++; continue; }
    else if (!strcmp(arg[i], "--objective") && i + 1 < narg) { cfg.objective_name = arg[++i];
    } else if (!strcmp(arg[i], "--shape") && i + 1 < narg) { cfg.shape_name = arg[++i];
    } else if (!strcmp(arg[i], "--alpha") && i + 1 < narg) {
      cfg.alpha = parse_doubles(arg[++i]);
      if (!cfg.alpha.empty()) cfg.alpha_set = true;
    } else if (!strcmp(arg[i], "--seeds") && i + 1 < narg) {
      cfg.seeds = parse_ints(arg[++i]);
    } else if (!strcmp(arg[i], "--max-iter") && i + 1 < narg) { cfg.max_iter = atoi(arg[++i]);
    } else if (!strcmp(arg[i], "--tol") && i + 1 < narg) { cfg.tol = atof(arg[++i]);
    } else if (!strcmp(arg[i], "--acceptable-tol") && i + 1 < narg) { cfg.acc_tol = atof(arg[++i]);
    } else if (!strcmp(arg[i], "--acceptable-iter") && i + 1 < narg) { cfg.acc_iter = atoi(arg[++i]);
    } else if (!strcmp(arg[i], "--h") && i + 1 < narg) { cfg.h = atof(arg[++i]);
    } else if (!strcmp(arg[i], "--experiment") && i + 1 < narg) { cfg.experiment = arg[++i];
    } else if (!strcmp(arg[i], "--min-area") && i + 1 < narg) {
      cfg.min_area = atof(arg[++i]); cfg.min_area_set = true;
    } else if (!strcmp(arg[i], "--nseg") && i + 1 < narg) { cfg.nseg = atoi(arg[++i]);
    } else if (!strcmp(arg[i], "--chord") && i + 1 < narg) { cfg.chord = atof(arg[++i]);
    } else if (!strcmp(arg[i], "--pl-L") && i + 1 < narg) { cfg.pl_L = atof(arg[++i]);
    } else if (!strcmp(arg[i], "--pl-rmax") && i + 1 < narg) { cfg.pl_Rmax = atof(arg[++i]);
    } else if (!strcmp(arg[i], "--pl-nx") && i + 1 < narg) { cfg.pl_nx = atoi(arg[++i]);
    } else if (!strcmp(arg[i], "--pl-ntheta") && i + 1 < narg) { cfg.pl_ntheta = atoi(arg[++i]);
    } else if (!strcmp(arg[i], "--nsettle") && i + 1 < narg) { cfg.run.nsettle = atoi(arg[++i]);
    } else if (!strcmp(arg[i], "--navg") && i + 1 < narg) { cfg.run.navg = atoi(arg[++i]);
    } else if (!strcmp(arg[i], "--vstream") && i + 1 < narg) { cfg.run.vstream = atof(arg[++i]);
    } else if (!strcmp(arg[i], "--wall-temp") && i + 1 < narg) { cfg.run.wall_temp = atof(arg[++i]);
    } else if (!strcmp(arg[i], "--wall-accom") && i + 1 < narg) { cfg.run.wall_accom = atof(arg[++i]);
    } else if (!strcmp(arg[i], "--specular")) { cfg.run.specular = 1;
    } else if (!strcmp(arg[i], "--score-correction")) { cfg.run.score_correction = true;
    } else if (!strcmp(arg[i], "--nocoll")) { cfg.run.collisions = 0;
    } else if (!strcmp(arg[i], "--verbose")) { cfg.run.verbose = 1;
    } else if (!strcmp(arg[i], "--help") || !strcmp(arg[i], "-h")) { printf("%s", USAGE); return 0;
    } else { fprintf(stderr, "unknown/incomplete arg: %s\n%s", arg[i], USAGE); return 1; }
  }

  if (!cfg.alpha_set) { fprintf(stderr, "ERROR: --alpha required\n%s", USAGE); return 1; }
  if (cfg.seeds.empty()) cfg.seeds = {12345, 67890, 13579};

  BezierShape bezier_shape(cfg.chord, cfg.nseg);
  PowerLawShape powerlaw_shape(cfg.pl_L, cfg.pl_Rmax, cfg.pl_nx, cfg.pl_ntheta);
  Shape *shape = 0;
  if (cfg.shape_name == "bezier") shape = &bezier_shape;
  else if (cfg.shape_name == "powerlaw") shape = &powerlaw_shape;
  else { fprintf(stderr, "unknown --shape '%s'\n", cfg.shape_name.c_str()); return 1; }

  if ((int) cfg.alpha.size() != shape->ndesign()) {
    fprintf(stderr, "--alpha needs %d values for shape '%s', got %zu\n",
            shape->ndesign(), cfg.shape_name.c_str(), cfg.alpha.size());
    return 1;
  }

  DragObjective drag_obj;
  HeatFluxObjective hf_obj;
  Objective *obj = 0;
  if (cfg.objective_name == "drag") obj = &drag_obj;
  else if (cfg.objective_name == "heatflux") obj = &hf_obj;
  else { fprintf(stderr, "unknown --objective '%s'\n", cfg.objective_name.c_str()); return 1; }

  if (cfg.objective_name == "heatflux" && cfg.run.wall_accom <= 0.0) {
    fprintf(stderr,
            "WARNING: --objective heatflux with wall_accom=%.3g (default 0)"
            " -- total heat flux is nearly degenerate without wall thermal"
            " accommodation. Pass --wall-accom 1.0 (or similar) for a"
            " usable gradient. Proceeding anyway.\n", cfg.run.wall_accom);
  }

  int n = shape->ndesign();
  const bool show_bar = (!cfg.run.verbose) && isatty(fileno(stderr));

  std::vector<double> xl(n), xu(n);
  shape->bounds(xl.data(), xu.data());

  // sanity on bounds
  for (int j = 0; j < n; j++) {
    if (xl[j] >= xu[j]) {
      fprintf(stderr, "ERROR: shape bounds invalid at component %d "
              "(lo=%g >= hi=%g)\n", j, xl[j], xu[j]);
      return 1;
    }
  }

  // Objective scaling: the objective can be as small as ~1e-21 (2D
  // drag) or as large as ~1e5 (3D heatflux), so without scaling IPOPT's
  // log-barrier term either swamps it or is swamped by it. Scale by
  // 1/|value(start)| so the scaled objective is O(1). One evaluation at
  // the start point.
  double f0 = evaluate_avg(*shape, *obj, cfg.alpha.data(), cfg.seeds.data(),
                           (int) cfg.seeds.size(), cfg.run);
  double obj_scale = 1.0 / std::max(fabs(f0), 1e-300);

  // output folder + config echo. Prefix reflects the actual gradient
  // method so AD and FD runs land in distinguishable, non-colliding
  // places.
#ifdef SPARTA_AD
  std::string dir = make_output_dir("ad", cfg.experiment);
#else
  std::string dir = make_output_dir("fd", cfg.experiment);
#endif

  char timestamp[64];
  { time_t t = time(0); strftime(timestamp, sizeof(timestamp),
                                 "%Y-%m-%d %H:%M:%S", localtime(&t)); }

  {
    std::ofstream cf((dir + "/config.txt").c_str());
    cf << "timestamp    = " << timestamp << "\n";
    cf << "experiment   = " << cfg.experiment << "\n";
    cf << "objective    = " << cfg.objective_name << "\n";
    cf << "shape        = " << cfg.shape_name << "\n";
#ifdef SPARTA_AD
    if (cfg.run.score_correction)
      cf << "gradient     = AD (forward-mode, score-corrected -- see "
            "docs/AD_GRADIENTS.md)\n";
    else
      cf << "gradient     = AD (forward-mode, UNCORRECTED -- see "
            "docs/AD_GRADIENTS.md)\n";
#else
    cf << "gradient     = finite-difference (CRN)\n";
#endif
    cf << "value(start) = " << f0 << "\n";
    cf << "obj_scaling  = " << obj_scale << "   (1/|value(start)|)\n";
    cf << "alpha_start  = ";
    for (int j = 0; j < n; j++) cf << cfg.alpha[j] << (j + 1 < n ? ", " : "\n");
    cf << "bounds       = ";
    for (int j = 0; j < n; j++)
      cf << "[" << xl[j] << "," << xu[j] << "]" << (j + 1 < n ? "  " : "\n");
    cf << "seeds        = ";
    for (size_t k = 0; k < cfg.seeds.size(); k++)
      cf << cfg.seeds[k] << (k + 1 < cfg.seeds.size() ? ", " : "");
    cf << "\n";
    cf << "max_iter     = " << cfg.max_iter << "\n";
    cf << "tol          = " << cfg.tol << "\n";
    cf << "acceptable_tol  = " << cfg.acc_tol << "\n";
    cf << "acceptable_iter = " << cfg.acc_iter << "\n";
    cf << "h (FD step)  = " << cfg.h << "\n";
    if (cfg.shape_name == "bezier") {
      cf << "chord        = " << cfg.chord << "\n";
      cf << "nseg         = " << cfg.nseg << "\n";
    } else if (cfg.shape_name == "powerlaw") {
      cf << "pl_L         = " << cfg.pl_L << "\n";
      cf << "pl_Rmax      = " << cfg.pl_Rmax << "\n";
      cf << "pl_nx, ntheta= " << cfg.pl_nx << ", " << cfg.pl_ntheta << "\n";
    }
    cf << "nsettle      = " << cfg.run.nsettle << "\n";
    cf << "navg         = " << cfg.run.navg << "\n";
    cf << "vstream      = " << cfg.run.vstream << "\n";
    cf << "wall_temp    = " << cfg.run.wall_temp << "\n";
    cf << "wall_accom   = " << cfg.run.wall_accom << "\n";
    cf << "specular     = " << cfg.run.specular << "\n";
    cf << "collisions   = " << cfg.run.collisions << "\n";
    cf << "score_correction = " << (cfg.run.score_correction ? 1 : 0) << "\n";
    if (cfg.min_area_set) cf << "min_area     = " << cfg.min_area << "\n";
    else cf << "min_area     = (none)\n";
  }

  printf("%s optimization (%s gradient)\n", cfg.objective_name.c_str(),
#ifdef SPARTA_AD
        "AD"
#else
        "finite-difference"
#endif
        );
  printf("  start alpha = ");
  for (int j = 0; j < n; j++) printf("%.5g ", cfg.alpha[j]);
  printf("\n  output dir  = %s\n", dir.c_str());
  fflush(stdout);

  MinSizeConstraint min_size_constraint(cfg.min_area);
  std::vector<const Constraint *> constraints;
  if (cfg.min_area_set) constraints.push_back(&min_size_constraint);

  SmartPtr<ShapeTNLP> nlp = new ShapeTNLP(*shape, *obj, cfg.alpha.data(),
                                          xl.data(), xu.data(),
                                          cfg.seeds.data(), (int) cfg.seeds.size(),
                                          cfg.run, cfg.h, cfg.max_iter, show_bar,
                                          obj_scale, constraints);

  SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
  app->Options()->SetStringValue("hessian_approximation", "limited-memory");
  app->Options()->SetStringValue("mu_strategy", "adaptive");
  app->Options()->SetStringValue("nlp_scaling_method", "user-scaling");
  app->Options()->SetNumericValue("tol", cfg.tol);
  app->Options()->SetNumericValue("acceptable_tol", cfg.acc_tol);
  app->Options()->SetIntegerValue("acceptable_iter", cfg.acc_iter);
  app->Options()->SetIntegerValue("max_iter", cfg.max_iter);
  app->Options()->SetStringValue("sb", "yes");                  // no banner
  app->Options()->SetIntegerValue("print_level", cfg.run.verbose ? 5 : 0);
  app->Options()->SetStringValue("output_file", dir + "/ipopt.log");
  app->Options()->SetIntegerValue("file_print_level", 5);

  ApplicationReturnStatus st = app->Initialize();
  if (st != Solve_Succeeded) {
    fprintf(stderr, "ERROR: IPOPT initialization failed (%s)\n", status_str(st));
    return 1;
  }

  st = app->OptimizeTNLP(nlp);

  // trajectory.csv
  {
    std::ofstream tf((dir + "/trajectory.csv").c_str());
    tf << "iter";
    for (int j = 0; j < n; j++) tf << ",alpha" << j;
    tf << ",value,inf_pr,inf_du\n";
    tf.setf(std::ios::scientific);
    tf.precision(8);
    for (const TrajPoint &p : nlp->traj) {
      tf << p.iter;
      for (int j = 0; j < n; j++) tf << "," << p.alpha[j];
      tf << "," << p.value << "," << p.inf_pr << "," << p.inf_du << "\n";
    }
  }

  write_mesh_txt(dir + "/shape_initial.txt", *shape, nlp->init_alpha.data());
  write_mesh_txt(dir + "/shape_final.txt", *shape, nlp->final_alpha.data());

  double size_measure = min_size_constraint.eval(*shape, nlp->final_alpha.data());
  const char *size_label = (shape->dim() == 2) ? "body area" : "body volume";

  // result.txt
  {
    std::ofstream rf((dir + "/result.txt").c_str());
    rf << "experiment    : " << cfg.experiment << "\n";
    rf << "objective     : " << cfg.objective_name << "\n";
    rf << "shape         : " << cfg.shape_name << "\n";
#ifdef SPARTA_AD
    rf << "gradient      : AD (forward-mode, "
       << (cfg.run.score_correction ? "score-corrected" : "uncorrected") << ")\n";
#else
    rf << "gradient      : finite-difference\n";
#endif
    rf << "ipopt status  : " << status_str(st) << "\n";
    rf << "iterations    : " << (nlp->traj.empty() ? 0 : nlp->traj.back().iter) << "\n";
    rf.setf(std::ios::scientific); rf.precision(8);
    rf << "alpha start   : ";
    for (int j = 0; j < n; j++) rf << cfg.alpha[j] << (j + 1 < n ? " " : "\n");
    rf << "alpha final   : ";
    for (int j = 0; j < n; j++) rf << nlp->final_alpha[j] << (j + 1 < n ? " " : "\n");
    rf << "value start   : " << nlp->init_value << "\n";
    rf << "value final   : " << nlp->final_value << "\n";
    rf << "value reduction: " << (nlp->init_value - nlp->final_value) << "\n";
    rf << size_label << "     : " << size_measure << "\n";
    if (cfg.min_area_set) {
      rf << "min_size bound : " << cfg.min_area
         << (size_measure + 1e-9 >= cfg.min_area ? "  (satisfied)" : "  (VIOLATED)")
         << "\n";
    }
  }

  // shapes.svg -- Bezier-specific renderer (see svg_shape.h), so only
  // meaningful for a 2D shape.
  if (shape->dim() == 2) {
    int svg = write_shapes_svg(dir + "/shapes.svg",
                               nlp->init_alpha.data(), nlp->init_value,
                               nlp->final_alpha.data(), nlp->final_value,
                               cfg.chord, cfg.nseg, cfg.objective_name);
    if (svg != 0)
      fprintf(stderr, "WARNING: failed to write shapes.svg\n");
  }

  printf("\ndone: %s\n", status_str(st));
  printf("  value  %.6e  ->  %.6e   (reduction %.3e)\n",
         nlp->init_value, nlp->final_value, nlp->init_value - nlp->final_value);
  printf("  alpha final = ");
  for (int j = 0; j < n; j++) printf("%.5g ", nlp->final_alpha[j]);
  printf("\n  results in  %s/\n", dir.c_str());

  return (st == Solve_Succeeded || st == Solved_To_Acceptable_Level) ? 0 : 1;
}

int main(int narg, char **arg)
{
  try {
    return main_impl(narg, arg);
  } catch (const std::exception &e) {
    fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
}
