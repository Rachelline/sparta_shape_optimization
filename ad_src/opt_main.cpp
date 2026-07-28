/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer

   opt_main: IPOPT-driven shape-optimization driver. Generalizes the
   reference ad_src's opt_ipopt.cpp: --objective drag|heatflux --shape
   bezier selection (mirrors shape_main.cpp exactly), gradients always
   from grad_fd() (see shape_tnlp.h for why there's no AD branch this
   milestone).

   Bounds come from Parametrization::bounds() (built in Milestone 1), not
   CLI flags -- the reference's --x-lo/--x-hi/--y-lo/--y-hi are
   Bezier-specific naming baked into what's meant to be a generic driver.

   Results are written to a dated output folder:
     output/fd_<YYYY-MM-DD>_<experiment>[_N]/
       config.txt      resolved parameters
       ipopt.log       IPOPT's own iteration log
       trajectory.csv  per-iteration alpha / value / infeasibility
       result.txt      start & final alpha, value, status, area
       shapes.svg      initial (faint) vs optimized (bold) outline
                        -- only written for --shape bezier

   IMPORTANT for --objective heatflux: total heat flux needs a
   thermally-accommodating wall (--wall-accom > 0) and enough averaging
   to see past its noise floor (energy is a higher/quadratic velocity
   moment than momentum, so it's inherently noisier than drag at the
   same sample size -- see docs/PLAN.md's Milestone-1 section). This
   driver does not auto-inject those; pass them explicitly, e.g.
   --wall-accom 1.0 --navg 4000.

   Run from ad_src/ (where N.species/N.vss live):
     ./build/shape_opt --alpha 1.3,1.0,2.7,0.8 --experiment baseline
------------------------------------------------------------------------- */

#include "bezier_geom.h"
#include "bezier_parametrization.h"
#include "drag_objective.h"
#include "heat_flux_objective.h"
#include "shape_case.h"
#include "shape_tnlp.h"
#include "svg_shape.h"

#include "IpIpoptApplication.hpp"

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <vector>

using namespace Ipopt;

static const char *USAGE =
  "usage: shape_opt [--input FILE] [options]\n"
  "  --objective drag|heatflux   (default: drag)\n"
  "  --shape bezier              (default: bezier; only value today)\n"
  "  --alpha x1,y1,x2,y2         (required unless set in --input)\n"
  "  --seeds s1,s2,...           (default: 12345,67890,13579 -- more\n"
  "                               than shape_main's single-seed default,\n"
  "                               since a real optimization loop needs a\n"
  "                               less noisy gradient than a one-off eval)\n"
  "  --max-iter N  --tol T  --acceptable-tol T  --acceptable-iter N\n"
  "  --h H                       FD step (default: 0.05)\n"
  "  --experiment NAME\n"
  "  --nsettle N  --navg N  --nseg N  --chord L  --vstream V\n"
  "  --wall-temp T  --wall-accom A   (diffuse wall; see note above for\n"
  "                                   --objective heatflux)\n"
  "  --specular  --nocoll  --verbose\n"
  "\n"
  "  Minimizes the chosen objective with IPOPT, subject to the shape's\n"
  "  own box bounds (Parametrization::bounds()). Run from a dir with\n"
  "  N.species/N.vss. Output written to output/fd_<date>_<experiment>/.\n";

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
  RunConfig c;
  bool alpha_set = false;
};

// ---- small string helpers -----------------------------------------------

static std::string trim(const std::string &s)
{
  size_t a = 0, b = s.size();
  while (a < b && isspace((unsigned char) s[a])) a++;
  while (b > a && isspace((unsigned char) s[b - 1])) b--;
  return s.substr(a, b - a);
}

static int parse_doubles(const std::string &v, std::vector<double> &out)
{
  out.clear();
  std::string tmp = v;
  for (char &ch : tmp) if (ch == ',') ch = ' ';
  std::istringstream ss(tmp);
  double d;
  while (ss >> d) out.push_back(d);
  return (int) out.size();
}

static void parse_ints(const std::string &v, std::vector<int> &out)
{
  out.clear();
  std::string tmp = v;
  for (char &ch : tmp) if (ch == ',') ch = ' ';
  std::istringstream ss(tmp);
  int d;
  while (ss >> d) out.push_back(d);
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
      if (parse_doubles(val, cfg.alpha) > 0) cfg.alpha_set = true;
    } else if (key == "seeds") { parse_ints(val, cfg.seeds);
    } else if (key == "max_iter") { cfg.max_iter = atoi(val.c_str());
    } else if (key == "tol") { cfg.tol = atof(val.c_str());
    } else if (key == "acceptable_tol") { cfg.acc_tol = atof(val.c_str());
    } else if (key == "acceptable_iter") { cfg.acc_iter = atoi(val.c_str());
    } else if (key == "h") { cfg.h = atof(val.c_str());
    } else if (key == "experiment") { cfg.experiment = val;
    } else if (key == "nseg") { cfg.c.nseg = atoi(val.c_str());
    } else if (key == "chord") { cfg.c.chord = atof(val.c_str());
    } else if (key == "nsettle") { cfg.c.nsettle = atoi(val.c_str());
    } else if (key == "navg") { cfg.c.navg = atoi(val.c_str());
    } else if (key == "vstream") { cfg.c.vstream = atof(val.c_str());
    } else if (key == "wall_temp") { cfg.c.wall_temp = atof(val.c_str());
    } else if (key == "wall_accom") { cfg.c.wall_accom = atof(val.c_str());
    } else if (key == "specular") { cfg.c.specular = atoi(val.c_str());
    } else if (key == "collisions") { cfg.c.collisions = atoi(val.c_str());
    } else if (key == "verbose") { cfg.c.verbose = atoi(val.c_str());
    }
    // unknown keys ignored
  }
  return true;
}

// ---- filesystem ----------------------------------------------------------

static bool dir_exists(const std::string &p)
{
  struct stat st;
  return stat(p.c_str(), &st) == 0 && S_ISDIR(st.st_mode);
}

// Create output/fd_<date>_<experiment>[_N]/ and return its path (empty on
// fail). "fd_" prefix (not "ad_") so an eventual AD-gradient variant's
// runs land in a distinguishable, non-colliding place later.
static std::string make_output_dir(const std::string &experiment)
{
  mkdir("output", 0755);   // ok if it already exists
  char datestr[32];
  time_t t = time(0);
  strftime(datestr, sizeof(datestr), "%Y-%m-%d", localtime(&t));
  std::string base = std::string("output/fd_") + datestr + "_" + experiment;
  std::string dir = base;
  int suffix = 2;
  while (dir_exists(dir)) {
    dir = base + "_" + std::to_string(suffix++);
  }
  if (mkdir(dir.c_str(), 0755) != 0 && errno != EEXIST) {
    fprintf(stderr, "ERROR: cannot create output dir '%s': %s\n",
            dir.c_str(), strerror(errno));
    return std::string();
  }
  return dir;
}

// ---- geometry: signed area of the final body (bezier only) --------------

static double body_area(const double alpha[4], double chord, int nseg)
{
  int npt = 2 * nseg + 1;
  std::vector<double> pts(2 * npt), norms(2 * (2 * nseg));
  BezierGeom::symmetric_body_to_lines(alpha, chord, nseg, pts.data(),
                                      norms.data());
  return BezierGeom::signed_area(2 * nseg, pts.data());
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

int main(int narg, char **arg)
{
  Config cfg;

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
      if (parse_doubles(arg[++i], cfg.alpha) > 0) cfg.alpha_set = true;
    } else if (!strcmp(arg[i], "--seeds") && i + 1 < narg) {
      parse_ints(arg[++i], cfg.seeds);
    } else if (!strcmp(arg[i], "--max-iter") && i + 1 < narg) { cfg.max_iter = atoi(arg[++i]);
    } else if (!strcmp(arg[i], "--tol") && i + 1 < narg) { cfg.tol = atof(arg[++i]);
    } else if (!strcmp(arg[i], "--acceptable-tol") && i + 1 < narg) { cfg.acc_tol = atof(arg[++i]);
    } else if (!strcmp(arg[i], "--acceptable-iter") && i + 1 < narg) { cfg.acc_iter = atoi(arg[++i]);
    } else if (!strcmp(arg[i], "--h") && i + 1 < narg) { cfg.h = atof(arg[++i]);
    } else if (!strcmp(arg[i], "--experiment") && i + 1 < narg) { cfg.experiment = arg[++i];
    } else if (!strcmp(arg[i], "--nseg") && i + 1 < narg) { cfg.c.nseg = atoi(arg[++i]);
    } else if (!strcmp(arg[i], "--chord") && i + 1 < narg) { cfg.c.chord = atof(arg[++i]);
    } else if (!strcmp(arg[i], "--nsettle") && i + 1 < narg) { cfg.c.nsettle = atoi(arg[++i]);
    } else if (!strcmp(arg[i], "--navg") && i + 1 < narg) { cfg.c.navg = atoi(arg[++i]);
    } else if (!strcmp(arg[i], "--vstream") && i + 1 < narg) { cfg.c.vstream = atof(arg[++i]);
    } else if (!strcmp(arg[i], "--wall-temp") && i + 1 < narg) { cfg.c.wall_temp = atof(arg[++i]);
    } else if (!strcmp(arg[i], "--wall-accom") && i + 1 < narg) { cfg.c.wall_accom = atof(arg[++i]);
    } else if (!strcmp(arg[i], "--specular")) { cfg.c.specular = 1;
    } else if (!strcmp(arg[i], "--nocoll")) { cfg.c.collisions = 0;
    } else if (!strcmp(arg[i], "--verbose")) { cfg.c.verbose = 1;
    } else if (!strcmp(arg[i], "--help") || !strcmp(arg[i], "-h")) { printf("%s", USAGE); return 0;
    } else { fprintf(stderr, "unknown/incomplete arg: %s\n%s", arg[i], USAGE); return 1; }
  }

  if (!cfg.alpha_set) { fprintf(stderr, "ERROR: --alpha required\n%s", USAGE); return 1; }
  if (cfg.seeds.empty()) cfg.seeds = {12345, 67890, 13579};

  BezierParametrization bezier_shape;
  Parametrization *shape = 0;
  if (cfg.shape_name == "bezier") shape = &bezier_shape;
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

  if (cfg.objective_name == "heatflux" && cfg.c.wall_accom <= 0.0) {
    fprintf(stderr,
            "WARNING: --objective heatflux with wall_accom=%.3g (default 0)"
            " -- total heat flux is nearly degenerate without wall thermal"
            " accommodation. Pass --wall-accom 1.0 (or similar) for a"
            " usable gradient. Proceeding anyway.\n", cfg.c.wall_accom);
  }

  int n = shape->ndesign();
  const bool show_bar = (!cfg.c.verbose) && isatty(fileno(stderr));

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

  // Objective scaling: the objective is ~1e-18 to ~1e-21, so without
  // scaling IPOPT's log-barrier term swamps it and the iterate drifts to
  // the center of the box instead of minimizing the objective. Scale by
  // 1/|value(start)| so the scaled objective is O(1). One evaluation at
  // the start point.
  double f0 = evaluate_avg(*shape, *obj, cfg.alpha.data(), cfg.seeds.data(),
                           (int) cfg.seeds.size(), cfg.c);
  double obj_scale = 1.0 / std::max(fabs(f0), 1e-300);

  // output folder + config echo
  std::string dir = make_output_dir(cfg.experiment);
  if (dir.empty()) return 1;

  char timestamp[64];
  { time_t t = time(0); strftime(timestamp, sizeof(timestamp),
                                 "%Y-%m-%d %H:%M:%S", localtime(&t)); }

  {
    std::ofstream cf((dir + "/config.txt").c_str());
    cf << "timestamp    = " << timestamp << "\n";
    cf << "experiment   = " << cfg.experiment << "\n";
    cf << "objective    = " << cfg.objective_name << "\n";
    cf << "shape        = " << cfg.shape_name << "\n";
    cf << "gradient     = finite-difference (CRN)\n";
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
    cf << "chord        = " << cfg.c.chord << "\n";
    cf << "nseg         = " << cfg.c.nseg << "\n";
    cf << "nsettle      = " << cfg.c.nsettle << "\n";
    cf << "navg         = " << cfg.c.navg << "\n";
    cf << "vstream      = " << cfg.c.vstream << "\n";
    cf << "wall_temp    = " << cfg.c.wall_temp << "\n";
    cf << "wall_accom   = " << cfg.c.wall_accom << "\n";
    cf << "specular     = " << cfg.c.specular << "\n";
    cf << "collisions   = " << cfg.c.collisions << "\n";
  }

  printf("%s optimization (finite-difference gradient)\n",
         cfg.objective_name.c_str());
  printf("  start alpha = ");
  for (int j = 0; j < n; j++) printf("%.5g ", cfg.alpha[j]);
  printf("\n  output dir  = %s\n", dir.c_str());
  fflush(stdout);

  SmartPtr<ShapeTNLP> nlp = new ShapeTNLP(*shape, *obj, cfg.alpha.data(),
                                          xl.data(), xu.data(),
                                          cfg.seeds.data(), (int) cfg.seeds.size(),
                                          cfg.c, cfg.h, cfg.max_iter, show_bar,
                                          obj_scale);

  SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
  app->Options()->SetStringValue("hessian_approximation", "limited-memory");
  app->Options()->SetStringValue("mu_strategy", "adaptive");
  app->Options()->SetStringValue("nlp_scaling_method", "user-scaling");
  app->Options()->SetNumericValue("tol", cfg.tol);
  app->Options()->SetNumericValue("acceptable_tol", cfg.acc_tol);
  app->Options()->SetIntegerValue("acceptable_iter", cfg.acc_iter);
  app->Options()->SetIntegerValue("max_iter", cfg.max_iter);
  app->Options()->SetStringValue("sb", "yes");                  // no banner
  app->Options()->SetIntegerValue("print_level", cfg.c.verbose ? 5 : 0);
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

  double area = -1.0;
  if (cfg.shape_name == "bezier" && n == 4)
    area = body_area(nlp->final_alpha.data(), cfg.c.chord, cfg.c.nseg);

  // result.txt
  {
    std::ofstream rf((dir + "/result.txt").c_str());
    rf << "experiment    : " << cfg.experiment << "\n";
    rf << "objective     : " << cfg.objective_name << "\n";
    rf << "shape         : " << cfg.shape_name << "\n";
    rf << "gradient      : finite-difference\n";
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
    if (area >= 0.0) rf << "body area      : " << area << "\n";
  }

  // shapes.svg -- Bezier-specific (see svg_shape.h)
  if (cfg.shape_name == "bezier" && n == 4) {
    int svg = write_shapes_svg(dir + "/shapes.svg",
                               nlp->init_alpha.data(), nlp->init_value,
                               nlp->final_alpha.data(), nlp->final_value,
                               cfg.c.chord, cfg.c.nseg, cfg.objective_name);
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
