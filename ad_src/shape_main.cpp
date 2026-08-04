/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer

   shape_main.cpp: CLI entry point. Evaluates a shape's Objective across
   seeds, prints per-seed value/mean/stderr, and the grad_fd() gradient.

   Generalizes the reference ad_src's drag_main.cpp: --objective selects
   the Objective, --shape selects the Shape ("bezier" today; the flag
   exists so the seam is visible for a second shape family). Stock
   builds print grad_fd()'s finite-difference gradient; AD builds print
   the gradient that comes free with evaluate_avg() -- see shape_case.h.
------------------------------------------------------------------------- */

#include "shape_case.h"
#include "bezier_shape.h"
#include "cli.h"
#include "drag_objective.h"
#include "heat_flux_objective.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <vector>

namespace {

void usage()
{
  std::fprintf(stderr,
    "usage: shape_main --alpha x1,y1,x2,y2 [options]\n"
    "  --objective drag|heatflux   (default: drag)\n"
    "  --shape bezier              (default: bezier; only value today)\n"
    "  --seeds s1,s2,...           (default: 12345)\n"
    "  --h H                       FD step (default: 0.05)\n"
    "  --nsettle N  --navg N  --nseg N  --chord L  --vstream V\n"
    "  --wall-temp T  --wall-accom A   (diffuse wall; A in [0,1], default\n"
    "                                   0.0 matches reference_gradient.txt;\n"
    "                                   HeatFluxObjective needs A > 0, e.g.\n"
    "                                   --wall-accom 1.0, or its signal is\n"
    "                                   nearly degenerate -- see shape_case.h)\n"
    "  --specular  --nocoll  --verbose\n"
    "  --score-correction / --no-score-correction   (default: off)\n"
    "                       AD build only: enable/disable the flux-measure\n"
    "                       score-function correction -- see docs/AD_GRADIENTS.md.\n"
    "                       No effect on stock builds (warns if passed there).\n"
    "                       Raw escape hatch: SPARTA_AD_SCORE_CORRECTION=1 env var.\n");
}

}  // namespace

int main_impl(int argc, char **argv)
{
  std::vector<double> alpha_v;
  std::vector<int> seeds_v;
  const char *objective_name = "drag";
  const char *shape_name = "bezier";
  double h = 0.05;
  double chord = 4.0;
  int nseg = 25;
  RunConfig c;

  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "--alpha") && i + 1 < argc) alpha_v = parse_doubles(argv[++i]);
    else if (!strcmp(argv[i], "--seeds") && i + 1 < argc) seeds_v = parse_ints(argv[++i]);
    else if (!strcmp(argv[i], "--objective") && i + 1 < argc) objective_name = argv[++i];
    else if (!strcmp(argv[i], "--shape") && i + 1 < argc) shape_name = argv[++i];
    else if (!strcmp(argv[i], "--h") && i + 1 < argc) h = atof(argv[++i]);
    else if (!strcmp(argv[i], "--nsettle") && i + 1 < argc) c.nsettle = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--navg") && i + 1 < argc) c.navg = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--nseg") && i + 1 < argc) nseg = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--chord") && i + 1 < argc) chord = atof(argv[++i]);
    else if (!strcmp(argv[i], "--vstream") && i + 1 < argc) c.vstream = atof(argv[++i]);
    else if (!strcmp(argv[i], "--wall-temp") && i + 1 < argc) c.wall_temp = atof(argv[++i]);
    else if (!strcmp(argv[i], "--wall-accom") && i + 1 < argc) c.wall_accom = atof(argv[++i]);
    else if (!strcmp(argv[i], "--specular")) c.specular = 1;
    else if (!strcmp(argv[i], "--score-correction")) c.corrections.flux_measure = true;
    else if (!strcmp(argv[i], "--no-score-correction")) c.corrections.flux_measure = false;
    else if (!strcmp(argv[i], "--nocoll")) c.collisions = 0;
    else if (!strcmp(argv[i], "--verbose")) c.verbose = 1;
    else { usage(); return 1; }
  }

  if (alpha_v.empty()) { usage(); return 1; }
  if (seeds_v.empty()) seeds_v.push_back(12345);
  c.corrections.warn_if_noop();

  BezierShape bezier_shape(chord, nseg);
  Shape *shape = 0;
  if (!strcmp(shape_name, "bezier")) shape = &bezier_shape;
  else {
    std::fprintf(stderr, "unknown --shape '%s'\n", shape_name);
    return 1;
  }

  if ((int) alpha_v.size() != shape->ndesign()) {
    std::fprintf(stderr, "--alpha needs %d values for shape '%s', got %zu\n",
                shape->ndesign(), shape_name, alpha_v.size());
    return 1;
  }

  DragObjective drag_obj;
  HeatFluxObjective hf_obj;
  Objective *obj = 0;
  if (!strcmp(objective_name, "drag")) obj = &drag_obj;
  else if (!strcmp(objective_name, "heatflux")) obj = &hf_obj;
  else {
    std::fprintf(stderr, "unknown --objective '%s'\n", objective_name);
    return 1;
  }

  int ndesign = shape->ndesign();
  int nseeds = (int) seeds_v.size();

  std::printf("objective=%s shape=%s alpha=", objective_name, shape_name);
  for (int j = 0; j < ndesign; j++) std::printf("%s%.6g", j ? "," : "", alpha_v[j]);
  std::printf("\n");
  std::printf("score_correction=%s\n", c.corrections.flux_measure ? "on" : "off");

  std::vector<double> vals(nseeds);
  for (int k = 0; k < nseeds; k++) {
    vals[k] = evaluate(*shape, *obj, alpha_v.data(), seeds_v[k], c);
    std::printf("  seed=%d value=%.10g\n", seeds_v[k], vals[k]);
  }

  double mean = 0.0;
  for (double v : vals) mean += v;
  mean /= nseeds;

  double stderr_val = 0.0;
  if (nseeds > 1) {
    double var = 0.0;
    for (double v : vals) var += (v - mean) * (v - mean);
    var /= (nseeds - 1);
    stderr_val = std::sqrt(var / nseeds);
  }
  std::printf("mean=%.10g", mean);
  if (nseeds > 1) std::printf("  stderr=%.3g", stderr_val);
  std::printf("\n");

  std::vector<double> grad(ndesign);
#ifdef SPARTA_AD
  // One SPARTA run yields value AND gradient together (see shape_case.cpp /
  // read_surf.cpp's SPARTA_AD_SEED_JACFILE hook). See docs/AD_GRADIENTS.md
  // for the known low bias this gradient carries and what
  // --score-correction does and doesn't fix.
  double base = evaluate_avg(*shape, *obj, alpha_v.data(), seeds_v.data(),
                             nseeds, c, grad.data());
  std::printf("grad_ad (%s, see docs/AD_GRADIENTS.md): [",
             c.corrections.flux_measure ? "score-corrected" : "uncorrected");
#else
  double base = grad_fd(*shape, *obj, alpha_v.data(), seeds_v.data(), nseeds,
                        c, h, grad.data());
  std::printf("grad_fd (h=%.4g): [", h);
#endif
  for (int j = 0; j < ndesign; j++) std::printf("%s%.6g", j ? ", " : "", grad[j]);
  std::printf("]  (base value used: %.10g)\n", base);

  return 0;
}

int main(int argc, char **argv)
{
  try {
    return main_impl(argc, argv);
  } catch (const std::exception &e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
}
