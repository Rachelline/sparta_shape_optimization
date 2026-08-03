/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer

   shape_case.h: the shared runner. Replaces the reference ad_src's
   drag_objective.h/.cpp runner half -- RunConfig was DragCase, minus
   anything drag-specific; evaluate()/evaluate_avg()/grad_fd() are its
   drag()/drag_avg()/grad_fd(), generalized to take a Parametrization and
   an Objective instead of hardcoding Bezier + drag.

   Deviations from docs/PLAN.md's terse RunConfig/evaluate sketch, and
   why: RunConfig here adds `species_names` (this repo's own single-
   species tools/ad_verify/N.species convention, vs. the reference's
   hardcoded 2-species "air") and `tinf` (mixture's freestream
   temperature -- the reference's deck never sets it explicitly and
   relies on a global fallback; setting it here explicitly avoids that
   ambiguity, and HeatFluxObjective's etot needs a sane thermal energy
   scale to be meaningful regardless).
------------------------------------------------------------------------- */

#ifndef SPARTA_SHAPE_CASE_H
#define SPARTA_SHAPE_CASE_H

#include "objective.h"
#include "parametrization.h"

struct RunConfig {
  double chord = 4.0;
  int nseg = 25;                    // meaning is shape-family-defined

  double origin[2] = {3.0, 5.0};
  double boxhi[2]  = {10.0, 10.0};
  int    grid[2]   = {20, 20};

  double vstream = 100.0;
  double tinf = 300.0;
  double nrho = 1.0, fnum = 0.001, tstep = 0.0001;
  const char *species_file = "N.species";
  const char *vss_file     = "N.vss";
  const char *species_names = "N";  // space-separated, as in `species` cmd

  int nsettle = 500;
  int navg    = 500;

  int verbose  = 0;
  int progress = 0;

  int specular   = 0;    // 1 = specular wall, 0 = diffuse
  double wall_temp  = 300.0;
  double wall_accom = 0.0;  // diffuse-wall accommodation coeff, [0,1].
                            // Default matches the reference ad_src's
                            // DragCase (hardcoded 0.0 there) so DragObjective
                            // reproduces reference_gradient.txt out of the
                            // box. IMPORTANT for HeatFluxObjective: acc=0
                            // means the wall does NOT thermally accommodate
                            // the gas, which makes total heat flux nearly
                            // degenerate by construction (barely any energy
                            // exchange to measure) -- callers driving
                            // HeatFluxObjective should override this to
                            // something meaningfully nonzero (1.0 = full
                            // accommodation was used for this milestone's
                            // own HeatFlux verification).
  int collisions = 1;   // 0 = free-molecular (no `collide vss`)

  bool score_correction = false;  // AD build only: sets
                                   // SPARTA_AD_SEED_JACFILE-sibling env var
                                   // SPARTA_AD_SCORE_CORRECTION before each
                                   // SPARTA invocation, enabling the
                                   // flux-measure score-function correction
                                   // in compute_surf.cpp (FINDINGS.md,
                                   // FINDING 2). No effect on stock builds
                                   // or on tallied VALUES; changes AD
                                   // gradients only. Off by default.
};

// Single realization: builds the shape, writes a tmp surf file, runs
// SPARTA, returns the Objective's value. grad (length shape.ndesign()),
// if non-NULL, is filled by the AD build only -- see objective.h.
double evaluate(const Parametrization &shape, const Objective &obj,
               const double *alpha, int seed, const RunConfig &c,
               double *grad = 0);

// Plain mean of evaluate() over seeds; grad (if requested) is the
// seed-averaged gradient.
double evaluate_avg(const Parametrization &shape, const Objective &obj,
                    const double *alpha, const int *seeds, int nseeds,
                    const RunConfig &c, double *grad = 0);

// Common-random-numbers central finite difference: for each of
// shape.ndesign() components, two evaluate() calls at alpha +/- h
// (same seed per pair, for noise correlation), averaged over `nseeds`.
// Returns the unperturbed evaluate_avg() value; fills grad (length
// shape.ndesign()).
double grad_fd(const Parametrization &shape, const Objective &obj,
              const double *alpha, const int *seeds, int nseeds,
              const RunConfig &c, double h, double *grad);

#endif
