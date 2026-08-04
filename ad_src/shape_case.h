/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer

   shape_case.h: the shared runner. evaluate()/evaluate_avg()/grad_fd()
   are dimension-aware -- one code path drives both a 2D Bezier body and
   a 3D power-law body of revolution, since the underlying SPARTA deck
   differs only in `dimension`, `boundary`, `create_box`, `create_grid`,
   and the surf/jac file element kind (lines vs triangles); everything
   else (species/mixture/emit/collide/run) was already identical.
   Replaces the previous 2D-only evaluate() plus power_law_main.cpp's
   own ~100-line duplicate of the same function.

   RunConfig's origin/boxlo/boxhi/grid are 3-element so one struct covers
   both dimensionalities: for a 2D shape only indices [0,1] (and z as a
   thin symmetric slab, boxlo[2]/boxhi[2]) matter; a 3D shape uses all
   three. The box-bounds check is uniform: for each dim d, translate by
   origin[d], then require boxlo[d] < pts[d] < boxhi[d] -- this covers
   both the 2D case (translate into a corner, [0,boxhi] per axis) and the
   3D case (translate only the apex in x, [boxlo,boxhi] symmetric about
   0 in y/z) as the same rule with different constants, not a branch.
------------------------------------------------------------------------- */

#ifndef SPARTA_SHAPE_CASE_H
#define SPARTA_SHAPE_CASE_H

#include "objective.h"
#include "shape.h"

// AD gradients carry a known low bias on DSMC surface tallies (see
// docs/AD_GRADIENTS.md). This struct is the single on/off switch for
// the runtime correction, factored out of RunConfig so a driver can
// echo/toggle it uniformly without reaching into RunConfig's other
// fields, and so the deferred nfactor_inverse correction (force/torque
// tallies -- not yet covered) has an obvious place to add a second flag
// later without changing every call site that constructs a RunConfig.
struct AdCorrections {
  // Corrects the fluxscale-weighted tallies (press/shear/heat-flux) --
  // sets/unsets SPARTA_AD_SCORE_CORRECTION in compute_surf.cpp. Does
  // NOT cover nfactor_inverse-weighted tallies (force/torque, e.g.
  // DragObjective) -- see docs/AD_GRADIENTS.md's coverage table. No
  // effect on stock (FD) builds or on tallied VALUES; changes AD
  // gradients only. Off by default.
  bool flux_measure = false;

  // Sets/unsets the env var explicitly (never just "set when true") so
  // a process calling evaluate() with different RunConfigs across calls
  // can't inherit a stale value from an earlier one. Safe to call in a
  // stock (non-AD) build too: compute_surf.cpp only reads this env var
  // from its SPARTA_AD-compiled code path, so it's inert there, not
  // wrong -- callers don't need to guard the call with #ifdef SPARTA_AD.
  void apply_to_env() const;

  // Stderr warning if the correction was requested in a build where it
  // can't do anything (today, any stock/FD build -- AD_BUILD reflects
  // whether this translation unit was compiled with -DSPARTA_AD).
  // Purely a UX nicety: callers decide when to invoke it (typically
  // right after CLI parsing), it doesn't gate or change behavior.
  void warn_if_noop() const;
};

struct RunConfig {
  double origin[3] = {3.0, 5.0, 0.0};
  double boxlo[3]  = {0.0, 0.0, -0.5};
  double boxhi[3]  = {10.0, 10.0, 0.5};
  int    grid[3]   = {20, 20, 1};
  const char *boundary = "o r p";

  double vstream = 100.0;
  double tinf = 300.0;
  double nrho = 1.0, fnum = 0.001, tstep = 0.0001;
  const char *species_file = "N.species";
  const char *vss_file     = "N.vss";
  const char *species_names = "N";  // space-separated, as in `species` cmd

  int nsettle = 500;
  int navg    = 500;
  int stats_every = 100;   // cadence for `stats N` when verbose

  int verbose  = 0;

  int specular   = 0;    // 1 = specular wall, 0 = diffuse
  double wall_temp  = 300.0;
  double wall_accom = 0.0;  // diffuse-wall accommodation coeff, [0,1].
                            // IMPORTANT for HeatFluxObjective: acc=0 means
                            // the wall does NOT thermally accommodate the
                            // gas, which makes total heat flux nearly
                            // degenerate by construction -- callers
                            // driving HeatFluxObjective should override
                            // this to something meaningfully nonzero
                            // (1.0 = full accommodation).
  int collisions = 1;   // 0 = free-molecular (no `collide vss`)

  AdCorrections corrections;
};

// Single realization: builds the shape, writes a tmp surf file, runs
// SPARTA, returns the Objective's value. grad (length shape.ndesign()),
// if non-NULL, is filled by the AD build only -- see objective.h.
double evaluate(const Shape &shape, const Objective &obj,
               const double *alpha, int seed, const RunConfig &c,
               double *grad = 0);

// Plain mean of evaluate() over seeds; grad (if requested) is the
// seed-averaged gradient.
double evaluate_avg(const Shape &shape, const Objective &obj,
                    const double *alpha, const int *seeds, int nseeds,
                    const RunConfig &c, double *grad = 0);

// Common-random-numbers central finite difference: for each of
// shape.ndesign() components, two evaluate() calls at alpha +/- h
// (same seed per pair, for noise correlation), averaged over `nseeds`.
// Returns the unperturbed evaluate_avg() value; fills grad (length
// shape.ndesign()).
double grad_fd(const Shape &shape, const Objective &obj,
              const double *alpha, const int *seeds, int nseeds,
              const RunConfig &c, double h, double *grad);

#endif
