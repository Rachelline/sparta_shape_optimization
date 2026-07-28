/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer

   ad_seed.h: the ONE seam an eventual AD-engine swap touches. Every
   Parametrization/Objective/shape_case call site that needs to poke a
   derivative into (or read one out of) an sfloat goes through these two
   functions -- never .fastAccessDx()/.dx() directly.

   ad_seed(x, drow, n):    seed x's derivative slots from drow (length n,
                           a row of a Jacobian d(x)/d(alpha)).
   ad_extract(x, drow, n): read x's derivative slots into drow (length n).

   Both are no-ops (nothing to seed/extract) in the stock build, where
   sfloat == double.

   Landmine note (see docs/PLAN.md's Phase-2 section, and
   tools/ad_verify/sacado_seed_selftest.cpp, which validates this
   empirically against this fork's actual installed Sacado): a
   default-constructed Sacado::Fad::SFad can have size()==0 in general,
   which would make a raw .fastAccessDx(j) poke into an as-yet-
   unconstructed field a silent no-op. That landmine does NOT apply to
   this fork's actual storage policy (SFad/StaticFixedStorage, Trilinos
   17.0.0) -- size() is always the compile-time SPARTA_AD_NDIR,
   confirmed by sacado_seed_selftest.cpp's 11/11 passing checks -- so a
   direct poke after the value is already set (the pattern below) is
   correct here. If a future engine swap reintroduces that landmine,
   this is the only function that needs to change.

   Not yet exercised with real derivatives: this milestone only builds
   the stock (FD-gradient) path. See docs/ad_phase_c_investigation/
   FINDINGS.md for the separate, already-diagnosed reason AD gradients
   aren't trustworthy yet even once this is wired up (the flux-measure
   derivative) -- unrelated to this file, which is pure plumbing.
------------------------------------------------------------------------- */

#ifndef SPARTA_AD_SEED_H
#define SPARTA_AD_SEED_H

#include "sfloat.h"

#ifdef SPARTA_AD

inline void ad_seed(sfloat &x, const double *drow, int n)
{
  for (int j = 0; j < n; j++) x.fastAccessDx(j) = drow[j];
}

inline void ad_extract(const sfloat &x, double *drow, int n)
{
  for (int j = 0; j < n; j++) drow[j] = x.fastAccessDx(j);
}

#else

inline void ad_seed(sfloat & /*x*/, const double * /*drow*/, int /*n*/) {}

inline void ad_extract(const sfloat & /*x*/, double *drow, int n)
{
  for (int j = 0; j < n; j++) drow[j] = 0.0;
}

#endif

#endif
