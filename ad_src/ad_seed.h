/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer

   ad_seed.h: the AD-engine read seam. Every Objective::extract() call
   site that needs to read a derivative out of an sfloat goes through
   ad_extract() -- never .fastAccessDx()/.dx() directly.

   ad_extract(x, drow, n): read x's derivative slots into drow (length n).
   A no-op (fills zeros) in the stock build, where sfloat == double.

   The seeding side of this seam lives in src/read_surf.cpp's
   SPARTA_AD_SEED_JACFILE hook (writes derivatives into surf points
   read from a Jacobian side-file written by shape_case.cpp) -- there is
   no seeding function here to call; SPARTA reads the seed file itself.

   Exercised with real derivatives: shape_main_ad / shape_opt_ad /
   power_law_opt_ad (see docs/ad_phase_c_investigation/FINDINGS.md for
   the known low bias forward-mode AD carries on DSMC surface tallies,
   and the SPARTA_AD_SCORE_CORRECTION runtime toggle that corrects it
   for fluxscale-weighted tallies).
------------------------------------------------------------------------- */

#ifndef SPARTA_AD_SEED_H
#define SPARTA_AD_SEED_H

#include "sfloat.h"

#ifdef SPARTA_AD

inline void ad_extract(const sfloat &x, double *drow, int n)
{
  for (int j = 0; j < n; j++) drow[j] = x.fastAccessDx(j);
}

#else

inline void ad_extract(const sfloat & /*x*/, double *drow, int n)
{
  for (int j = 0; j < n; j++) drow[j] = 0.0;
}

#endif

#endif
