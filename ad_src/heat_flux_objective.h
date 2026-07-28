/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer

   heat_flux_objective.h: total surface heat flux, summed over the whole
   body. Same shape as DragObjective, swapping `fx` for the already-
   existing `etot` compute-surf keyword -- no new SPARTA-side capability
   needed (confirmed: compute_surf.cpp's ETOT case is unmodified core
   SPARTA). This is the concrete proof docs/PLAN.md's Objective split
   holds for more than one QOI.
------------------------------------------------------------------------- */

#ifndef SPARTA_HEAT_FLUX_OBJECTIVE_H
#define SPARTA_HEAT_FLUX_OBJECTIVE_H

#include "objective.h"

class HeatFluxObjective : public Objective {
 public:
  void setup(void *spa, int navg) const override;
  double extract(void *spa, int ndesign, double *grad) const override;
};

#endif
