/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer

   heat_flux_objective.h: total surface heat flux, summed over the whole
   body. A thin named subclass of SurfSumObjective (keyword "etot", tag
   "hf") -- see surf_sum_objective.h for the shared measurement chain
   that this and DragObjective both use. `etot` needs no new SPARTA-side
   capability (compute_surf.cpp's ETOT case is unmodified core SPARTA).
------------------------------------------------------------------------- */

#ifndef SPARTA_HEAT_FLUX_OBJECTIVE_H
#define SPARTA_HEAT_FLUX_OBJECTIVE_H

#include "surf_sum_objective.h"

class HeatFluxObjective : public SurfSumObjective {
 public:
  HeatFluxObjective() : SurfSumObjective("etot", "hf") {}
};

#endif
