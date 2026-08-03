/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer

   drag_objective.h: x-component of net surface force (drag), summed
   over the whole body. A thin named subclass of SurfSumObjective
   (keyword "fx", tag "drag") -- see surf_sum_objective.h for the shared
   measurement chain (compute surf .../ fix ave/surf / compute reduce
   sum) that this and HeatFluxObjective both use.
------------------------------------------------------------------------- */

#ifndef SPARTA_DRAG_OBJECTIVE_H
#define SPARTA_DRAG_OBJECTIVE_H

#include "surf_sum_objective.h"

class DragObjective : public SurfSumObjective {
 public:
  DragObjective() : SurfSumObjective("fx", "drag") {}
};

#endif
