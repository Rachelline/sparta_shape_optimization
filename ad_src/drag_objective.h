/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer

   drag_objective.h: x-component of net surface force (drag), summed
   over the whole body. Ports the reference ad_src's measurement chain
   (compute forces surf ... fx / fix ave/surf / compute reduce sum)
   unchanged in spirit -- see drag_objective.cpp for the exact commands.
------------------------------------------------------------------------- */

#ifndef SPARTA_DRAG_OBJECTIVE_H
#define SPARTA_DRAG_OBJECTIVE_H

#include "objective.h"

class DragObjective : public Objective {
 public:
  void setup(void *spa, int navg) const override;
  double extract(void *spa, int ndesign, double *grad) const override;
};

#endif
