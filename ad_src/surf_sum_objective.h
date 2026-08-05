/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer

   surf_sum_objective.h: shared implementation behind DragObjective and
   HeatFluxObjective, which differ only in which compute_surf keyword
   they sum (fx / etot) and the fix/compute id strings they use. Ports
   the reference ad_src's measurement chain unchanged in spirit:
       compute <tag>_surf surf all all <keyword>
       fix <tag>_favg ave/surf all 1 navg navg c_<tag>_surf
       compute <tag>_result reduce sum f_<tag>_favg
   DragObjective/HeatFluxObjective stay as concrete named subclasses
   (default-constructible, same class names as before this refactor) so
   driver code and --objective strings are unaffected.
------------------------------------------------------------------------- */

#ifndef SPARTA_SURF_SUM_OBJECTIVE_H
#define SPARTA_SURF_SUM_OBJECTIVE_H

#include "objective.h"

#include <string>

class SurfSumObjective : public Objective {
 public:
  // keyword: a compute-surf per-surf-element value (e.g. "fx", "etot").
  // tag: short id prefix for this objective's private compute/fix
  // names (e.g. "drag", "hf") -- keeps two SurfSumObjective instances
  // usable on the same SPARTA instance without name collisions.
  SurfSumObjective(std::string keyword, std::string tag);

  void setup(void *spa, int navg) const override;
  double extract(void *spa, int ndesign, double *grad) const override;

 private:
  std::string keyword_;
  std::string tag_;
};

#endif
