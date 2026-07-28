/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   drag_objective.cpp -- see drag_objective.h
------------------------------------------------------------------------- */

#include "drag_objective.h"
#include "ad_seed.h"

#include "library.h"
#include "sfloat.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace {

void cmd(void *spa, const char *str)
{
  char buf[1024];
  std::snprintf(buf, sizeof(buf), "%s", str);
  sparta_command(spa, buf);
}

}  // namespace

void DragObjective::setup(void *spa, int navg) const
{
  char line[256];

  cmd(spa, "compute drag_forces surf all all fx");

  std::snprintf(line, sizeof(line),
               "fix drag_favg ave/surf all 1 %d %d c_drag_forces[*]",
               navg, navg);
  cmd(spa, line);

  // fix ave/surf with a single tracked value (compute drag_forces has
  // only "fx") produces a per-surf VECTOR (fix_ave_surf.cpp: nvalues==1
  // -> size_per_surf_cols=0), addressed as f_ID with no bracket -- an
  // array-style f_ID[1] fails with "does not calculate a per-surf array".
  cmd(spa, "compute drag_result reduce sum f_drag_favg");
}

double DragObjective::extract(void *spa, int ndesign, double *grad) const
{
  sfloat *s = (sfloat *) sparta_extract_compute(spa, (char *) "drag_result", 0, 0);
  if (!s) {
    std::fprintf(stderr, "DragObjective::extract: compute drag_result "
                         "not found\n");
    std::exit(1);
  }

  if (grad) {
#ifdef SPARTA_AD
    ad_extract(*s, grad, ndesign);
#else
    static bool warned = false;
    if (!warned) {
      std::fprintf(stderr, "DragObjective::extract: no gradient available "
                           "in the stock build (use grad_fd instead)\n");
      warned = true;
    }
    for (int j = 0; j < ndesign; j++) grad[j] = 0.0;
#endif
  }

  return spval(*s);
}
