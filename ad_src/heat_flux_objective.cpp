/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   heat_flux_objective.cpp -- see heat_flux_objective.h
------------------------------------------------------------------------- */

#include "heat_flux_objective.h"
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

void HeatFluxObjective::setup(void *spa, int navg) const
{
  char line[256];

  cmd(spa, "compute hf_surf surf all all etot");

  std::snprintf(line, sizeof(line),
               "fix hf_favg ave/surf all 1 %d %d c_hf_surf[*]",
               navg, navg);
  cmd(spa, line);

  // see drag_objective.cpp: a single-value fix ave/surf is a per-surf
  // VECTOR (f_ID, no bracket), not an array (f_ID[1]).
  cmd(spa, "compute hf_result reduce sum f_hf_favg");
}

double HeatFluxObjective::extract(void *spa, int ndesign, double *grad) const
{
  sfloat *s = (sfloat *) sparta_extract_compute(spa, (char *) "hf_result", 0, 0);
  if (!s) {
    std::fprintf(stderr, "HeatFluxObjective::extract: compute hf_result "
                         "not found\n");
    std::exit(1);
  }

  if (grad) {
#ifdef SPARTA_AD
    ad_extract(*s, grad, ndesign);
#else
    static bool warned = false;
    if (!warned) {
      std::fprintf(stderr, "HeatFluxObjective::extract: no gradient "
                           "available in the stock build (use grad_fd "
                           "instead)\n");
      warned = true;
    }
    for (int j = 0; j < ndesign; j++) grad[j] = 0.0;
#endif
  }

  return spval(*s);
}
