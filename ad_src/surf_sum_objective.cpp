/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   surf_sum_objective.cpp -- see surf_sum_objective.h
------------------------------------------------------------------------- */

#include "surf_sum_objective.h"
#include "ad_seed.h"
#include "sparta_util.h"

#include "sfloat.h"

#include <cstdio>
#include <utility>

SurfSumObjective::SurfSumObjective(std::string keyword, std::string tag)
  : keyword_(std::move(keyword)), tag_(std::move(tag))
{
}

void SurfSumObjective::setup(void *spa, int navg) const
{
  char line[256];

  std::snprintf(line, sizeof(line), "compute %s_surf surf all all %s",
               tag_.c_str(), keyword_.c_str());
  cmd(spa, line);

  std::snprintf(line, sizeof(line),
               "fix %s_favg ave/surf all 1 %d %d c_%s_surf[*]",
               tag_.c_str(), navg, navg, tag_.c_str());
  cmd(spa, line);

  // fix ave/surf with a single tracked value (compute <tag>_surf has
  // only one keyword) produces a per-surf VECTOR (fix_ave_surf.cpp:
  // nvalues==1 -> size_per_surf_cols=0), addressed as f_ID with no
  // bracket -- an array-style f_ID[1] fails with "does not calculate a
  // per-surf array".
  std::snprintf(line, sizeof(line), "compute %s_result reduce sum f_%s_favg",
               tag_.c_str(), tag_.c_str());
  cmd(spa, line);
}

double SurfSumObjective::extract(void *spa, int ndesign, double *grad) const
{
  std::string result_name = tag_ + "_result";
  sfloat *s = (sfloat *) extract_compute(spa, result_name.c_str(), 0, 0);
  if (!s)
    die(tag_ + "_objective::extract: compute " + result_name + " not found");

  if (grad) {
#ifdef SPARTA_AD
    ad_extract(*s, grad, ndesign);
#else
    static bool warned = false;
    if (!warned) {
      std::fprintf(stderr, "%s_objective::extract: no gradient available "
                          "in the stock build (use grad_fd instead)\n",
                  tag_.c_str());
      warned = true;
    }
    for (int j = 0; j < ndesign; j++) grad[j] = 0.0;
#endif
  }

  return spval(*s);
}
