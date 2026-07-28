/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer

   objective.h: the Objective interface. A quantity of interest (QOI)
   computed from a converged DSMC run on a given surface.

   Implementations issue whatever compute/fix commands they need in
   setup() -- called once surf + flow are set up (right after
   read_surf/surf_modify), before the settle/averaging run window -- then
   read the result (and, in the AD build, its gradient) back out in
   extract(), called after the averaging run has completed.

   Each Objective owns its own compute/fix ids internally (private
   literal strings): shape_case.cpp never hardcodes a compute name, and
   since only one Objective drives any single SPARTA instance at a time,
   name collisions across Objectives are a non-issue.

   Deviation from docs/PLAN.md's terse interface sketch: setup() there
   was shown as setup(void *spa) with no other arguments. In practice an
   Objective needs `navg` to build its own `fix ave/surf all 1 navg navg
   ...` averaging window -- that's the one piece of RunConfig every
   Objective structurally needs, so it's threaded through explicitly
   rather than the Objective reaching into a RunConfig it otherwise has
   no business knowing about.
------------------------------------------------------------------------- */

#ifndef SPARTA_OBJECTIVE_H
#define SPARTA_OBJECTIVE_H

class Objective {
 public:
  virtual ~Objective() {}

  // spa: opaque SPARTA library handle (void* from sparta_open_no_mpi),
  // passed straight to sparta_command(). navg: length of the averaging
  // window this Objective's own `fix ave/surf` should use.
  virtual void setup(void *spa, int navg) const = 0;

  // After the averaging window has run: return the QOI's value.
  // ndesign: length of grad (may be NULL if the caller doesn't want a
  // gradient -- every call site in this milestone, since only FD
  // gradients are wired up). In the AD build, fills
  // grad[j] = d(QOI)/d(alpha[j]) via ad_extract(); in the stock build,
  // fills grad with 0.0 (once, with a one-time warning) since no
  // gradient is available that way.
  virtual double extract(void *spa, int ndesign, double *grad) const = 0;
};

#endif
