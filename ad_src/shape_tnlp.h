/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer

   shape_tnlp: an IPOPT TNLP that minimizes obj(alpha) for a given
   Parametrization/Objective pair, subject to box bounds:

       shape.bounds(lo,hi)[j] <= alpha[j] <= shape.bounds(lo,hi)[j]

   plus an optional list of general inequality constraints (Constraint,
   e.g. MinAreaConstraint) -- empty by default, unconstrained behavior
   unchanged when no constraints are passed.

   Generalizes the reference ad_src's DragTNLP (fixed double[4], hardcoded
   Bezier + drag) to a runtime-sized alpha (length shape.ndesign()) and an
   arbitrary Parametrization/Objective pair -- the same generalization
   shape_main.cpp already did for single evaluations, applied to the
   optimizer loop.

   Gradient source: stock build always uses shape_case.h's grad_fd()
   (common-random-number central finite difference). AD build (compiled
   with -DSPARTA_AD, linked against an AD-configured SPARTA library) gets
   the whole gradient from a single SPARTA run via evaluate_avg(...,grad),
   using src/read_surf.cpp's SPARTA_AD_SEED_JACFILE hook to seed all
   ndesign directions at once. IMPORTANT: AD gradients are used exactly as
   the solver returns them -- no flux-measure-derivative correction is
   applied (explicit instruction). They carry a known, well-characterized
   bias (~50% low for specular reflection, ~40% for diffuse -- see
   docs/ad_phase_c_investigation/FINDINGS.md); since that bias is a strict
   positive multiplier, not a sign flip, the *direction* AD reports should
   still usually point downhill, but do not expect the same trajectory,
   iteration count, or convergence behavior as the FD-driven build.

   Each evaluation is expensive (nseeds SPARTA runs for a value; stock
   build: (2*ndesign+1)*nseeds more for a gradient; AD build: nseeds more,
   not (2*ndesign+1)*nseeds, since gradient comes free with the value), so
   results are cached per design point to avoid recomputation between the
   paired eval_f / eval_grad_f calls IPOPT makes at the same x.
------------------------------------------------------------------------- */

#ifndef SPARTA_SHAPE_TNLP_H
#define SPARTA_SHAPE_TNLP_H

#include "IpTNLP.hpp"
#include "constraint.h"
#include "objective.h"
#include "parametrization.h"
#include "shape_case.h"

#include <vector>

// One recorded IPOPT iterate (for the output trajectory + progress bar).
struct TrajPoint {
  int    iter;
  std::vector<double> alpha;
  double value;
  double inf_pr;   // primal infeasibility (0 here: bounds only)
  double inf_du;   // dual infeasibility (optimality measure)
};

class ShapeTNLP : public Ipopt::TNLP {
 public:
  // constraints: empty by default (today's unconstrained behavior,
  // unchanged); each is g_lo <= constraint->eval(...) <= g_hi. Pointers
  // are not owned -- caller (opt_main.cpp) keeps them alive for the
  // TNLP's lifetime.
  ShapeTNLP(const Parametrization &shape, const Objective &obj,
           const double *alpha0, const double *xlo, const double *xhi,
           const int *seeds, int nseeds,
           const RunConfig &c, double h,
           int max_iter, bool show_bar, double obj_scale,
           const std::vector<const Constraint *> &constraints =
             std::vector<const Constraint *>());

  // --- IPOPT TNLP interface -------------------------------------------
  bool get_nlp_info(Ipopt::Index &n, Ipopt::Index &m,
                    Ipopt::Index &nnz_jac_g, Ipopt::Index &nnz_h_lag,
                    IndexStyleEnum &index_style) override;

  bool get_bounds_info(Ipopt::Index n, Ipopt::Number *x_l, Ipopt::Number *x_u,
                       Ipopt::Index m, Ipopt::Number *g_l,
                       Ipopt::Number *g_u) override;

  // Scale the (tiny, ~1e-21) objective to O(1) so it is not swamped by
  // IPOPT's log-barrier term. obj_scaling = 1/|value(start)|.
  bool get_scaling_parameters(Ipopt::Number &obj_scaling,
                              bool &use_x_scaling, Ipopt::Index n,
                              Ipopt::Number *x_scaling, bool &use_g_scaling,
                              Ipopt::Index m, Ipopt::Number *g_scaling) override;

  bool get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number *x,
                          bool init_z, Ipopt::Number *z_L, Ipopt::Number *z_U,
                          Ipopt::Index m, bool init_lambda,
                          Ipopt::Number *lambda) override;

  bool eval_f(Ipopt::Index n, const Ipopt::Number *x, bool new_x,
              Ipopt::Number &obj_value) override;

  bool eval_grad_f(Ipopt::Index n, const Ipopt::Number *x, bool new_x,
                   Ipopt::Number *grad_f) override;

  bool eval_g(Ipopt::Index n, const Ipopt::Number *x, bool new_x,
              Ipopt::Index m, Ipopt::Number *g) override;

  bool eval_jac_g(Ipopt::Index n, const Ipopt::Number *x, bool new_x,
                  Ipopt::Index m, Ipopt::Index nele_jac, Ipopt::Index *iRow,
                  Ipopt::Index *jCol, Ipopt::Number *values) override;

  bool intermediate_callback(Ipopt::AlgorithmMode mode, Ipopt::Index iter,
                             Ipopt::Number obj_value, Ipopt::Number inf_pr,
                             Ipopt::Number inf_du, Ipopt::Number mu,
                             Ipopt::Number d_norm,
                             Ipopt::Number regularization_size,
                             Ipopt::Number alpha_du, Ipopt::Number alpha_pr,
                             Ipopt::Index ls_trials,
                             const Ipopt::IpoptData *ip_data,
                             Ipopt::IpoptCalculatedQuantities *ip_cq) override;

  void finalize_solution(Ipopt::SolverReturn status, Ipopt::Index n,
                         const Ipopt::Number *x, const Ipopt::Number *z_L,
                         const Ipopt::Number *z_U, Ipopt::Index m,
                         const Ipopt::Number *g, const Ipopt::Number *lambda,
                         Ipopt::Number obj_value, const Ipopt::IpoptData *ip_data,
                         Ipopt::IpoptCalculatedQuantities *ip_cq) override;

  // --- results (read after OptimizeTNLP returns) ----------------------
  std::vector<double> init_alpha;
  double init_value;
  std::vector<double> final_alpha;
  double final_value;
  int    solve_status;          // Ipopt::SolverReturn as int
  std::vector<TrajPoint> traj;

 private:
  // Objective + gradient at x, cached. Returns value(x); fills grad (if
  // want_grad; grad may be NULL). Reuses the cache when x is unchanged.
  double evaluate(const double *x, bool want_grad, double *grad);

  const Parametrization &shape_;
  const Objective &obj_;
  int ndesign_;
  std::vector<double> a0_;
  std::vector<double> xl_, xu_;
  std::vector<int> seeds_;
  RunConfig c_;
  double h_;
  int  max_iter_;
  bool show_bar_;
  double obj_scale_;
  std::vector<const Constraint *> constraints_;

  // value/gradient cache keyed on the design point
  bool   fcache_valid_;
  std::vector<double> fcache_x_;
  double fcache_f_;
  bool   gcache_valid_;
  std::vector<double> gcache_x_;
  std::vector<double> gcache_g_;

  bool got_init_;               // captured init_value on first eval
};

#endif
