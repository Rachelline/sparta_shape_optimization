/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   shape_tnlp: IPOPT TNLP implementation. See shape_tnlp.h for the design.
------------------------------------------------------------------------- */

#include "shape_tnlp.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <unistd.h>   // isatty

using namespace Ipopt;

// ---- colorful progress bar (rendered on stderr) -----------------------
// Pure terminal-escape-code logic, no Parametrization/Objective
// dependency at all -- ported verbatim from the reference's drag_tnlp.cpp.

static int bar_color(double t)
{
  static const int ramp[] = {
    51, 50, 49, 48, 47, 46, 82, 118, 154, 190, 226, 220, 214, 208, 202, 196
  };
  const int n = (int)(sizeof(ramp) / sizeof(ramp[0]));
  int i = (int)(t * (n - 1) + 0.5);
  if (i < 0) i = 0;
  if (i >= n) i = n - 1;
  return ramp[i];
}

static void draw_bar(int iter, int max_iter, double value, double inf_du)
{
  const int width = 30;
  double frac = (max_iter > 0) ? (double)iter / (double)max_iter : 0.0;
  if (frac > 1.0) frac = 1.0;
  int fill = (int)(frac * width + 0.5);

  fprintf(stderr, "\r\033[1m[\033[0m");
  for (int i = 0; i < width; i++) {
    if (i < fill) {
      double t = (width > 1) ? (double)i / (double)(width - 1) : 0.0;
      fprintf(stderr, "\033[38;5;%dm\xE2\x96\x88", bar_color(t));  // U+2588 block
    } else {
      fprintf(stderr, "\033[38;5;240m\xE2\x94\x80");               // U+2500 dash
    }
  }
  fprintf(stderr, "\033[0m\033[1m]\033[0m %3d%%  iter %d/%d  "
                  "\033[36mvalue=%.4e\033[0m  opt=%.1e   ",
          (int)(frac * 100.0), iter, max_iter, value, inf_du);
  fflush(stderr);
}

// ---- helpers ------------------------------------------------------------

static bool same_vec(const double *a, const std::vector<double> &b, int n)
{
  for (int i = 0; i < n; i++) if (a[i] != b[i]) return false;
  return true;
}

// ---- construction ---------------------------------------------------------

ShapeTNLP::ShapeTNLP(const Parametrization &shape, const Objective &obj,
                     const double *alpha0, const double *xlo, const double *xhi,
                     const int *seeds, int nseeds,
                     const RunConfig &c, double h,
                     int max_iter, bool show_bar, double obj_scale,
                     const std::vector<const Constraint *> &constraints)
  : init_value(0.0), final_value(0.0),
    shape_(shape), obj_(obj), ndesign_(shape.ndesign()),
    c_(c), h_(h), max_iter_(max_iter), show_bar_(show_bar),
    obj_scale_(obj_scale), constraints_(constraints),
    fcache_valid_(false), fcache_f_(0.0), gcache_valid_(false),
    got_init_(false)
{
  a0_.assign(alpha0, alpha0 + ndesign_);
  xl_.assign(xlo, xlo + ndesign_);
  xu_.assign(xhi, xhi + ndesign_);
  init_alpha  = a0_;
  final_alpha = a0_;
  for (int k = 0; k < nseeds; k++) seeds_.push_back(seeds[k]);
}

// ---- objective + gradient with caching -------------------------------

double ShapeTNLP::evaluate(const double *x, bool want_grad, double *grad)
{
  if (want_grad && gcache_valid_ && same_vec(x, gcache_x_, ndesign_)) {
    for (int j = 0; j < ndesign_; j++) grad[j] = gcache_g_[j];
    if (fcache_valid_ && same_vec(x, fcache_x_, ndesign_)) return fcache_f_;
  }
  if (!want_grad && fcache_valid_ && same_vec(x, fcache_x_, ndesign_))
    return fcache_f_;

  const int ns = (int) seeds_.size();
  double v;

#ifdef SPARTA_AD
  // AD build: one SPARTA run yields value AND gradient together (the
  // Jacobian side-file / SPARTA_AD_SEED_JACFILE hook seeds all ndesign_
  // directions in that single run -- see shape_case.cpp / read_surf.cpp).
  // Gradients here are used exactly as the solver returns them: no
  // flux-measure-derivative correction, per explicit instruction -- see
  // docs/ad_phase_c_investigation/FINDINGS.md for the known ~50%
  // (specular) / ~40% (diffuse) low bias this carries.
  {
    std::vector<double> g(ndesign_);
    v = evaluate_avg(shape_, obj_, x, seeds_.data(), ns, c_,
                     want_grad ? g.data() : 0);
    fcache_x_.assign(x, x + ndesign_); fcache_f_ = v; fcache_valid_ = true;
    if (want_grad) {
      gcache_x_.assign(x, x + ndesign_);
      gcache_g_ = g; gcache_valid_ = true;
      for (int j = 0; j < ndesign_; j++) grad[j] = g[j];
    }
  }
#else
  // Stock build: gradient from the CRN central finite difference.
  if (want_grad) {
    std::vector<double> g(ndesign_);
    v = grad_fd(shape_, obj_, x, seeds_.data(), ns, c_, h_, g.data());
    fcache_x_.assign(x, x + ndesign_); fcache_f_ = v; fcache_valid_ = true;
    gcache_x_.assign(x, x + ndesign_);
    gcache_g_ = g; gcache_valid_ = true;
    for (int j = 0; j < ndesign_; j++) grad[j] = g[j];
  } else {
    v = evaluate_avg(shape_, obj_, x, seeds_.data(), ns, c_);
    fcache_x_.assign(x, x + ndesign_); fcache_f_ = v; fcache_valid_ = true;
  }
#endif

  if (!got_init_) { init_value = v; got_init_ = true; }
  return v;
}

// ---- TNLP interface ---------------------------------------------------

bool ShapeTNLP::get_nlp_info(Index &n, Index &m, Index &nnz_jac_g,
                             Index &nnz_h_lag, IndexStyleEnum &index_style)
{
  n = ndesign_;
  m = (Index) constraints_.size();
  nnz_jac_g = n * m;            // dense: trivially small (few constraints x
                                // small ndesign), no sparse structure needed
  nnz_h_lag = 0;               // Hessian via limited-memory, not supplied
  index_style = TNLP::C_STYLE;
  return true;
}

bool ShapeTNLP::get_bounds_info(Index n, Number *x_l, Number *x_u,
                                Index m, Number *g_l, Number *g_u)
{
  for (Index i = 0; i < n; i++) { x_l[i] = xl_[i]; x_u[i] = xu_[i]; }
  for (Index i = 0; i < m; i++) {
    g_l[i] = constraints_[i]->lower();
    g_u[i] = constraints_[i]->upper();
  }
  return true;
}

bool ShapeTNLP::get_scaling_parameters(Number &obj_scaling, bool &use_x_scaling,
                                       Index n, Number *x_scaling,
                                       bool &use_g_scaling, Index m,
                                       Number *g_scaling)
{
  (void) n; (void) x_scaling; (void) m; (void) g_scaling;
  obj_scaling   = obj_scale_;   // 1/|value(start)| -> scaled objective ~ O(1)
  use_x_scaling = false;        // design vars are already O(1)
  use_g_scaling = false;        // constraint values (e.g. area) are already
                                // O(1)-ish; no scaling needed for them either
  return true;
}

bool ShapeTNLP::get_starting_point(Index n, bool init_x, Number *x,
                                   bool init_z, Number *z_L, Number *z_U,
                                   Index m, bool init_lambda, Number *lambda)
{
  (void) init_z; (void) z_L; (void) z_U; (void) m; (void) init_lambda; (void) lambda;
  if (init_x) for (Index i = 0; i < n; i++) x[i] = a0_[i];
  return true;
}

bool ShapeTNLP::eval_f(Index n, const Number *x, bool new_x, Number &obj_value)
{
  (void) n; (void) new_x;
  obj_value = evaluate(x, /*want_grad=*/false, 0);
  return true;
}

bool ShapeTNLP::eval_grad_f(Index n, const Number *x, bool new_x, Number *grad_f)
{
  (void) new_x;
  std::vector<double> g(n);
  evaluate(x, /*want_grad=*/true, g.data());
  for (Index i = 0; i < n; i++) grad_f[i] = g[i];
  return true;
}

bool ShapeTNLP::eval_g(Index n, const Number *x, bool new_x, Index m, Number *g)
{
  (void) new_x;
  std::vector<double> alpha(x, x + n);
  for (Index i = 0; i < m; i++)
    g[i] = constraints_[i]->eval(shape_, alpha.data(), c_.chord, c_.nseg);
  return true;
}

bool ShapeTNLP::eval_jac_g(Index n, const Number *x, bool new_x, Index m,
                           Index nele_jac, Index *iRow, Index *jCol,
                           Number *values)
{
  (void) new_x; (void) nele_jac;
  if (values == NULL) {
    // dense triplet structure: row i (constraint i), col j (design var j)
    int k = 0;
    for (Index i = 0; i < m; i++)
      for (Index j = 0; j < n; j++) { iRow[k] = i; jCol[k] = j; k++; }
  } else {
    std::vector<double> alpha(x, x + n);
    std::vector<double> g(n);
    int k = 0;
    for (Index i = 0; i < m; i++) {
      constraints_[i]->grad(shape_, alpha.data(), c_.chord, c_.nseg, g.data());
      for (Index j = 0; j < n; j++) values[k++] = g[j];
    }
  }
  return true;
}

bool ShapeTNLP::intermediate_callback(AlgorithmMode mode, Index iter,
                                      Number obj_value, Number inf_pr,
                                      Number inf_du, Number mu, Number d_norm,
                                      Number regularization_size, Number alpha_du,
                                      Number alpha_pr, Index ls_trials,
                                      const IpoptData *ip_data,
                                      IpoptCalculatedQuantities *ip_cq)
{
  (void) mode; (void) mu; (void) d_norm; (void) regularization_size;
  (void) alpha_du; (void) alpha_pr; (void) ls_trials;

  TrajPoint p;
  p.iter   = (int) iter;
  p.value  = obj_value;
  p.inf_pr = inf_pr;
  p.inf_du = inf_du;
  p.alpha.resize(ndesign_);

  std::vector<double> x(ndesign_);
  bool ok = get_curr_iterate(ip_data, ip_cq, /*scaled=*/false, ndesign_,
                             x.data(), 0, 0, 0, 0, 0);
  if (ok) p.alpha = x;
  else if (!traj.empty()) p.alpha = traj.back().alpha;
  else p.alpha = a0_;
  traj.push_back(p);

  if (show_bar_) draw_bar((int) iter, max_iter_, obj_value, inf_du);
  return true;
}

void ShapeTNLP::finalize_solution(SolverReturn status, Index n, const Number *x,
                                  const Number *z_L, const Number *z_U, Index m,
                                  const Number *g, const Number *lambda,
                                  Number obj_value, const IpoptData *ip_data,
                                  IpoptCalculatedQuantities *ip_cq)
{
  (void) status; (void) z_L; (void) z_U; (void) m; (void) g; (void) lambda;
  (void) ip_data; (void) ip_cq;
  final_value = obj_value;
  final_alpha.resize(n);
  for (Index i = 0; i < n; i++) final_alpha[i] = x[i];
  if (show_bar_) { fprintf(stderr, "\n"); fflush(stderr); }
}
