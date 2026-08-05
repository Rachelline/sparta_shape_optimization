/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer

   sfloat: project-wide scalar type for the physics core (AD seam).

   Stock build (default)   : sfloat == double, zero overhead, behavior
                             bit-identical to upstream SPARTA.
   AD build (-DSPARTA_AD)   : sfloat == Sacado forward-mode dual number
                             (Sacado::Fad::SFad<double, SPARTA_AD_NDIR>)
                             with SPARTA_AD_NDIR derivative directions
                             (default 4 = the Bezier design parameters
                             alpha = x1,y1,x2,y2). SPARTA_AD_NDIR is set
                             by the build (see cmake/common/process/
                             sparta_build_options.cmake).

   Why Sacado (vs a hand-rolled dual): Sacado's SFad is a mature,
   Kokkos-aware forward-AD type. When SPARTA_KOKKOS is defined we include
   the Kokkos-aware umbrella (Sacado.hpp) so SFad's methods are device-
   callable; otherwise the lighter Sacado_No_Kokkos.hpp. Sacado and
   SPARTA MUST share ONE Kokkos of the same version -- see the one-Kokkos
   guard in sparta_build_options.cmake.

   Public contract relied on by the (to-be-converted) SPARTA core:
   - implicit construction FROM double (Sacado: derivative = 0)
   - value extraction via spval(x) / x.val(); derivatives END at spval
   - comparisons compare VALUES only (branches follow the primal flow;
     Sacado's Fad comparison operators already do this)
   - spval(x): identity for passive types, x.val() for sfloat
   - MPI_SFLOAT selects the datatype per build (see STUBS/mpi.h). Under
     Sacado the layout is NOT "5 doubles"; the MPI stub datatype must be
     sized as sizeof(sfloat) -- see STUBS/mpi.* (follow-up).

   Passive by design (never sfloat): RNG internals/outputs, timers,
   binary file formats, ubuf bit-punning unions.

   NOTE (follow-ups before the full double->sfloat conversion links):
   - Seeding a design-variable derivative must go through Sacado's sized
     constructor (a default-constructed SFad has size()==0); see the
     ad_seed() helper (ad_src) rather than poking .dx(j) into a bare SFad.
   - Math-function coverage below is the delta over what Sacado already
     provides (erf/erfc/tgamma/lgamma etc.); verified by compiling this
     header against the Sacado headers.
------------------------------------------------------------------------- */

#ifndef SPARTA_SFLOAT_H
#define SPARTA_SFLOAT_H

#include <cmath>

#ifndef SPARTA_AD

/* ---------------- stock build: plain double ---------------- */

typedef double sfloat;
template <typename T> inline T spval(T x) { return x; }   // identity
#define MPI_SFLOAT MPI_DOUBLE

#else

/* ---------------- AD build: Sacado forward-mode SFad ---------------- */

#ifndef SPARTA_AD_NDIR
#define SPARTA_AD_NDIR 4
#endif

#include <type_traits>
#ifdef SPARTA_KOKKOS
#include "Sacado.hpp"            // Kokkos-aware: SFad usable in device kernels
#else
#include "Sacado_No_Kokkos.hpp"  // no Kokkos dependency
#endif

typedef Sacado::Fad::SFad<double, SPARTA_AD_NDIR> sfloat;

// value extraction: identity for passive (arithmetic) types, .val() for
// sfloat AND any Sacado Fad EXPRESSION (a-b, a*c, ... are expression-template
// types, not concrete sfloat; they convert to sfloat here so .val() runs).
// The passive template is constrained to arithmetic types so it does NOT
// greedily capture Sacado expressions (which would return the expression
// instead of its double value -- the reference's concrete sfloat never had
// expression templates, so this constraint is Sacado-specific).
// derivatives deliberately end here, mirroring the stock spval contract.
template <typename T>
inline typename std::enable_if<std::is_arithmetic<T>::value, T>::type
spval(T x) { return x; }
inline double spval(const sfloat &x) { return x.val(); }

// ---- math functions Sacado does NOT overload for Fad types ----
// Sacado already provides sqrt, fabs, abs, exp, log, log10, sin, cos,
// tan, asin, acos, atan, atan2, pow, ceil, and the comparison ops; do
// NOT redefine those. Only the chain-rule cases below are added, written
// against Sacado's public API (.val()/.size()/.fastAccessDx()).

inline sfloat erf(const sfloat &a) {
  sfloat r(a);
  const double v = a.val();
  const double g = 2.0/std::sqrt(M_PI)*std::exp(-v*v);
  r.val() = std::erf(v);
  for (int i = 0; i < a.size(); i++) r.fastAccessDx(i) = g*a.fastAccessDx(i);
  return r;
}
inline sfloat erfc(const sfloat &a) {
  sfloat r(a);
  const double v = a.val();
  const double g = -2.0/std::sqrt(M_PI)*std::exp(-v*v);
  r.val() = std::erfc(v);
  for (int i = 0; i < a.size(); i++) r.fastAccessDx(i) = g*a.fastAccessDx(i);
  return r;
}

// digamma via recurrence + asymptotic series (for tgamma/lgamma derivs)
inline double sfloat_digamma(double x) {
  double r = 0.0;
  while (x < 6.0) { r -= 1.0/x; x += 1.0; }
  double f = 1.0/(x*x);
  return r + std::log(x) - 0.5/x
         - f*(1.0/12.0 - f*(1.0/120.0 - f/252.0));
}
inline sfloat tgamma(const sfloat &a) {
  sfloat r(a);
  r.val() = std::tgamma(a.val());
  const double g = r.val()*sfloat_digamma(a.val());
  for (int i = 0; i < a.size(); i++) r.fastAccessDx(i) = g*a.fastAccessDx(i);
  return r;
}
inline sfloat lgamma(const sfloat &a) {
  sfloat r(a);
  r.val() = std::lgamma(a.val());
  const double g = sfloat_digamma(a.val());
  for (int i = 0; i < a.size(); i++) r.fastAccessDx(i) = g*a.fastAccessDx(i);
  return r;
}

// derivative-zero-a.e. functions: derivative deliberately dropped
inline sfloat floor(const sfloat &a) { return sfloat(std::floor(a.val())); }
inline sfloat ceil(const sfloat &a)  { return sfloat(std::ceil(a.val())); }
inline sfloat round(const sfloat &a) { return sfloat(std::round(a.val())); }
inline sfloat trunc(const sfloat &a) { return sfloat(std::trunc(a.val())); }
inline sfloat fmod(const sfloat &a, const sfloat &b) {
  sfloat r(a);                       // d/da fmod = 1 a.e.; keep a's derivs
  r.val() = std::fmod(a.val(), b.val());
  return r;
}

// value-based selects / predicates (Sacado does not overload these for Fad)
inline sfloat fmax(const sfloat &a, const sfloat &b) { return a.val() >= b.val() ? a : b; }
inline sfloat fmin(const sfloat &a, const sfloat &b) { return a.val() <= b.val() ? a : b; }
inline sfloat abs(const sfloat &a)  { return a.val() >= 0.0 ? a : -a; }
inline bool isnan(const sfloat &a)  { return std::isnan(a.val()); }
inline bool isinf(const sfloat &a)  { return std::isinf(a.val()); }

#define MPI_SFLOAT MPI_5DOUBLE

#endif  /* SPARTA_AD */

#endif
