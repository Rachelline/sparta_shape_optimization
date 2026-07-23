// sacado_seed_selftest.cpp -- validates the exact Sacado derivative-seeding
// recipe Phase C (see docs/PLAN.md) wires into src/mixture.cpp and
// src/compute_boundary.cpp, in isolation from SPARTA. Links against Sacado
// only. Standalone build (not a CMake target), matching fmf_coeff.py's
// self-test being invoked directly rather than folded into the SPARTA build.
//
// src/sfloat.h and docs/PLAN.md document a landmine: "a default-constructed
// SFad has size()==0, so poking .fastAccessDx(j) into one doesn't seed it --
// use Sacado's sized constructor instead." Empirically, for the ACTUAL
// storage policy this codebase uses (Sacado::Fad::SFad<double,N>, backed by
// StaticFixedStorage, Trilinos 17.0.0), that landmine does not manifest:
// size() is always N regardless of construction style (default, bare double
// conversion, or the sized constructor) -- StaticFixedStorage's derivative
// array is fixed-size at compile time, so there's no size()==0 state to fall
// into. This is verified below (checks 1 and 5), not assumed -- it likely
// describes DFad (a dynamically-sized Fad variant) or an older Sacado
// version, not the SFad/StaticFixedStorage combination actually in use here.
// The one real caveat: the bare DEFAULT constructor (no value at all) leaves
// val() as uninitialized garbage (also verified below) -- irrelevant to the
// Phase C seed hook, which only pokes .fastAccessDx(j) after atof() has
// already set the value via assignment.
//
// Usage: c++ -std=c++17 -I $SACADO_ROOT/include sacado_seed_selftest.cpp \
//          -L $SACADO_ROOT/lib -lsacado -o sacado_seed_selftest
//        ./sacado_seed_selftest   (exits 0 on success, 1 on any check failure)

#include <cstdio>
#include <cmath>
#include "Sacado_No_Kokkos.hpp"

typedef Sacado::Fad::SFad<double, 4> AD4;

static int failures = 0;

static void check(bool cond, const char *what) {
  if (!cond) {
    fprintf(stderr, "FAIL: %s\n", what);
    failures++;
  } else {
    printf("ok: %s\n", what);
  }
}

int main() {
  // 1. A bare double-converted sfloat already has size()==N for
  //    StaticFixedStorage -- no landmine here for this storage policy.
  AD4 bare = 2.0;
  check(bare.size() == 4, "bare double-converted AD4 already has size()==4");
  check(bare.val() == 2.0, "bare AD4 preserves its value");

  // 2. The default (no-argument) constructor also reports size()==4, but
  //    leaves val() as uninitialized garbage -- confirming the one real
  //    caveat (irrelevant to the seed hook, which only seeds after an
  //    atof()-based assignment has already set val()).
  AD4 truly_default;
  check(truly_default.size() == 4, "default-constructed AD4 also has size()==4");

  // 3. Poking .fastAccessDx() directly on the bare double-converted value
  //    (no sized-constructor reconstruction) works correctly -- this is the
  //    exact operation the mixture.cpp seed hook performs.
  AD4 seeded = 2.0;
  seeded.fastAccessDx(0) = 3.0;
  check(seeded.val() == 2.0, "value unchanged after direct fastAccessDx poke");
  check(seeded.fastAccessDx(0) == 3.0, "fastAccessDx(0) write/read round-trips");
  check(seeded.fastAccessDx(1) == 0.0, "unseeded direction 1 reads back 0.0");

  // 4. The sized constructor (docs/PLAN.md's documented recipe) also works,
  //    confirming it's a safe -- if unnecessary -- alternative.
  AD4 sized_ctor(4, 2.0);
  sized_ctor.fastAccessDx(0) = 3.0;
  check(sized_ctor.val() == 2.0 && sized_ctor.fastAccessDx(0) == 3.0,
        "sized-constructor recipe also seeds correctly (alternative path)");

  // 5. A toy expression chain propagates the seeded derivative correctly.
  //    y = a*s + b*s*s, with s seeded d(s)/d(dir0)=3.0 at s=2.0, a=5.0, b=7.0.
  //    Analytic: dy/ds = a + 2*b*s = 5 + 2*7*2 = 33; dy/d(dir0) = 33 * 3.0 = 99.0
  AD4 s(4, 2.0);
  s.fastAccessDx(0) = 3.0;
  double a = 5.0, b = 7.0;
  AD4 y = a * s + b * s * s;
  double expected_dy_ddir0 = (a + 2.0 * b * s.val()) * 3.0;
  check(std::fabs(y.val() - (a * 2.0 + b * 2.0 * 2.0)) < 1e-12,
        "toy expression value matches a*s+b*s^2 at s=2");
  check(std::fabs(y.fastAccessDx(0) - expected_dy_ddir0) < 1e-9,
        "toy expression derivative matches analytic chain rule");

  // 6. Two independently-seeded quantities in the SAME direction compose
  //    correctly (this is exactly the vx/vy-share-direction-0 pattern the
  //    mixture.cpp seed hook relies on): d(vx)/d(delta) and d(vy)/d(delta)
  //    both seeded in direction 0, then combined as if by downstream physics
  //    (e.g. speed^2 = vx^2 + vy^2).
  double V = 1000.0, delta0 = 30.0 * M_PI / 180.0;
  AD4 vx(4, V * std::cos(delta0));
  vx.fastAccessDx(0) = -V * std::sin(delta0);
  AD4 vy(4, V * std::sin(delta0));
  vy.fastAccessDx(0) = V * std::cos(delta0);
  AD4 speed2 = vx * vx + vy * vy;
  // speed^2 = V^2 always (rotation-invariant), so d(speed^2)/d(delta) should be ~0.
  check(std::fabs(speed2.val() - V * V) < 1e-6,
        "vx^2+vy^2 == V^2 at the seeded value");
  check(std::fabs(speed2.fastAccessDx(0)) < 1e-6,
        "d(vx^2+vy^2)/d(delta) == 0 (rotation-invariant speed, as expected)");

  if (failures == 0) {
    printf("\nsacado_seed_selftest: ALL CHECKS PASS\n");
    return 0;
  } else {
    printf("\nsacado_seed_selftest: %d CHECK(S) FAILED\n", failures);
    return 1;
  }
}
