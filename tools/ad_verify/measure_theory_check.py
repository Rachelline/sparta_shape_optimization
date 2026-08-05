#!/usr/bin/env python3
"""
measure_theory_check.py -- self-test for the FLUX-MEASURE DERIVATIVE that
pathwise (forward-mode) AD structurally cannot see.

WHY THIS EXISTS
---------------
Phase C found SPARTA's Sacado derivative d(press)/d(alpha) to be exactly
HALF the true value: ratios 0.4979 / 0.5012 / 0.5025 / 0.5120 across every
angle, both plate lengths, and a 5-seed ensemble (sigma ~ 0.0016 -- a bias,
not noise). This file explains and pins that factor, with NO SPARTA
involved: just sampling from a drifting Maxwellian.

THE ARGUMENT
------------
One specular hit deposits normal momentum 2m|u|, u = v.n (see compute_surf.cpp
case PRESS, and MathExtra::reflect3: v -= 2(v.n)n). But hits are NOT drawn
uniformly from the freestream -- a surface is struck at a rate proportional to
the incoming flux |u|, so the hit population carries density ~ |u| f(v):

    E[press] ~ INT_{u<0} |u| f(v) * 2m|u| dv  =  2m INT u^2 f dv
                         ^^^^^^^^ sampling weight
                                   ^^^^^^^ momentum per hit

BOTH factors depend on alpha, since both contain n(alpha). Forward AD
differentiates only the second: in a single realization the hit set is already
drawn and frozen, so the sampling weight's alpha-dependence is invisible to it.
The integrand is QUADRATIC in |u| with each factor LINEAR, so the two
contributions are exactly equal:

    pathwise / true  ==  1/2      (exactly; any speed ratio, any angle)

This is the standard score-function / likelihood-ratio term:

    d/da E_p(a)[g(a)]  =  E[dg/da]   +   E[g * dlog p/da]
                          ^^^^^^^^^ AD    ^^^^^^^^^^^^^^^ AD is blind to this

CONSEQUENCE (see docs/ad_phase_c_investigation/FINDINGS.md)
    Doubling AD's answer recovers the truth to <2.4% on every case measured --
    but the clean factor of 2 is SPECIFIC to specular + normal pressure, where
    momentum-per-hit is exactly linear in |u|. For diffuse walls (thermal
    re-emission adds an |u|-independent piece) or for shear, the two terms
    differ and the general score term is required. Do NOT hardcode a 2 in the
    solver; this test exists partly to catch anyone who tries.

Standard library only, so it runs in CI with no install step. Running this
file executes the self-test and exits 0 (pass) / 1 (fail).
"""

import argparse
import math
import random
import sys

KB = 1.380649e-23
M = 2.325e-26
TINF = 200.0
V = 1000.0

# Predicted ratio, and how far it may drift before we call it a failure.
# Tolerance is dominated by Monte-Carlo sampling error in the central
# difference, not by anything physical.
EXPECTED_RATIO = 0.5
RATIO_TOL = 5e-3
ANGLES_DEG = (15.0, 30.0, 45.0, 60.0)


def moments(pop, alpha):
    """Given a FIXED particle population (one 'realization'), return
    (press_integrand, pathwise_derivative_integrand) at tilt `alpha`.

    The population is deliberately reused across alpha values -- that is
    exactly what a single AD trace does, and it is what makes the sampling
    weight's alpha-dependence invisible to it."""
    ca, sa = math.cos(alpha), math.sin(alpha)
    nx, ny = -ca, -sa            # outward normal, facing the +x freestream
    tx, ty = sa, -ca             # dn/dalpha (unit tangent)
    press = 0.0
    ad = 0.0
    for vx, vy in pop:
        u = vx * nx + vy * ny
        if u >= 0.0:
            continue             # receding: never strikes this face
        zeta = vx * tx + vy * ty         # = du/dalpha
        w = -u                           # |u| : the flux sampling weight
        press += w * (2.0 * w)           # weight * momentum-per-hit
        ad += w * (-2.0 * zeta)          # weight * d(momentum-per-hit)/dalpha
    n = len(pop)
    return press / n, ad / n


def _self_test(nsamples, seed, verbose=True):
    vmp = math.sqrt(2 * KB * TINF / M)
    sigma = vmp / math.sqrt(2.0)         # per-component thermal std dev
    s = V / vmp
    rng = random.Random(seed)
    pop = [(V + rng.gauss(0, sigma), rng.gauss(0, sigma))
           for _ in range(nsamples)]

    if verbose:
        print(f"speed ratio s = {s:.4f}   N = {nsamples:,} particles   "
              f"seed = {seed}")
        print(f"{'alpha':>7} {'press':>14} {'pathwise(AD)':>14} "
              f"{'true(FD)':>14} {'AD/true':>9}  result")

    ok = True
    h = 1e-4
    for deg in ANGLES_DEG:
        a = math.radians(deg)
        p0, ad = moments(pop, a)
        pp, _ = moments(pop, a + h)
        pm, _ = moments(pop, a - h)
        # true derivative: the WHOLE integral moves, sampling weight included
        true = (pp - pm) / (2 * h)
        ratio = ad / true
        good = abs(ratio - EXPECTED_RATIO) <= RATIO_TOL
        ok = ok and good
        if verbose:
            print(f"{deg:7.1f} {p0:14.5g} {ad:14.5g} {true:14.5g} "
                  f"{ratio:9.4f}  {'PASS' if good else 'FAIL'}")
    return ok


def main():
    ap = argparse.ArgumentParser(
        description="Self-test: pathwise AD recovers exactly half of a "
                    "flux-weighted surface derivative.")
    ap.add_argument("-n", "--nsamples", type=int, default=2_000_000,
                    help="particles sampled (default 2000000)")
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()

    print("flux-measure derivative self-test "
          f"(expect AD/true = {EXPECTED_RATIO} +/- {RATIO_TOL}):")
    ok = _self_test(args.nsamples, args.seed)
    print("\nresult:", "ALL PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
