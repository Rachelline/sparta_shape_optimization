#!/usr/bin/env python3
"""
fmf_coeff.py -- closed-form free-molecular-flow (FMF) aerodynamic coefficients
for a flat surface panel, used as a ZERO-NOISE analytic ground truth for
validating SPARTA's AD (Sacado) drag derivatives.

WHY THIS EXISTS
---------------
In the free-molecular limit (no gas-gas collisions) the force on a flat panel
is an exact integral over the incoming Maxwellian, with a closed form. That
gives an analytic reference for both the drag VALUE and, by differentiation,
the drag DERIVATIVE d(drag)/d(design) -- something finite-differencing a noisy
DSMC run cannot provide. See docs/PLAN.md ("AD-vs-analytic-FMF" test ladder).
This module is Phase A: the analytic reference itself, self-verified, with no
dependency on SPARTA. Phases B-D compare it against stock/AD SPARTA runs.

WORKFLOW INTENT
---------------
Standard library only (math) so it runs in CI with no install step. Running
this file executes its self-test and exits 0 (pass) / 1 (fail); it is intended
to be registered as a CI step so the ground truth can never silently drift.

FORMULAS (verbatim structure from ADBSat, an open-source, peer-reviewed FMF
panel-method tool from the University of Manchester, GPL-3.0):
  - Schaaf & Chambre model: coeff_schaaf.m
  - permalink (commit e40b10d):
    https://github.com/nhcrisp/ADBSat/blob/e40b10def4705a9d36e43c316cd9b4cd2030558f/toolbox/fmf_eq/coeff_schaaf.m
Original derivations (cited in that file):
  - Schaaf, S.A. & Chambre, P.L., "Flow of Rarefied Gases", Princeton Univ.
    Press, 1961 (general sigmaN/sigmaT form).
  - Sentman, L.H., "Free Molecule Flow Theory and Its Application to the
    Determination of Aerodynamic Forces", Lockheed report LMSC-448514, 1961.
Independent V&V of the ADBSat implementation:
  - Crisp et al., "ADBSat: Verification and validation of a novel panel method
    for quick aerodynamic analysis of satellites", Comput. Phys. Commun. 275
    (2022) 108309. https://doi.org/10.1016/j.cpc.2022.108309

CONVENTIONS
-----------
  delta : angle between the freestream velocity and the OUTWARD surface normal
          [rad]. delta = 0 is normal incidence (flow hits the panel head-on).
  s     : molecular speed ratio, s = V_inf / v_mp, with the most-probable speed
          v_mp = sqrt(2 k T_inf / m). Computable from a RunConfig's vstream,
          species mass, and freestream temperature.
  sigmaN, sigmaT : normal / tangential momentum accommodation coefficients.
          Specular reflection = (0, 0); fully diffuse = (1, 1).
  TwTi  : wall-to-freestream temperature ratio Tw / T_inf (diffuse only; the
          specular coefficients do not depend on it).
Coefficients cp (pressure, along -normal), ctau (shear, along the surface),
cd (drag, along the flow) are per the ADBSat sign conventions.
"""

import math

SQRT_PI = math.sqrt(math.pi)


def cp_schaaf(delta, s, sigmaN=0.0, TwTi=1.0):
    """Pressure coefficient (Schaaf & Chambre). Specular: sigmaN=0."""
    c = math.cos(delta)
    e = math.exp(-s * s * c * c)
    ep = 1.0 + math.erf(s * c)
    return (1.0 / s**2) * (
        ((2 - sigmaN) * s / SQRT_PI * c + sigmaN / 2 * math.sqrt(TwTi)) * e
        + ((2 - sigmaN) * (0.5 + s**2 * c**2)
           + sigmaN / 2 * math.sqrt(TwTi) * SQRT_PI * s * c) * ep)


def ctau_schaaf(delta, s, sigmaT=0.0):
    """Shear (tangential) coefficient. Specular: sigmaT=0 -> ctau == 0."""
    c = math.cos(delta)
    e = math.exp(-s * s * c * c)
    ep = 1.0 + math.erf(s * c)
    return sigmaT * math.sin(delta) / (s * SQRT_PI) * (e + s * SQRT_PI * c * ep)


def cd_coeff(delta, s, sigmaN=0.0, sigmaT=0.0, TwTi=1.0):
    """Drag coefficient (component along the freestream)."""
    return (cp_schaaf(delta, s, sigmaN, TwTi) * math.cos(delta)
            + ctau_schaaf(delta, s, sigmaT) * math.sin(delta))


def cl_coeff(delta, s, sigmaN=0.0, sigmaT=0.0, TwTi=1.0):
    """Lift coefficient (component perpendicular to the freestream)."""
    return (cp_schaaf(delta, s, sigmaN, TwTi) * math.sin(delta)
            - ctau_schaaf(delta, s, sigmaT) * math.cos(delta))


def dcd_ddelta_specular(delta, s):
    """
    Analytic d(cd)/d(delta) for the SPECULAR case (sigmaN=sigmaT=0).
    Derived symbolically (sympy) and transcribed; the self-test below verifies
    it against a central finite difference of cd_coeff(). This is the term the
    Phase C/D chain rule multiplies by d(delta)/d(design).
    """
    c = math.cos(delta)
    sd = math.sin(delta)
    e = math.exp(-s * s * c * c)
    ep = 1.0 + math.erf(s * c)
    term = 2 * s * e * c / SQRT_PI + (2 * s**2 * c**2 + 1) * ep
    dterm = (4 * s**3 * e * sd * c**2 / SQRT_PI
             - 4 * s**2 * ep * sd * c
             - 2 * s * (2 * s**2 * c**2 + 1) * e * sd / SQRT_PI
             - 2 * s * e * sd / SQRT_PI)
    return -term * sd / s**2 + dterm * c / s**2


# ------------------------------------------------------------------ self-test
def _self_test():
    ok = True

    def approx(a, b, tol, what):
        nonlocal ok
        good = abs(a - b) <= tol
        print(f"  [{'PASS' if good else 'FAIL'}] {what}: {a:.10f} vs {b:.10f} "
              f"(|d|={abs(a-b):.2e}, tol={tol:g})")
        ok = ok and good

    # 1. specular shear must be identically zero
    approx(ctau_schaaf(0.9, 2.0, sigmaT=0.0), 0.0, 0.0, "specular ctau == 0")

    # 2. analytic derivative vs central FD of the closed form (the permanent
    #    guard on the transcription); noise-free because cd is closed form
    h = 1e-6
    maxerr = 0.0
    for s in (0.5, 2.0, 5.0):
        for i in range(1, 6):          # 15 deg .. 75 deg
            d = i * math.pi / 12
            fd = (cd_coeff(d + h, s) - cd_coeff(d - h, s)) / (2 * h)
            an = dcd_ddelta_specular(d, s)
            maxerr = max(maxerr, abs(fd - an))
    approx(maxerr, 0.0, 1e-6, "analytic dcd/ddelta vs FD-of-formula (max err)")

    # 3. hardcoded regression anchors (independently computed via sympy)
    approx(cd_coeff(0.0, 2.0),                 4.499808588970, 1e-9, "cd(0, s=2)")
    approx(cd_coeff(math.pi/4, 2.0),           1.765727400077, 1e-9, "cd(pi/4, s=2)")
    approx(dcd_ddelta_specular(0.0, 2.0),      0.0,            1e-9, "dcd(0, s=2) == 0")
    approx(dcd_ddelta_specular(math.pi/4, 2.0), -4.606162191618, 1e-8, "dcd(pi/4, s=2)")
    approx(dcd_ddelta_specular(math.pi/6, 5.0), -4.540000000034, 1e-8, "dcd(pi/6, s=5)")
    # diffuse anchor (sigmaN=sigmaT=1, Tw==Tinf) -- future diffuse-wall tests
    approx(cd_coeff(math.pi/4, 2.0, 1.0, 1.0, 1.0), 2.037967032969, 1e-9,
           "cd diffuse(pi/4, s=2)")

    return ok


if __name__ == "__main__":
    import sys
    print("FMF closed-form coefficient self-test:")
    sys.exit(0 if _self_test() else 1)
