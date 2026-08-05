#!/usr/bin/env python3
"""
flatplate_value_match.py -- Phase B of the AD-vs-analytic FMF ladder.

Validates that STOCK SPARTA (no AD needed) reproduces the closed-form
free-molecular-flow (FMF) panel pressure. A uniform single-species freestream
strikes an axis-aligned SPECULAR wall (one-sided, so it matches the one-sided
FMF formula); the incidence angle delta is set by tilting the freestream, and
the wall's normal pressure is compared to cp(delta, s) * q_inf from
fmf_coeff.py. Specular => shear must be ~0 at every angle.

WHY THIS MATTERS
  This is the VALUE-match gate: it pins down units, the q_inf normalization,
  the speed-ratio convention, and the angle convention against a zero-noise
  analytic reference -- WITHOUT any AD. Only once this passes is the Phase C/D
  derivative comparison (Sacado d(drag)/d(alpha) vs analytic) meaningful.
  See docs/PLAN.md, "AD-vs-analytic-FMF verification ladder".

USAGE
  SPARTA=/path/to/spa_<machine>  python3 flatplate_value_match.py
  (exit 0 = all angles within tolerance; 1 = mismatch or setup error)

WORKFLOW INTENT
  Standard library only. Registerable as a CI job on the stock build:
  build spa, then run this. Fast (a few short single-species runs).
"""
import os
import re
import sys
import math
import shutil
import tempfile
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fmf_coeff  # noqa: E402

# ---- physical constants / flow conditions (single source of truth) ----
KB = 1.380649e-23           # Boltzmann [J/K]
V = 1000.0                  # freestream speed [m/s]
TINF = 200.0                # freestream temperature [K]
NRHO = 1.0e20               # number density [1/m^3]
FNUM = 5.0e14               # real/sim particle ratio
ANGLES_DEG = [0.0, 30.0, 45.0, 60.0]   # incidence angles to test
PRESS_RTOL = 0.02           # 2% (Monte-Carlo sampling floor over the window)
SHEAR_ATOL_FRAC = 0.02      # |shear| must be < 2% of q_inf (specular => ~0)


def species_mass(path):
    """Parse the single-species Molmass (kg) from an N.species file."""
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                return float(line.split()[2])
    raise RuntimeError("no species line found in " + path)


DECK = """seed            12345
dimension       2
global          gridcut 0.0 comm/sort yes
boundary        os pp pp
create_box      0 1 0 1 -0.5 0.5
create_grid     20 20 1
balance_grid    rcb cell
global          nrho {nrho} fnum {fnum}
species         {species} N
mixture         gas N vstream {vx} {vy} 0 temp {tinf}
surf_collide    1 specular
bound_modify    xhi collide 1
fix             in emit/face gas xlo twopass
compute         cb boundary all press shx shy
variable        pxhi equal c_cb[2][1]
variable        sxhi equal c_cb[2][2]
timestep        1e-6
stats           2000
stats_style     step np ncoll v_pxhi
run             6000
fix             fp ave/time 1 2000 2000 v_pxhi v_sxhi
stats_style     step np ncoll v_pxhi f_fp[1] f_fp[2]
run             8000
"""


def run_case(spa, workdir, species_path, delta):
    vx = V * math.cos(delta)
    vy = V * math.sin(delta)
    deck = DECK.format(nrho=NRHO, fnum=FNUM, species=os.path.basename(species_path),
                       vx=repr(vx), vy=repr(vy), tinf=TINF)
    infile = os.path.join(workdir, "in.case")
    with open(infile, "w") as f:
        f.write(deck)
    out = subprocess.run([spa, "-in", infile], cwd=workdir,
                         capture_output=True, text=True)
    if out.returncode != 0:
        raise RuntimeError("SPARTA failed:\n" + out.stdout[-2000:] + out.stderr[-2000:])
    # parse rows of the averaging run: cols step np ncoll v_pxhi f_fp[1] f_fp[2]
    press, shear, ncoll_max = [], [], 0
    for line in out.stdout.splitlines():
        c = line.split()
        if len(c) == 6 and re.fullmatch(r"\d+", c[0]):
            step = int(c[0])
            ncoll_max = max(ncoll_max, int(c[2]))
            if step >= 8000:                 # steady, averaging window
                press.append(float(c[4]))
                shear.append(float(c[5]))
    if not press:
        raise RuntimeError("no averaging-window rows parsed")
    return sum(press) / len(press), sum(shear) / len(shear), ncoll_max


def main():
    spa = os.environ.get("SPARTA")
    if not spa or not shutil.which(spa) and not os.path.exists(spa):
        print("ERROR: set SPARTA=/path/to/spa_<machine> (stock build)")
        return 1
    species_path = os.path.join(HERE, "N.species")
    m = species_mass(species_path)
    vmp = math.sqrt(2 * KB * TINF / m)
    s = V / vmp
    q = 0.5 * NRHO * m * V * V
    print(f"speed ratio s = {s:.4f},  q_inf = {q:.6g} Pa,  m = {m:g} kg")
    print(f"{'delta':>6} {'press_DSMC':>11} {'press_ana':>11} {'rel_err':>9} "
          f"{'shear/q':>9} {'ncoll':>6}  result")

    ok = True
    workdir = tempfile.mkdtemp(prefix="fmf_plate_")
    shutil.copy(species_path, workdir)
    try:
        for deg in ANGLES_DEG:
            d = math.radians(deg)
            pd, sh, ncoll = run_case(spa, workdir, species_path, d)
            pa = fmf_coeff.cp_schaaf(d, s) * q
            rel = abs(pd - pa) / pa
            shq = abs(sh) / q
            good = (rel < PRESS_RTOL and shq < SHEAR_ATOL_FRAC and ncoll == 0)
            ok = ok and good
            print(f"{deg:6.1f} {pd:11.5f} {pa:11.5f} {rel:9.2%} {shq:9.2%} "
                  f"{ncoll:6d}  {'PASS' if good else 'FAIL'}")
    finally:
        shutil.rmtree(workdir, ignore_errors=True)

    print("\nPhase B flat-plate value match:", "ALL PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
