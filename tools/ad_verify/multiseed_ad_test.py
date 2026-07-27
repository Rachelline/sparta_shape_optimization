#!/usr/bin/env python3
"""Multi-SEED averaging of the AD derivative (as opposed to the TIME
averaging within a single run that thinplate_derivative_match.py already
does). Runs the full AD ensemble test at ONE fixed angle across several
independent RNG seeds and averages dP/dalpha across seeds.

Purpose: distinguish NOISE from BIAS in the ~50% AD-vs-FD gap. If the gap is
statistical noise, multi-seed averaging should pull the AD mean toward the
FD/analytic value as more seeds are added. If the gap is a structural bias
(e.g. AD missing the particle-hit/miss "capture-count" derivative term --
see conversation), multi-seed averaging will NOT close the gap: the AD mean
will converge tightly, just to the same ~50%-of-truth value, with shrinking
seed-to-seed spread but an unchanged offset from FD.
"""
import os, re, sys, math, shutil, tempfile, subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fmf_coeff

KB = 1.380649e-23
V = 1000.0
TINF = 200.0
NRHO = 1.0e20
FNUM = 5.0e14
M = 2.325e-26
ALPHA_DEG = 30.0
L = 0.3
W = 0.02
FRONT_LINE_ID = 4
SEEDS = [12345, 54321, 98765, 24681, 13579]

DECK = """seed            {seed}
dimension       2
global          gridcut 0.0 comm/sort yes
boundary        oo pp pp
create_box      -1 1 -1 1 -0.5 0.5
create_grid     21 23 1
balance_grid    rcb cell
global          nrho {nrho} fnum {fnum}
species         N.species N
mixture         gas N vstream {v} 0 0 temp {tinf}
read_surf       data.thinplate
surf_collide    1 specular
surf_modify     all collide 1
fix             in emit/face gas xlo twopass
compute         s1 surf all all press shx shy
fix             f2 ave/surf all 1 4000 4000 c_s1[*]
timestep        1e-6
stats           4000
stats_style     step np ncoll
run             8000
dump            d1 surf all 4000 dump.surf id f_f2[1] f_f2[2]
dump_modify     d1 buffer no
run             24000
"""


def thin_plate(alpha, L=L, w=W):
    corners = list(reversed([(-w/2, -L/2), (w/2, -L/2), (w/2, L/2), (-w/2, L/2)]))
    rot = []
    for x, y in corners:
        xr = x*math.cos(alpha) - y*math.sin(alpha)
        yr = x*math.sin(alpha) + y*math.cos(alpha)
        rot.append((xr, yr))
    pts = "\n".join(f"{i+1} {x!r} {y!r}" for i, (x, y) in enumerate(rot))
    lines = "\n".join(f"{i+1} {i+1} {(i+1)%4+1}" for i in range(4))
    return f"""small closed thin-rectangle plate, tilt alpha

4 points
4 lines

Points

{pts}

Lines

{lines}
"""


def run_ad_derivative(spa_ad, alpha, seed, workdir):
    data_path = os.path.join(workdir, "data.thinplate")
    with open(data_path, "w") as f:
        f.write(thin_plate(alpha))
    deck = DECK.format(seed=seed, nrho=NRHO, fnum=FNUM, v=V, tinf=TINF)
    infile = os.path.join(workdir, "in.case")
    with open(infile, "w") as f:
        f.write(deck)
    env = dict(os.environ)
    env["SPARTA_AD_SEED_ALPHA"] = "1"
    env["SPARTA_AD_SEED_DIR"] = "0"
    env["SPARTA_AD_DUMP_DX"] = "0"
    out = subprocess.run([spa_ad, "-in", infile], cwd=workdir, capture_output=True,
                          text=True, env=env)
    if out.returncode != 0:
        raise RuntimeError("SPARTA(AD) failed:\n" + out.stdout[-2000:] + out.stderr[-2000:])

    blocks = []
    cur = []
    for line in out.stderr.splitlines():
        m = re.match(r"AD_DX row=(\d+) col=(\d+) val=([\-0-9.eE+]+) dx=([\-0-9.eE+]+)", line)
        if not m:
            continue
        row, col, val, dx = int(m.group(1)), int(m.group(2)), float(m.group(3)), float(m.group(4))
        if row == 0 and col == 0 and cur:
            blocks.append(cur)
            cur = []
        cur.append((row, col, val, dx))
    if cur:
        blocks.append(cur)
    if not blocks:
        raise RuntimeError("No AD_DX lines found in stderr:\n" + out.stderr[-3000:])

    vals, dxs = [], []
    for block in blocks:
        id_val = None
        for row, col, val, dx in block:
            if col == 0 and int(val) == FRONT_LINE_ID:
                id_val = row
        for row, col, val, dx in block:
            if row == id_val and col == 1:
                vals.append(val)
                dxs.append(dx)
    if not vals:
        raise RuntimeError(f"front line id {FRONT_LINE_ID} not found in any AD_DX block")
    return sum(vals)/len(vals), sum(dxs)/len(dxs)


def main():
    spa_ad = os.environ.get("SPARTA_AD")
    if not spa_ad:
        print("ERROR: set SPARTA_AD")
        return 1

    vmp = math.sqrt(2*KB*TINF/M)
    s = V/vmp
    q = 0.5*NRHO*M*V*V
    alpha = math.radians(ALPHA_DEG)
    h = 1e-6
    dcp_fd = (fmf_coeff.cp_schaaf(alpha+h, s) - fmf_coeff.cp_schaaf(alpha-h, s)) / (2*h)
    ana_ref = dcp_fd * q
    # reference AD_vs_FD numbers already measured (single seed 12345, from
    # thinplate_derivative_match.py): dP/da(FD) = -3.95496 at alpha=30deg
    fd_ref = -3.95496

    workdir = tempfile.mkdtemp(prefix="multiseed_ad_")
    shutil.copy(os.path.join(HERE, "N.species"), workdir)

    print(f"alpha={ALPHA_DEG} deg   analytic ref dP/da={ana_ref:.5f}   "
          f"single-seed FD ref dP/da={fd_ref:.5f}")
    print(f"{'seed':>8} {'press':>10} {'dP/da (AD)':>12} {'running_mean':>13} {'running_std':>12}")

    dxs = []
    for seed in SEEDS:
        press, dx = run_ad_derivative(spa_ad, alpha, seed, workdir)
        dxs.append(dx)
        mean = sum(dxs)/len(dxs)
        std = (sum((x-mean)**2 for x in dxs)/len(dxs))**0.5 if len(dxs) > 1 else 0.0
        print(f"{seed:8d} {press:10.5f} {dx:12.5f} {mean:13.5f} {std:12.5f}")

    mean_ad = sum(dxs)/len(dxs)
    print(f"\nfinal {len(dxs)}-seed AD mean dP/da = {mean_ad:.5f}")
    print(f"  vs single-seed FD  = {fd_ref:.5f}  (ratio AD/FD = {mean_ad/fd_ref:.4f})")
    print(f"  vs analytic        = {ana_ref:.5f}  (ratio AD/ana = {mean_ad/ana_ref:.4f})")

    shutil.rmtree(workdir, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
