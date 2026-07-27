#!/usr/bin/env python3
"""Deterministic, single-particle synthetic test isolating the AD chain
through surf geometry (read_surf point seeding) -> specular reflection ->
compute_surf/fix_ave_surf tally -> dump extraction, with ZERO randomness
(no emission, no particle-particle collisions, exactly one particle placed
by hand via `create_particles ... single`). This directly tests whether the
AD *machinery* is correct, independent of any grid-topology / statistical
noise questions raised for the full DSMC ensemble test.

Since there is no RNG anywhere in this deck, a plain central finite
difference on the STOCK binary with a TINY eps (no statistical floor) is an
essentially exact ground truth -- if AD disagrees with FD here, the bug is
in the AD chain itself (reflection formula, normal computation, or tally),
not in grid/emission-rate discreteness.
"""
import os, sys, math, shutil, tempfile, subprocess

HERE = os.path.dirname(os.path.abspath(__file__))

L = 0.3
W = 0.02
FRONT_LINE_ID = 4
V = 1000.0
UPSTREAM_DIST = 0.5
NSTEPS = 12
DT = 1e-4  # particle travels V*DT = 0.1 per step -> covers UPSTREAM_DIST*~2 over NSTEPS

DECK = """seed            12345
dimension       2
global          gridcut 0.0 comm/sort yes
boundary        oo pp pp
create_box      -1 1 -1 1 -0.5 0.5
create_grid     21 23 1
balance_grid    rcb cell
global          nrho 1.0e20 fnum 1.0
species         N.species N
mixture         gas N
read_surf       data.thinplate
surf_collide    1 specular
surf_modify     all collide 1
create_particles gas single N {x0!r} {y0!r} 0.0 {vx!r} 0.0 0.0
compute         s1 surf all all press shx shy
fix             f2 ave/surf all 1 {nsteps} {nsteps} c_s1[*]
timestep        {dt!r}
stats           1
stats_style     step np ncoll
run             {nsteps}
dump            d1 surf all {nsteps} dump.surf id f_f2[1] f_f2[2]
dump_modify     d1 buffer no
run             0
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
    return f"""single-particle synthetic test plate, tilt alpha

4 points
4 lines

Points

{pts}

Lines

{lines}
""", rot


def front_line_midpoint(rot):
    # line4 connects point4 -> point1 in the 1-indexed point list (see
    # thin_plate()/lines construction): rot[3] -> rot[0]
    x1, y1 = rot[3]
    x2, y2 = rot[0]
    return (x1+x2)/2.0, (y1+y2)/2.0


def build_deck(alpha, workdir):
    plate_txt, rot = thin_plate(alpha)
    with open(os.path.join(workdir, "data.thinplate"), "w") as f:
        f.write(plate_txt)
    mx, my = front_line_midpoint(rot)
    x0 = mx - UPSTREAM_DIST
    y0 = my
    deck = DECK.format(x0=x0, y0=y0, vx=V, nsteps=NSTEPS, dt=DT)
    with open(os.path.join(workdir, "in.case"), "w") as f:
        f.write(deck)
    return x0, y0


def parse_dump(path, front_id=FRONT_LINE_ID):
    press = None
    with open(path) as f:
        lines = f.read().splitlines()
    i = 0
    while i < len(lines):
        if lines[i].startswith("ITEM: SURFS"):
            i += 1
            for _ in range(4):
                parts = lines[i].split()
                if int(parts[0]) == front_id:
                    press = float(parts[1])
                i += 1
        else:
            i += 1
    return press


def run_value(spa, alpha, workdir):
    build_deck(alpha, workdir)
    out = subprocess.run([spa, "-in", "in.case"], cwd=workdir, capture_output=True, text=True)
    if out.returncode != 0:
        raise RuntimeError("SPARTA failed:\n" + out.stdout[-3000:] + out.stderr[-3000:])
    ncoll = None
    for line in out.stdout.splitlines():
        c = line.split()
        if len(c) == 3 and c[0].lstrip('-').isdigit():
            ncoll = int(c[2])
    press = parse_dump(os.path.join(workdir, "dump.surf"))
    return press, ncoll


def run_ad(spa_ad, alpha, workdir):
    build_deck(alpha, workdir)
    env = dict(os.environ)
    env["SPARTA_AD_SEED_ALPHA"] = "1"
    env["SPARTA_AD_SEED_DIR"] = "0"
    env["SPARTA_AD_DUMP_DX"] = "0"
    out = subprocess.run([spa_ad, "-in", "in.case"], cwd=workdir, capture_output=True,
                          text=True, env=env)
    if out.returncode != 0:
        raise RuntimeError("SPARTA(AD) failed:\n" + out.stdout[-3000:] + out.stderr[-3000:])
    import re
    press_val = press_dx = None
    for line in out.stderr.splitlines():
        m = re.match(r"AD_DX row=(\d+) col=(\d+) val=([\-0-9.eE+]+) dx=([\-0-9.eE+]+)", line)
        if not m:
            continue
        row, col, val, dx = int(m.group(1)), int(m.group(2)), float(m.group(3)), float(m.group(4))
        if col == 0 and int(val) == FRONT_LINE_ID:
            front_row = row
    for line in out.stderr.splitlines():
        m = re.match(r"AD_DX row=(\d+) col=(\d+) val=([\-0-9.eE+]+) dx=([\-0-9.eE+]+)", line)
        if not m:
            continue
        row, col, val, dx = int(m.group(1)), int(m.group(2)), float(m.group(3)), float(m.group(4))
        if row == front_row and col == 1:
            press_val, press_dx = val, dx
    ncoll = None
    for line in out.stdout.splitlines():
        c = line.split()
        if len(c) == 3 and c[0].lstrip('-').isdigit():
            ncoll = int(c[2])
    return press_val, press_dx, ncoll


def main():
    spa = os.environ.get("SPARTA_STOCK")
    spa_ad = os.environ.get("SPARTA_AD")
    if not spa or not spa_ad:
        print("ERROR: set SPARTA_STOCK and SPARTA_AD")
        return 1

    alpha0_deg = 30.0
    alpha0 = math.radians(alpha0_deg)
    workdir = tempfile.mkdtemp(prefix="synthetic1p_")
    shutil.copy(os.path.join(HERE, "N.species"), workdir)

    press_ad, dpress_ad, ncoll_ad = run_ad(spa_ad, alpha0, workdir)
    print(f"AD:  alpha={alpha0_deg} press={press_ad:.10g} dP/dalpha={dpress_ad:.10g} ncoll={ncoll_ad}")

    print(f"{'eps(rad)':>12} {'dP/dalpha (FD)':>16} {'rel_err_vs_AD':>14} {'ncoll+':>7} {'ncoll-':>7}")
    for eps in [1e-2, 1e-3, 1e-4, 1e-5, 1e-6]:
        p_plus, nc_plus = run_value(spa, alpha0+eps, workdir)
        p_minus, nc_minus = run_value(spa, alpha0-eps, workdir)
        if p_plus is None or p_minus is None:
            print(f"{eps:12.1e}  MISS (no collision recorded on front line at this eps)")
            continue
        dfd = (p_plus - p_minus) / (2*eps)
        rel = abs(dfd - dpress_ad) / (abs(dpress_ad) + 1e-30)
        print(f"{eps:12.1e} {dfd:16.10g} {rel:14.2%} {nc_plus:7} {nc_minus:7}")

    shutil.rmtree(workdir, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
