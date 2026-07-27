#!/usr/bin/env python3
"""Derivative-match check (Phase C, plate-rotation variant) for the small
closed thin-rectangle plate geometry. Freestream stays FIXED along +x;
alpha tilts the plate itself via a rigid-body rotation of its 4 corner
points. This sidesteps the emission-rate-truncation bug found in the
earlier vstream-seeded (delta) version entirely: ntarget in fix emit/face
never depends on alpha, since the freestream is never seeded -- only the
plate's own point coordinates carry a derivative.

AD side: SPARTA_AD_SEED_ALPHA (any nonempty value) + SPARTA_AD_SEED_DIR
trigger read_surf.cpp's point-coordinate seed hook (see src/read_surf.cpp,
read_points()): every point (x,y) read from the surf file gets
d(x,y)/dalpha = (-y,x), i.e. the point rotated a further 90 degrees. This
derivative survives untouched through fix ave/surf's pure-sfloat
accumulation and is extracted at the last possible point -- inside
dump_surf.cpp's write_text(), gated by SPARTA_AD_DUMP_DX -- printed to
stderr as AD_DX lines, bypassing the normal spval()-only file output.

Reference side: central finite difference of the STOCK (non-AD) binary
across alpha, and the closed-form analytic derivative
d(cp_schaaf)/ddelta (with delta == alpha, confirmed by the Phase-value
match script using the same convention) from fmf_coeff.dcd_ddelta_specular.
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
ANGLES_DEG = [15.0, 30.0, 45.0]
L = 0.3
W = 0.02
FRONT_LINE_ID = 4
FD_EPS_DEG = 0.5
REL_TOL = 0.15

DECK = """seed            12345
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


def parse_dump(path, front_id=FRONT_LINE_ID):
    press = []
    with open(path) as f:
        lines = f.read().splitlines()
    i = 0
    while i < len(lines):
        if lines[i].startswith("ITEM: SURFS"):
            i += 1
            for _ in range(4):
                parts = lines[i].split()
                if int(parts[0]) == front_id:
                    press.append(float(parts[1]))
                i += 1
        else:
            i += 1
    return press


def run_value(spa, alpha, workdir):
    data_path = os.path.join(workdir, "data.thinplate")
    with open(data_path, "w") as f:
        f.write(thin_plate(alpha))
    deck = DECK.format(nrho=NRHO, fnum=FNUM, v=V, tinf=TINF)
    infile = os.path.join(workdir, "in.case")
    with open(infile, "w") as f:
        f.write(deck)
    out = subprocess.run([spa, "-in", infile], cwd=workdir, capture_output=True, text=True)
    if out.returncode != 0:
        raise RuntimeError("SPARTA failed:\n" + out.stdout[-2000:] + out.stderr[-2000:])
    press = parse_dump(os.path.join(workdir, "dump.surf"))
    return sum(press)/len(press)


def run_ad_derivative(spa_ad, alpha, workdir):
    data_path = os.path.join(workdir, "data.thinplate")
    with open(data_path, "w") as f:
        f.write(thin_plate(alpha))
    deck = DECK.format(nrho=NRHO, fnum=FNUM, v=V, tinf=TINF)
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

    # AD_DX row=<i> col=<j> val=<v> dx=<d> ; col 1 = press (id,press,shear)
    # last block of AD_DX lines corresponds to the FINAL dump (after run 24000)
    press_val = press_dx = None
    id_val = None
    blocks = []
    cur = []
    prev_row = -1
    for line in out.stderr.splitlines():
        m = re.match(r"AD_DX row=(\d+) col=(\d+) val=([\-0-9.eE+]+) dx=([\-0-9.eE+]+)", line)
        if not m:
            continue
        row, col, val, dx = int(m.group(1)), int(m.group(2)), float(m.group(3)), float(m.group(4))
        if row == 0 and col == 0 and cur:
            blocks.append(cur)
            cur = []
        cur.append((row, col, val, dx))
        prev_row = row
    if cur:
        blocks.append(cur)
    if not blocks:
        raise RuntimeError("No AD_DX lines found in stderr:\n" + out.stderr[-3000:])

    # average over ALL dump snapshots (each is itself a fix ave/surf window),
    # matching run_value()/parse_dump()'s averaging over every ITEM: SURFS
    # block in dump.surf, to keep noise levels comparable between AD and FD.
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
    press_val = sum(vals)/len(vals)
    press_dx = sum(dxs)/len(dxs)
    return press_val, press_dx


def main():
    spa = os.environ.get("SPARTA_STOCK")
    spa_ad = os.environ.get("SPARTA_AD")
    if not spa or not spa_ad:
        print("ERROR: set SPARTA_STOCK and SPARTA_AD")
        return 1

    vmp = math.sqrt(2*KB*TINF/M)
    s = V/vmp
    q = 0.5*NRHO*M*V*V
    print(f"speed ratio s = {s:.4f}  q_inf = {q:.6g} Pa")
    print(f"{'alpha':>6} {'dP/da (AD)':>12} {'dP/da (FD)':>12} {'dP/da (ana)':>12} "
          f"{'AD_vs_FD':>9} {'AD_vs_ana':>10}  result")

    workdir = tempfile.mkdtemp(prefix="thinplate_deriv_")
    shutil.copy(os.path.join(HERE, "N.species"), workdir)
    ok = True
    for deg in ANGLES_DEG:
        alpha = math.radians(deg)
        eps = math.radians(FD_EPS_DEG)

        press_ad, dpress_ad = run_ad_derivative(spa_ad, alpha, workdir)

        p_plus = run_value(spa, alpha + eps, workdir)
        p_minus = run_value(spa, alpha - eps, workdir)
        dpress_fd = (p_plus - p_minus) / (2*eps)

        dpress_ana = fmf_coeff.dcd_ddelta_specular(alpha, s) * q

        rel_fd = abs(dpress_ad - dpress_fd) / abs(dpress_fd)
        rel_ana = abs(dpress_ad - dpress_ana) / abs(dpress_ana)
        good = rel_fd < REL_TOL and rel_ana < REL_TOL
        ok = ok and good
        print(f"{deg:6.1f} {dpress_ad:12.5f} {dpress_fd:12.5f} {dpress_ana:12.5f} "
              f"{rel_fd:9.2%} {rel_ana:10.2%}  {'PASS' if good else 'FAIL'}")
        print(f"       (AD value check: press={press_ad:.5f} vs analytic "
              f"{fmf_coeff.cp_schaaf(alpha, s)*q:.5f})")

    shutil.rmtree(workdir, ignore_errors=True)
    print("\nthin-plate derivative match (alpha-seeded, plate rotation):",
          "ALL PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
