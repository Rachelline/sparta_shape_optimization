#!/usr/bin/env python3
"""Direct test of the 'hit-count sensitivity' hypothesis: does the NUMBER
of particles that geometrically hit the front line (id=4) change with
alpha, at a magnitude comparable to the ~50% of d(press)/dalpha that AD is
missing? Uses compute_surf's raw NUM tally (a plain per-window hit count,
NOT normalized by area/fluxscale -- see compute_surf.cpp case NUM: `vec[k++]
+= 1.0`), finite-differenced across alpha on the STOCK binary (no AD
needed). Same full-ensemble deck/parameters as thinplate_derivative_match.py
so results are directly comparable to that script's FD/AD numbers.
"""
import os, re, sys, math, shutil, tempfile, subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fmf_coeff

V = 1000.0
TINF = 200.0
NRHO = 1.0e20
FNUM = 5.0e14
M = 2.325e-26
KB = 1.380649e-23
L = 0.3
W = 0.02
FRONT_LINE_ID = 4
ANGLES_DEG = [15.0, 30.0, 45.0]
FD_EPS_DEG = 0.5

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
compute         s1 surf all all press shx shy n
fix             f2 ave/surf all 1 4000 4000 c_s1[*]
timestep        1e-6
stats           4000
stats_style     step np ncoll
run             8000
dump            d1 surf all 4000 dump.surf id f_f2[1] f_f2[2] f_f2[3] f_f2[4]
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
    press, shear, num = [], [], []
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
                    shear.append(float(parts[2]))
                    num.append(float(parts[4]))
                i += 1
        else:
            i += 1
    return press, num


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
    press, num = parse_dump(os.path.join(workdir, "dump.surf"))
    return sum(press)/len(press), sum(num)/len(num)


def main():
    spa = os.environ.get("SPARTA_STOCK")
    if not spa:
        print("ERROR: set SPARTA_STOCK")
        return 1

    vmp = math.sqrt(2*KB*TINF/M)
    s = V/vmp
    q = 0.5*NRHO*M*V*V
    workdir = tempfile.mkdtemp(prefix="hitcount_")
    shutil.copy(os.path.join(HERE, "N.species"), workdir)

    print(f"{'alpha':>6} {'press(a)':>10} {'hits(a)':>10} {'dP/da(FD)':>11} "
          f"{'dN/da(FD)':>11} {'mean_p/hit':>11} {'dN/da*mp/hit':>13} {'implied_%_of_gap':>16}")
    for deg in ANGLES_DEG:
        alpha = math.radians(deg)
        eps = math.radians(FD_EPS_DEG)

        p0, n0 = run_value(spa, alpha, workdir)
        p_plus, n_plus = run_value(spa, alpha+eps, workdir)
        p_minus, n_minus = run_value(spa, alpha-eps, workdir)

        dpress_fd = (p_plus - p_minus)/(2*eps)
        dnum_fd = (n_plus - n_minus)/(2*eps)
        mean_p_per_hit = p0/n0 if n0 else float('nan')
        # naive product-rule decomposition: dP/da = dN/da*(P/N) + N*d(P/N)/da
        # the first term is the "count/capture" contribution the AD single-
        # realization derivative structurally cannot see; report it as a
        # fraction of the total FD derivative to see if it's ~half.
        count_term = dnum_fd * mean_p_per_hit
        frac_of_fd = count_term/dpress_fd if dpress_fd else float('nan')

        print(f"{deg:6.1f} {p0:10.5f} {n0:10.3f} {dpress_fd:11.5f} "
              f"{dnum_fd:11.5f} {mean_p_per_hit:11.6f} {count_term:13.5f} {frac_of_fd:16.2%}")

    shutil.rmtree(workdir, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
