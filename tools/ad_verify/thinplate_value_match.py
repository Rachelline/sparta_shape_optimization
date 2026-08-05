#!/usr/bin/env python3
"""Value-match check for the small CLOSED thin-rectangle plate geometry
(clear of all boundaries). Freestream stays FIXED along +x; alpha tilts
the plate itself. Extracts ONLY the front (windward) surf element's
press via a dump, NOT a reduce-sum across all 4 sides -- summing raw
(intensive, per-area) pressure across front/back/edge-caps is physically
meaningless, since they see completely different flow conditions.
Confirmed via a direct per-element dump: line id 4 is consistently the
front face (matches analytic almost exactly at alpha=0); id 2 is the
shadowed back (reads exactly 0); ids 1,3 are the thin edge caps (small
nonzero -- a real artifact of a finite-thickness rectangle, not part of
the idealized infinite-plate physics being tested)."""
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
ANGLES_DEG = [0.0, 30.0, 45.0, 60.0]
PRESS_RTOL = 0.05
SHEAR_ATOL_FRAC = 0.05
L = 0.3
W = 0.02
FRONT_LINE_ID = 4

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
read_surf       data.thinplate_{deg}
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
    press, shear = [], []
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
                i += 1
        else:
            i += 1
    return press, shear


def main():
    spa = os.environ.get("SPARTA_STOCK")
    if not spa:
        print("ERROR: set SPARTA_STOCK")
        return 1

    vmp = math.sqrt(2*KB*TINF/M)
    s = V/vmp
    q = 0.5*NRHO*M*V*V
    print(f"speed ratio s = {s:.4f}  q_inf = {q:.6g} Pa  plate L={L} w={W}  "
          f"front line id={FRONT_LINE_ID}")
    print(f"{'alpha':>6} {'press_DSMC':>11} {'press_ana':>11} {'rel_err':>9} "
          f"{'shear/q':>9} {'ncoll':>6}  result")

    workdir = tempfile.mkdtemp(prefix="thinplate_vm2_")
    shutil.copy(os.path.join(HERE, "N.species"), workdir)
    ok = True
    for deg in ANGLES_DEG:
        data_path = os.path.join(workdir, f"data.thinplate_{int(deg)}")
        with open(data_path, "w") as f:
            f.write(thin_plate(math.radians(deg)))
        deck = DECK.format(nrho=NRHO, fnum=FNUM, v=V, tinf=TINF, deg=int(deg))
        infile = os.path.join(workdir, "in.case")
        with open(infile, "w") as f:
            f.write(deck)
        out = subprocess.run([spa, "-in", infile], cwd=workdir, capture_output=True, text=True)
        if out.returncode != 0:
            raise RuntimeError("SPARTA failed:\n" + out.stdout[-2000:] + out.stderr[-2000:])

        ncoll_max = 0
        for line in out.stdout.splitlines():
            c = line.split()
            if c and re.fullmatch(r"\d+", c[0]) and len(c) == 3:
                ncoll_max = max(ncoll_max, int(c[2]))

        press, shear = parse_dump(os.path.join(workdir, "dump.surf"))
        pd = sum(press)/len(press)
        sh = sum(shear)/len(shear) if shear else 0.0
        d = math.radians(deg)
        pa = fmf_coeff.cp_schaaf(d, s) * q
        rel = abs(pd-pa)/pa
        shq = abs(sh)/q
        good = (rel < PRESS_RTOL and shq < SHEAR_ATOL_FRAC and ncoll_max == 0)
        ok = ok and good
        print(f"{deg:6.1f} {pd:11.5f} {pa:11.5f} {rel:9.2%} {shq:9.2%} "
              f"{ncoll_max:6d}  {'PASS' if good else 'FAIL'}")
    shutil.rmtree(workdir, ignore_errors=True)
    print("\nthin-plate value match (front face only):", "ALL PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
