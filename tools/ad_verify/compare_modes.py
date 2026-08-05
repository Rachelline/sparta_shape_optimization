#!/usr/bin/env python3
"""
compare_modes.py -- run one shape/objective/alpha config three ways
(finite-difference, AD-uncorrected, AD-corrected) through ad_src's
shape_main / shape_main_ad binaries, and tabulate value, gradient, and
the AD/FD ratio per component.

This is the direct answer to "run with and without the fix": one
command, one table, showing exactly what --score-correction changes (or
doesn't -- e.g. it's a documented no-op for --objective drag, since drag
tallies fx via nfactor_inverse, not fluxscale; see docs/AD_GRADIENTS.md).

Usage:
  cd ad_src && cmake --build build   # need both shape_main and shape_main_ad
  python3 tools/ad_verify/compare_modes.py --alpha 1.3,1.0,2.7,0.8

  # heatflux needs a thermally-accommodating wall to be non-degenerate:
  python3 tools/ad_verify/compare_modes.py --objective heatflux \
      --wall-accom 1.0 --alpha 1.3,1.0,2.7,0.8

Binary paths default to ad_src/build/{shape_main,shape_main_ad} relative
to this script; override with --stock-bin/--ad-bin if built elsewhere.
Extra flags after "--" are forwarded verbatim to both binaries (e.g.
--seeds, --nsettle, --navg, --h).
"""

import argparse
import re
import subprocess
import sys
from pathlib import Path

GRAD_RE = re.compile(
    r"grad_(?:fd|ad) \([^)]*\):\s*\[(?P<grad>.*?)\]\s*"
    r"\(base value used:\s*(?P<value>[^)]+)\)"
)


def run(binary, args, cwd):
    cmd = [str(binary)] + args
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=600, cwd=cwd)
    except FileNotFoundError:
        sys.exit(f"ERROR: binary not found: {binary}")
    if proc.returncode != 0:
        sys.exit(f"ERROR: {' '.join(cmd)} failed (exit {proc.returncode}):\n"
                 f"{proc.stdout}\n{proc.stderr}")
    m = GRAD_RE.search(proc.stdout)
    if not m:
        sys.exit(f"ERROR: could not parse gradient line from:\n{proc.stdout}")
    value = float(m.group("value"))
    grad = [float(x) for x in m.group("grad").split(",")]
    return value, grad


def fmt_vec(v):
    return "[" + ", ".join(f"{x: .6e}" for x in v) + "]"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--alpha", required=True)
    ap.add_argument("--objective", default="drag")
    ap.add_argument("--shape", default="bezier")
    ap.add_argument("--stock-bin", default=None,
                    help="path to shape_main (default: ad_src/build/shape_main)")
    ap.add_argument("--ad-bin", default=None,
                    help="path to shape_main_ad (default: ad_src/build/shape_main_ad)")
    ap.add_argument("extra", nargs="*",
                    help="extra flags forwarded to both binaries, e.g. "
                         "-- --seeds 12345,777 --wall-accom 1.0")
    args = ap.parse_args()

    ad_src = Path(__file__).resolve().parents[2] / "ad_src"
    stock_bin = Path(args.stock_bin) if args.stock_bin else ad_src / "build" / "shape_main"
    ad_bin = Path(args.ad_bin) if args.ad_bin else ad_src / "build" / "shape_main_ad"

    common = ["--alpha", args.alpha, "--objective", args.objective,
             "--shape", args.shape] + args.extra

    print(f"config: objective={args.objective} shape={args.shape} alpha={args.alpha}"
         + (f" extra={' '.join(args.extra)}" if args.extra else ""))
    print(f"  stock binary: {stock_bin}")
    print(f"  AD binary:    {ad_bin}\n")

    fd_value, fd_grad = run(stock_bin, common, ad_src)
    ad_unc_value, ad_unc_grad = run(ad_bin, common, ad_src)
    ad_cor_value, ad_cor_grad = run(ad_bin, common + ["--score-correction"], ad_src)

    rows = [
        ("FD (grad_fd)",         fd_value,     fd_grad),
        ("AD (uncorrected)",     ad_unc_value, ad_unc_grad),
        ("AD (score-corrected)", ad_cor_value, ad_cor_grad),
    ]

    print(f"{'mode':<22} {'value':>16}   gradient")
    for name, value, grad in rows:
        print(f"{name:<22} {value: .6e}   {fmt_vec(grad)}")

    print("\nAD/FD ratio per component (uncorrected, corrected):")
    for j in range(len(fd_grad)):
        r_unc = ad_unc_grad[j] / fd_grad[j] if fd_grad[j] != 0.0 else float("nan")
        r_cor = ad_cor_grad[j] / fd_grad[j] if fd_grad[j] != 0.0 else float("nan")
        print(f"  component {j}: uncorrected={r_unc: .4f}   corrected={r_cor: .4f}")

    if ad_unc_grad == ad_cor_grad:
        print("\nNOTE: uncorrected and corrected AD gradients are identical -- "
             f"--score-correction is a no-op for --objective {args.objective} "
             "(expected for tallies routed through nfactor_inverse, e.g. "
             "drag's fx; see docs/AD_GRADIENTS.md's coverage table).")


if __name__ == "__main__":
    main()
