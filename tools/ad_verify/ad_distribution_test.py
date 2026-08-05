#!/usr/bin/env python3
"""
ad_distribution_test.py -- two-sample distributional equivalence check
(permutation Kolmogorov-Smirnov) for the AD build's stochastic tests.

WHY THIS EXISTS
  ad_stochastic_equivalence.py (TOST-style mean equivalence) and
  ad_convergence_sweep.py (1/sqrt(N) trend) both test the MEAN of the AD
  ensemble against stock's. Neither would necessarily catch a bug that
  leaves the mean roughly unchanged but distorts the shape of the
  distribution -- e.g. a bug that only sometimes triggers a worse cascade
  after an RNG desync, producing a heavier tail or a bimodal split, without
  moving the average much. This script complements the mean-based tests by
  comparing the full empirical distributions.

  This is NOT a replacement for the mean-equivalence tests, and it is not,
  by itself, proof of equivalence -- failing to reject "same distribution"
  is exactly as unable to prove equivalence as failing to reject "same
  mean" was before TOST was adopted for that check (see
  ad_stochastic_equivalence.py). There is no widely-used equivalence
  formulation for a full-distribution test the way TOST exists for a mean.
  Its role here is a REGRESSION SCREEN: it is more sensitive than a mean
  check to shape/variance/tail differences, so it is the more likely of the
  three tests to catch a *future* AD-path change that breaks something in a
  way that does not show up as a simple mean shift.

METHOD
  Collect N independent stock samples and N independent AD samples per
  tracked quantity (reusing ad_convergence_sweep.run_pool). The seed VALUES
  are shared between the stock and AD pools purely for run-to-run
  reproducibility/bookkeeping -- they must NOT be treated as paired
  observations. Live-code tracing (see docs/PLAN.md) already established
  that AD and stock trajectories desync at the first decision-boundary
  flip and evolve on independent RNG streams after that, so stock[seed=X]
  and AD[seed=X] are no more correlated than stock[seed=X] and
  AD[seed=Y]. This script always treats the two samples as independent
  (unpaired) sets, which is why it draws on two DISTINCT seed blocks
  (stock uses seed_base+i, AD uses seed_base+AD_OFFSET+i) rather than
  literally reusing identical seed values -- removing any temptation to
  pair by seed index downstream.

  For each tracked quantity, compute the two-sample Kolmogorov-Smirnov
  statistic D (the max gap between the two empirical CDFs), then estimate
  its p-value by permutation: pool both samples, reshuffle into two groups
  of the original sizes many times, recompute D each time, and take
  p = fraction of shuffles with D >= the observed D. A permutation test
  (rather than the classical asymptotic KS p-value) is used because Ncoll
  is integer-valued with ties, which the asymptotic KS formula assumes
  away; the permutation approach makes no distributional assumptions and
  handles ties correctly.

  D itself is reported alongside the p-value as an effect size: D is the
  largest fraction of probability mass, at any point, that sits on the
  "wrong side" between the two empirical CDFs -- e.g. D=0.05 means the two
  distributions never disagree by more than 5 percentage points anywhere.

  Which columns get tested is auto-detected from the stock pool (every
  stats column except Step/CPU that isn't constant across that pool -- see
  select_tracked_columns() in ad_stochastic_equivalence.py), or set
  explicitly with --columns. Not hardcoded, since different decks report
  different columns (e.g. in.shocktube has no c_temp; in.surf_react_heatflux
  has no Natt).

USAGE
  SPARTA_STOCK=/path/to/stock/spa_<machine> \\
  SPARTA_AD=/path/to/ad/spa_<machine> \\
  python3 ad_distribution_test.py [--n 64] [--case examples/bfield/in.bfield] \\
      [--jobs 6] [--n-perm 10000] [--alpha 0.05] [--out-json dist.json]

WORKFLOW INTENT (not yet wired into CI)
  Same two-binary requirement and SPARTA_STOCK/SPARTA_AD interface as
  ad_stochastic_equivalence.py and ad_convergence_sweep.py. Standard
  library only (random.shuffle for the permutation loop, no third-party
  dependency), so it drops into the same future two-binary CI job as the
  other two scripts without extra installs.
"""
import os
import sys
import json
import random
import argparse

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from ad_convergence_sweep import run_pool  # noqa: E402
from ad_stochastic_equivalence import (  # noqa: E402
    transpose, select_tracked_columns, REPO_ROOT, DEFAULT_TIMEOUT)

AD_OFFSET = 5_000_000  # keeps the AD pool's seeds distinct from stock's --
                       # see METHOD above: samples are always unpaired


def ks_statistic(a, b):
    """Two-sample KS D: max gap between empirical CDFs of a and b."""
    vals = sorted(set(a) | set(b))
    na, nb = len(a), len(b)
    sa, sb = sorted(a), sorted(b)

    def cdf_at(sorted_vals, x):
        # fraction of sorted_vals <= x, via binary search
        import bisect
        return bisect.bisect_right(sorted_vals, x) / len(sorted_vals)

    d = 0.0
    for x in vals:
        d = max(d, abs(cdf_at(sa, x) - cdf_at(sb, x)))
    return d


def permutation_test(a, b, n_perm, rng):
    observed = ks_statistic(a, b)
    pooled = list(a) + list(b)
    na = len(a)
    count_ge = 0
    for _ in range(n_perm):
        rng.shuffle(pooled)
        pa, pb = pooled[:na], pooled[na:]
        if ks_statistic(pa, pb) >= observed:
            count_ge += 1
    p_value = count_ge / n_perm
    return observed, p_value


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=64, help="samples per side")
    ap.add_argument("--case", default="examples/bfield/in.bfield")
    ap.add_argument("--seed-base", type=int, default=8880001)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--n-perm", type=int, default=10000)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--rng-seed", type=int, default=12345,
                     help="seed for the permutation-test shuffler (reproducibility)")
    ap.add_argument("--out-json", default=None)
    ap.add_argument("--columns", default=None,
                     help="comma-separated column names to test (default: "
                          "auto-detect every non-constant column in the stats output)")
    ap.add_argument("--timeout", type=int, default=DEFAULT_TIMEOUT,
                     help=f"per-run subprocess timeout in seconds (default {DEFAULT_TIMEOUT})")
    args = ap.parse_args()

    stock = os.environ.get("SPARTA_STOCK")
    ad = os.environ.get("SPARTA_AD")
    if not stock or not ad:
        print("ERROR: set SPARTA_STOCK and SPARTA_AD to the two binaries to compare")
        return 1

    case_path = os.path.join(REPO_ROOT, args.case)
    case_dir = os.path.dirname(case_path)
    with open(case_path) as f:
        base_deck = f.read()

    stock_seeds = [args.seed_base + i for i in range(args.n)]
    ad_seeds = [args.seed_base + AD_OFFSET + i for i in range(args.n)]

    print(f"case={args.case}  n={args.n} per side  n_perm={args.n_perm}  alpha={args.alpha}  "
          f"timeout={args.timeout}s")
    print(f"running {args.n} stock seeds...")
    stock_r = run_pool(stock, case_dir, base_deck, stock_seeds, args.jobs, args.timeout)
    print(f"running {args.n} AD seeds...")
    ad_r = run_pool(ad, case_dir, base_deck, ad_seeds, args.jobs, args.timeout)

    requested = args.columns.split(",") if args.columns else None
    tracked = select_tracked_columns(transpose(stock_r), requested)
    print(f"tracked columns: {tracked}")

    rng = random.Random(args.rng_seed)
    ok = True
    results = {}
    print(f"\n{'quantity':<10} {'KS D':>8} {'p-value':>10}  result")
    for q in tracked:
        a = [r[q] for r in stock_r]
        b = [r[q] for r in ad_r]
        d, p = permutation_test(a, b, args.n_perm, rng)
        good = p >= args.alpha
        ok = ok and good
        results[q] = {"D": d, "p_value": p, "pass": good}
        print(f"{q:<10} {d:8.4f} {p:10.4f}  {'PASS' if good else 'FAIL'}")

    print(f"\ndistribution screen ({args.n} samples/side, alpha={args.alpha}): "
          f"{'ALL PASS' if ok else 'FAIL'}")
    print("Note: PASS is a regression screen (no distributional difference detected "
          "at this sample size), not a proof of equivalence. See "
          "ad_stochastic_equivalence.py / ad_convergence_sweep.py for the mean-"
          "equivalence tests this complements.")

    if args.out_json:
        out = {
            "case": args.case,
            "n": args.n,
            "n_perm": args.n_perm,
            "alpha": args.alpha,
            "tracked": tracked,
            "raw": {"stock": stock_r, "AD": ad_r},
            "results": results,
        }
        with open(args.out_json, "w") as f:
            json.dump(out, f, indent=2)
        print(f"wrote {args.out_json}")

    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
