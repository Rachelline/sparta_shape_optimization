#!/usr/bin/env python3
"""
ad_convergence_sweep.py -- 1/sqrt(N) convergence check for the AD build's
stochastic (RNG-branch-sensitive) tests.

WHY THIS EXISTS
  ad_stochastic_equivalence.py checks equivalence at a single, fixed N. That
  answers "is the AD-vs-stock mean difference within noise at this sample
  size" but not "is this noise, or a small systematic bias that a bigger N
  would have caught". By the Central Limit Theorem, the standard error of a
  sampling-noise mean difference shrinks as 1/sqrt(N) and the mean difference
  itself converges to 0. A systematic bias (wrong probability, dropped energy
  term, sign error) does NOT shrink -- it converges to a nonzero constant.
  Sweeping N and plotting the trend distinguishes the two cases directly,
  without needing to pick an equivalence margin up front.

METHOD
  Build three seed-indexed pools once, at the largest N in the sweep:
    stock_A -- stock binary,  seeds = seed_base + i
    stock_B -- stock binary,  seeds = seed_base + STOCK_B_OFFSET + i
    ad_A    -- AD binary,     seeds = seed_base + i          (same seeds as stock_A)
  For each sweep checkpoint N (nested/cumulative -- the first N seeds of each
  pool, not a fresh independent draw), compute two quantities per tracked
  field:
    noise_floor(N) = mean(stock_B[:N]) - mean(stock_A[:N])   -- stock vs itself,
                      the "meaningless difference" yardstick
    ad_bias(N)      = mean(ad_A[:N])    - mean(stock_A[:N])   -- the candidate
                      AD-vs-stock difference under test
  Using nested subsets (not fresh samples per checkpoint) keeps the sequence
  smooth as N grows instead of introducing independent sampling noise between
  checkpoints.

  |ad_bias(N)| trending toward 0 at roughly the same rate as |noise_floor(N)|
  (both ~ 1/sqrt(N)) is the signature of pure sampling noise. |ad_bias(N)|
  flattening at a level above |noise_floor(N)| as N grows is the signature of
  a systematic bias.

  PASS/FAIL: a z-test (|ad_bias(N_max)/sem_bias(N_max)| < Z_THRESHOLD) at
  the LARGEST N in the sweep -- the most statistically powerful point,
  since sem shrinks as 1/sqrt(N). A bias that still fails this after having
  had the most data to shrink into noise is the flattening-not-shrinking
  signature of a real bug. This is a genuine gate (exit 1 on failure), not
  only a printed diagnostic -- used as one of two required checks (with
  ad_stochastic_equivalence.py) by ctest_fallback_gate.py.

  Which columns get swept is auto-detected from stock_A (every stats column
  except Step/CPU that isn't constant across that pool -- see
  select_tracked_columns() in ad_stochastic_equivalence.py), or set
  explicitly with --columns. Not hardcoded, since different decks report
  different columns.

USAGE
  SPARTA_STOCK=/path/to/stock/spa_<machine> \\
  SPARTA_AD=/path/to/ad/spa_<machine> \\
  python3 ad_convergence_sweep.py [--sweep-ns 4,8,16,32,64] \\
      [--case examples/bfield/in.bfield] [--jobs 6] [--out-json sweep.json]

WORKFLOW INTENT
  Wired into CI via tools/ad_verify/ctest_fallback_gate.py, run from
  .github/workflows/ad.yml's mpi-stubs-ad-sacado job when ctest reports a
  failure: needs a stock build and an -DSPARTA_ENABLE_AD build in the same
  job (that job builds both). Standard library only (uses a thread pool
  for parallel runs, no third-party dependency). Emits JSON so a separate
  step can render it (e.g. as a build artifact) without re-running SPARTA.
"""
import os
import sys
import json
import argparse
import statistics as st
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from ad_stochastic_equivalence import (  # noqa: E402
    make_deck, run_one, transpose, select_tracked_columns, REPO_ROOT, Z_THRESHOLD,
    DEFAULT_TIMEOUT)

STOCK_B_OFFSET = 5_000_000  # keeps stock_B's seeds far from stock_A/ad_A's


def run_pool(binary, case_dir, base_deck, seeds, jobs, timeout=DEFAULT_TIMEOUT):
    """Run `binary` once per seed in `seeds`, in parallel. Returns a list of
    per-seed {quantity: value} dicts, in the same order as `seeds`."""
    import tempfile

    aux_files = [os.path.join(case_dir, f) for f in os.listdir(case_dir)
                 if not f.startswith("in.") and not f.startswith("log.")]

    def one(seed):
        with tempfile.TemporaryDirectory(prefix="ad_sweep_") as workdir:
            deck_text = make_deck(base_deck, seed)
            deck_path = os.path.join(workdir, "in.case")
            with open(deck_path, "w") as f:
                f.write(deck_text)
            return run_one(binary, deck_path, workdir, aux_files, timeout=timeout)

    with ThreadPoolExecutor(max_workers=jobs) as ex:
        return list(ex.map(one, seeds))


def mean_of(pool_slice, q):
    return st.mean(r[q] for r in pool_slice)


def sem_of(pool_slice, q):
    vals = [r[q] for r in pool_slice]
    if len(vals) < 2:
        return float("nan")
    return st.stdev(vals) / (len(vals) ** 0.5)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sweep-ns", default="4,8,16,32,64",
                     help="comma-separated, ascending sample sizes to check")
    ap.add_argument("--case", default="examples/bfield/in.bfield")
    ap.add_argument("--seed-base", type=int, default=9990001)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--out-json", default=None,
                     help="write full results (raw pools + per-N summary) here")
    ap.add_argument("--columns", default=None,
                     help="comma-separated column names to sweep (default: "
                          "auto-detect every non-constant column in the stats output)")
    ap.add_argument("--timeout", type=int, default=DEFAULT_TIMEOUT,
                     help=f"per-run subprocess timeout in seconds (default {DEFAULT_TIMEOUT})")
    args = ap.parse_args()

    sweep_ns = [int(x) for x in args.sweep_ns.split(",")]
    n_max = max(sweep_ns)

    stock = os.environ.get("SPARTA_STOCK")
    ad = os.environ.get("SPARTA_AD")
    if not stock or not ad:
        print("ERROR: set SPARTA_STOCK and SPARTA_AD to the two binaries to compare")
        return 1

    case_path = os.path.join(REPO_ROOT, args.case)
    case_dir = os.path.dirname(case_path)
    with open(case_path) as f:
        base_deck = f.read()

    seeds_a = [args.seed_base + i for i in range(n_max)]
    seeds_b = [args.seed_base + STOCK_B_OFFSET + i for i in range(n_max)]

    print(f"case={args.case}  sweep_ns={sweep_ns}  jobs={args.jobs}  timeout={args.timeout}s")
    print(f"running {n_max} stock_A seeds...")
    stock_a = run_pool(stock, case_dir, base_deck, seeds_a, args.jobs, args.timeout)
    print(f"running {n_max} stock_B seeds (noise-floor control)...")
    stock_b = run_pool(stock, case_dir, base_deck, seeds_b, args.jobs, args.timeout)
    print(f"running {n_max} AD seeds...")
    ad_a = run_pool(ad, case_dir, base_deck, seeds_a, args.jobs, args.timeout)

    requested = args.columns.split(",") if args.columns else None
    tracked = select_tracked_columns(transpose(stock_a), requested)
    print(f"tracked columns: {tracked}")

    summary = {q: [] for q in tracked}
    print()
    for q in tracked:
        print(f"=== {q} ===")
        print(f"{'N':>5} {'noise_floor':>14} {'ad_bias':>14} {'|bias|/|noise|':>16} "
              f"{'sem_noise':>12} {'sem_bias':>12}")
        for n in sweep_ns:
            noise = mean_of(stock_b[:n], q) - mean_of(stock_a[:n], q)
            bias = mean_of(ad_a[:n], q) - mean_of(stock_a[:n], q)
            sem_noise = (sem_of(stock_a[:n], q) ** 2 + sem_of(stock_b[:n], q) ** 2) ** 0.5
            sem_bias = (sem_of(stock_a[:n], q) ** 2 + sem_of(ad_a[:n], q) ** 2) ** 0.5
            ratio = abs(bias) / abs(noise) if noise != 0 else float("inf")
            print(f"{n:>5} {noise:>14.4f} {bias:>14.4f} {ratio:>16.2f} "
                  f"{sem_noise:>12.4f} {sem_bias:>12.4f}")
            summary[q].append({"N": n, "noise_floor": noise, "ad_bias": bias,
                                "sem_noise": sem_noise, "sem_bias": sem_bias})
        print()

    # Pass/fail: a z-test at the LARGEST N in the sweep (the most
    # statistically powerful point, since sem shrinks as 1/sqrt(N)) --
    # |ad_bias| must stay within Z_THRESHOLD sem_bias of zero, same
    # threshold/logic as ad_stochastic_equivalence.py's single-N check.
    # This is what actually gives the sweep a pass/fail verdict distinct
    # from being purely a printed diagnostic: a bias that FAILS this at
    # the largest N tested, after having had the most data to shrink into
    # noise, is the flattening-instead-of-shrinking signature of a
    # systematic bug, not sampling noise.
    ok = True
    print(f"{'quantity':<10} {'N':>5} {'ad_bias':>12} {'sem_bias':>12} {'z':>8}  result")
    for q in tracked:
        last = summary[q][-1]
        z = last["ad_bias"] / last["sem_bias"] if last["sem_bias"] > 0 else float("inf")
        good = abs(z) < Z_THRESHOLD
        ok = ok and good
        print(f"{q:<10} {last['N']:>5} {last['ad_bias']:>12.4f} {last['sem_bias']:>12.4f} "
              f"{z:>8.2f}  {'PASS' if good else 'FAIL'}")
    print(f"\nconvergence sweep (z-test at N={sweep_ns[-1]}, |z|<{Z_THRESHOLD}): "
          f"{'ALL PASS' if ok else 'FAIL'}")

    if args.out_json:
        out = {
            "case": args.case,
            "sweep_ns": sweep_ns,
            "seed_base": args.seed_base,
            "tracked": tracked,
            "raw": {"stock_A": stock_a, "stock_B": stock_b, "AD_A": ad_a},
            "summary": summary,
        }
        with open(args.out_json, "w") as f:
            json.dump(out, f, indent=2)
        print(f"wrote {args.out_json}")

    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
