#!/usr/bin/env python3
"""
ad_stochastic_equivalence.py -- statistical equivalence check for the AD
build's stochastic (RNG-branch-sensitive) behavior.

WHY THIS EXISTS
  A handful of gold-log regression tests are EXCLUDED from the AD build's
  ctest suite (see cmake/common/set/sparta_cmake_defaults.cmake, the
  SPARTA_ENABLE_AD block) because they were proven -- by cross-checking
  against an independent, previously-validated hand-rolled AD engine -- to
  diverge from a single fixed-seed trajectory for ANY AD scalar
  substitution, not just this fork's Sacado port. Root cause: a rare
  (~1-in-tens-of-thousands) 1-ULP rounding difference in mixed
  double/sfloat arithmetic occasionally flips a probabilistic accept/reject
  or randomized-rounding branch (e.g. collide_vss.cpp's
  "attempt += random->uniform()", or create_particles.cpp's per-cell
  particle-count accept test), desyncing the RNG stream for the rest of
  that run. Traced directly in a live run (docs/PLAN.md has the full
  trace), not assumed.

  A SINGLE fixed-seed trajectory can therefore never be relied on to match
  between stock and AD for these tests. But "different specific realization
  of the same random process" is a testable claim distinct from "silently
  wrong": if the AD build is correct, the DISTRIBUTION of outcomes across
  many seeds must be statistically indistinguishable from stock's
  distribution across the same seeds -- same mean, no systematic bias, just
  ordinary seed-to-seed scatter. A real bug (wrong probability formula,
  dropped energy term, sign error) would show up as a BIASED mean, not just
  scatter. This script checks exactly that, replacing the single-trajectory
  regression test these excluded cases can't use.

METHOD
  Run the same input (default: examples/bfield/in.bfield) N times each on
  a stock and an AD binary, varying only the seed. Every stats-line column
  except "Step" (loop index) and "CPU" (wall-clock timing) is a candidate
  tracked quantity -- which ones actually get tested is decided per case,
  not hardcoded, since different decks report different columns (e.g.
  in.bfield has Natt/Ncoll/c_temp, in.shocktube has Natt/Ncoll/Nscoll/
  Nscheck with no c_temp at all, in.surf_react_heatflux has Nreact/
  Nsreact/Ngrid/c_echem/Maxlevel with no Natt). See select_tracked_columns():
  a column is dropped if it's constant across the reference (stock)
  ensemble -- e.g. Np in in.bfield never varies (no particle
  insertion/removal in that deck), so there is nothing to test. Pass
  --columns to test an explicit subset instead of the auto-detected set.

  For each tracked quantity, compute:
    - mean and stdev of each ensemble
    - pooled standard error of the mean difference (SEM_pooled =
      sqrt(SEM_stock^2 + SEM_ad^2))
    - z = (mean_ad - mean_stock) / SEM_pooled
  PASS if |z| < Z_THRESHOLD for every tracked quantity (default 3.0, i.e.
  the AD ensemble mean must fall within stock's own expected seed-to-seed
  sampling noise -- a systematic bug would push |z| far past this).
  This is NOT a claim of proof; it is a statistical screen. Increasing
  N_SEEDS tightens it (SEM shrinks as 1/sqrt(N)).

USAGE
  SPARTA_STOCK=/path/to/stock/spa_<machine> \\
  SPARTA_AD=/path/to/ad/spa_<machine> \\
  python3 ad_stochastic_equivalence.py [--seeds N] [--case examples/bfield/in.bfield]

WORKFLOW INTENT (not yet wired into CI)
  Standard library only. Needs TWO binaries built in the same job (a stock
  build and an -DSPARTA_ENABLE_AD build), unlike the single-binary jobs in
  ad.yml today -- that's the reason this isn't wired in yet, not a
  limitation of the test itself. Runtime scales with --seeds (default 16
  seeds x 2 builds; each bfield run is a few seconds).
"""
import os
import re
import sys
import math
import shutil
import argparse
import tempfile
import subprocess
import statistics as st

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))

Z_THRESHOLD = 3.0
IGNORE_COLUMNS = {"Step", "CPU"}  # never physics: loop index and wall-clock timing

SEED_RE = re.compile(r"^(\s*seed\s+)\S+", re.M)


def make_deck(base_deck_text, seed):
    return SEED_RE.sub(lambda m: f"{m.group(1)}{seed}", base_deck_text, count=1)


def parse_stats_line(header, toks):
    """Turn one "Step CPU Np Natt Ncoll c_temp ..." stats row into
    {column_name: float}, skipping IGNORE_COLUMNS and any non-numeric
    column. Keeps the deck's own column names (e.g. "c_temp", "Natt") --
    no renaming -- so results are directly comparable to the gold log."""
    vals = {}
    for name, tok in zip(header, toks):
        if name in IGNORE_COLUMNS:
            continue
        try:
            vals[name] = float(tok)
        except ValueError:
            continue
    return vals


def run_one(binary, deck_path, workdir, aux_files):
    for f in aux_files:
        shutil.copy(f, workdir)
    out = subprocess.run([binary, "-in", os.path.basename(deck_path)],
                         cwd=workdir, capture_output=True, text=True, timeout=120)
    if out.returncode != 0:
        raise RuntimeError(f"SPARTA failed (exit {out.returncode}):\n{out.stdout[-1500:]}")
    # Read the LAST stats row of the LAST "run" block -- not a fixed step
    # number, since different decks run for different lengths, and some
    # decks issue multiple "run" commands (each prints its own "Step ..."
    # header followed by its own stats rows); tracking the most recent
    # header/row pair as we scan naturally lands on the final state of the
    # whole simulation, whichever block it came from. Header is whatever
    # this deck's "stats" command configured, e.g. "Step CPU Np Natt Ncoll
    # c_temp" -- different decks configure different columns.
    header = None
    last_header = None
    last_row = None
    for line in out.stdout.splitlines():
        toks = line.split()
        if toks[:1] == ["Step"]:
            header = toks
            continue
        if header and toks and toks[0].lstrip("-").isdigit():
            last_header, last_row = header, toks
    if last_row is None:
        raise RuntimeError("no stats rows found in SPARTA output")
    return parse_stats_line(last_header, last_row)


def transpose(pool):
    """[{col: val, ...}, ...] (one dict per seed) -> {col: [val, ...]}."""
    if not pool:
        return {}
    return {col: [r[col] for r in pool] for col in pool[0]}


def select_tracked_columns(reference_results, requested=None):
    """Decide which columns are actually worth testing. `reference_results`
    is a {column: [values]} dict (typically from the stock ensemble).
    A column with zero variance in the reference (e.g. Np in in.bfield,
    which never changes -- no particle insertion/removal in that deck) has
    nothing to test and is dropped, with a printed note, unless it was
    explicitly `requested`. Raises if a `requested` column doesn't exist
    in the stats output at all."""
    available = list(reference_results.keys())
    if requested:
        missing = [c for c in requested if c not in available]
        if missing:
            raise SystemExit(f"requested column(s) not found in stats output: {missing} "
                              f"(available: {sorted(available)})")
        return list(requested)
    tracked, dropped = [], []
    for c in available:
        vals = reference_results[c]
        (dropped if len(set(vals)) <= 1 else tracked).append(c)
    if dropped:
        print(f"dropping constant column(s) (no variance to test): {sorted(dropped)}")
    return tracked


def ensemble(binary, case_dir, base_deck, n_seeds, seed_base):
    pool = []
    aux_files = [os.path.join(case_dir, f) for f in os.listdir(case_dir)
                if not f.startswith("in.") and not f.startswith("log.")]
    with tempfile.TemporaryDirectory(prefix="ad_stoch_") as workdir:
        for i in range(n_seeds):
            seed = seed_base + i
            deck_text = make_deck(base_deck, seed)
            deck_path = os.path.join(workdir, "in.case")
            with open(deck_path, "w") as f:
                f.write(deck_text)
            pool.append(run_one(binary, deck_path, workdir, aux_files))
    return transpose(pool)


def z_score(a, b):
    ma, mb = st.mean(a), st.mean(b)
    sa, sb = st.stdev(a), st.stdev(b)
    sem = math.sqrt((sa / math.sqrt(len(a))) ** 2 + (sb / math.sqrt(len(b))) ** 2)
    return (mb - ma) / sem if sem > 0 else float("inf"), ma, mb, sa, sb


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=16)
    ap.add_argument("--case", default="examples/bfield/in.bfield")
    ap.add_argument("--seed-base", type=int, default=8880001)
    ap.add_argument("--columns", default=None,
                     help="comma-separated column names to test (default: "
                          "auto-detect every non-constant column in the stats output)")
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

    print(f"case={args.case}  seeds={args.seeds}  reading last stats row of each run")
    print("running stock ensemble...")
    stock_r = ensemble(stock, case_dir, base_deck, args.seeds, args.seed_base)
    print("running AD ensemble...")
    ad_r = ensemble(ad, case_dir, base_deck, args.seeds, args.seed_base)

    requested = args.columns.split(",") if args.columns else None
    tracked = select_tracked_columns(stock_r, requested)
    print(f"tracked columns: {tracked}")

    ok = True
    print(f"\n{'quantity':<10} {'stock mean':>12} {'AD mean':>12} {'stock std':>10} "
          f"{'AD std':>10} {'z':>8}  result")
    for q in tracked:
        z, ma, mb, sa, sb = z_score(stock_r[q], ad_r[q])
        good = abs(z) < Z_THRESHOLD
        ok = ok and good
        print(f"{q:<10} {ma:12.4f} {mb:12.4f} {sa:10.4f} {sb:10.4f} {z:8.2f}  "
              f"{'PASS' if good else 'FAIL'}")

    print(f"\nstochastic equivalence ({args.seeds} seeds, |z|<{Z_THRESHOLD}): "
          f"{'ALL PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
