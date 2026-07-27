#!/usr/bin/env python3
"""
ctest_fallback_gate.py -- statistical fallback gate for the AD build's
ctest regression suite (.github/workflows/ad.yml, job mpi-stubs-ad-sacado).

WHY THIS EXISTS
  A handful of gold-log regression inputs are known to sometimes diverge
  from a single fixed-seed trajectory under -DSPARTA_ENABLE_AD, not
  because the AD build is buggy, but because Sacado's expression-template
  arithmetic is not guaranteed bit-identical to plain double on the same
  inputs (confirmed semantic, not a compiler/codegen artifact -- see the
  comment block in cmake/common/set/sparta_cmake_defaults.cmake for the
  traced mechanism). Given enough particles and steps, that pervasive
  1-ULP noise eventually flips a probabilistic accept/reject branch or a
  truncation-based decision boundary, desyncing the RNG stream for the
  rest of that run.

  Rather than maintain a hand-curated, must-be-kept-in-sync list of which
  specific inputs are known to hit this (this repo's previous design),
  ctest now runs the FULL suite unconditionally, and this script provides
  an automatic, evidence-based fallback for whatever fails: if a failing
  input's AD build is still statistically equivalent to stock's across
  many seeds (unbiased mean, and that agreement holding as sample size
  grows), the divergence is ordinary AD-vs-double floating-point noise,
  not a real bug, and the job should not fail because of it. If the
  fallback ALSO fails, that's treated as a real failure.

METHOD
  Reads <build-dir>/Testing/Temporary/LastTestsFailed.log (a standard
  CTest artifact populated after `ctest` runs with at least one failure),
  extracts the deck name ("in.XXXX") from each failing test's name
  (stripping the machine prefix and optional ".mpi_N" suffix), de-
  duplicates across MPI-rank variants of the same deck, locates each
  deck's file under examples/, and runs BOTH
  ad_stochastic_equivalence.py (TOST-style mean-equivalence z-test) and
  ad_convergence_sweep.py (1/sqrt(N) convergence trend) against it.
  BOTH must pass for a deck to be treated as noise.

  Deliberately does NOT use ad_distribution_test.py here: running many
  independent hypothesis tests (one per tracked column, one per deck)
  makes occasional false positives from that test likely by chance alone
  (confirmed in practice -- see docs/PLAN.md), and it lacks the kind of
  physically-grounded equivalence margin ad_stochastic_equivalence.py's
  TOST framing has. It remains a useful, informative-only diagnostic --
  see .github/workflows/ad_distribution.yml, which reports it without
  ever failing a build on it -- but not a fallback gate.

  Fails safe: a ctest failure whose deck can't be identified, or whose
  file can't be found under examples/, is treated as a real failure (not
  silently passed through), since the whole point is to only wave through
  failures that are POSITIVELY confirmed to be AD-inherent noise.

USAGE
  SPARTA_STOCK=/path/to/stock/spa_<machine> \\
  SPARTA_AD=/path/to/ad/spa_<machine> \\
  python3 ctest_fallback_gate.py --build-dir build \\
      [--seeds 32] [--sweep-ns 4,8,16,32,64]

  Exit 0 if there were no ctest failures, or every failing deck's
  fallback passed both checks. Exit 1 otherwise.
"""
import os
import re
import sys
import argparse
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)
from ad_stochastic_equivalence import DEFAULT_TIMEOUT  # noqa: E402

# Matches sparta_add_test's naming convention (cmake/common/test/
# sparta_test_utils.cmake): "<machine>.<deck>[.mpi_<N>][<config>]", e.g.
# "serial.in.thermostat.mpi_1" or "serial.in.circle.constant.mpi_4".
# Non-greedy capture on the deck name so decks with dots in their own
# name (like "in.circle.constant") are captured whole, not truncated at
# the first dot.
TEST_NAME_RE = re.compile(r"^\w+\.(in\.\S+?)(?:\.mpi_\d+)?$")


def parse_failed_tests(build_dir):
    log_path = os.path.join(build_dir, "Testing", "Temporary", "LastTestsFailed.log")
    if not os.path.exists(log_path):
        return []
    names = []
    with open(log_path) as f:
        for line in f:
            line = line.strip()
            if not line or ":" not in line:
                continue
            names.append(line.split(":", 1)[1])
    return names


def deck_name_from_test(test_name):
    m = TEST_NAME_RE.match(test_name)
    return m.group(1) if m else None


def find_deck_path(deck_name):
    for root, _dirs, files in os.walk(os.path.join(REPO_ROOT, "examples")):
        if deck_name in files:
            return os.path.relpath(os.path.join(root, deck_name), REPO_ROOT)
    return None


def run_fallback(case_path, seeds, sweep_ns, timeout):
    """Streams both scripts' output directly to the CI log (not captured)
    so the underlying statistical detail is visible, not hidden behind a
    summary."""
    print(f"  -- ad_stochastic_equivalence.py --seeds {seeds} --case {case_path} "
          f"--timeout {timeout}")
    r1 = subprocess.run(
        [sys.executable, os.path.join(HERE, "ad_stochastic_equivalence.py"),
         "--seeds", str(seeds), "--case", case_path, "--timeout", str(timeout)],
        cwd=REPO_ROOT)
    equiv_ok = r1.returncode == 0

    print(f"  -- ad_convergence_sweep.py --sweep-ns {sweep_ns} --case {case_path} "
          f"--timeout {timeout}")
    r2 = subprocess.run(
        [sys.executable, os.path.join(HERE, "ad_convergence_sweep.py"),
         "--sweep-ns", sweep_ns, "--jobs", "4", "--case", case_path,
         "--timeout", str(timeout)],
        cwd=REPO_ROOT)
    sweep_ok = r2.returncode == 0

    return equiv_ok, sweep_ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--build-dir", default="build")
    ap.add_argument("--seeds", type=int, default=32)
    ap.add_argument("--sweep-ns", default="4,8,16,32,64")
    ap.add_argument("--timeout", type=int, default=DEFAULT_TIMEOUT,
                     help=f"per-run subprocess timeout in seconds, passed through to "
                          f"both fallback scripts (default {DEFAULT_TIMEOUT})")
    args = ap.parse_args()

    if not os.environ.get("SPARTA_STOCK") or not os.environ.get("SPARTA_AD"):
        print("ERROR: set SPARTA_STOCK and SPARTA_AD to the two binaries to compare")
        return 1

    failed_tests = parse_failed_tests(args.build_dir)
    if not failed_tests:
        print("no ctest failures recorded -- nothing for the fallback gate to check")
        return 0

    print(f"ctest reported {len(failed_tests)} failing test(s): {failed_tests}")

    decks, unmapped = {}, []
    for name in failed_tests:
        deck = deck_name_from_test(name)
        (unmapped.append(name) if deck is None
         else decks.setdefault(deck, []).append(name))

    if unmapped:
        print(f"\nERROR: could not parse a deck name out of: {unmapped} "
              f"-- failing safe (not silently passing an unrecognized failure)")
        return 1

    print(f"\n{len(decks)} unique deck(s) to check: {sorted(decks)}\n")

    all_ok = True
    for deck, test_names in sorted(decks.items()):
        case_path = find_deck_path(deck)
        print(f"=== {deck} (from {test_names}) ===")
        if case_path is None:
            print("  ERROR: could not find this deck under examples/ -- failing safe")
            all_ok = False
            continue
        print(f"  case={case_path}")
        equiv_ok, sweep_ok = run_fallback(case_path, args.seeds, args.sweep_ns, args.timeout)
        ok = equiv_ok and sweep_ok
        all_ok = all_ok and ok
        print(f"  mean-equivalence: {'PASS' if equiv_ok else 'FAIL'}   "
              f"convergence-sweep: {'PASS' if sweep_ok else 'FAIL'}   -> "
              f"{'fallback PASS (AD-inherent noise, not a bug)' if ok else 'fallback FAIL (real failure)'}\n")

    print(f"ctest fallback gate: {'ALL PASS' if all_ok else 'FAIL'}")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
