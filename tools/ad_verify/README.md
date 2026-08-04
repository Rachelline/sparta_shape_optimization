# tools/ad_verify: AD correctness verification scripts

An index of what each script proves, what it needs, and whether CI runs
it. Most of these predate `ad_src/` and validate SPARTA's AD build
directly (raw `spa_<machine>` binaries via `SPARTA_STOCK`/`SPARTA_AD` env
vars); `compare_modes.py` is the exception, built for `ad_src/`'s own
drivers specifically.

For the theory these scripts are validating (the flux-measure derivative
bias and what corrects it), see
[`docs/AD_GRADIENTS.md`](../../docs/AD_GRADIENTS.md).

## Env vars

Most scripts read `SPARTA_STOCK`/`SPARTA_AD` (paths to a stock and an
AD-enabled `spa_<machine>` binary); `flatplate_value_match.py` and
`fmf_coeff.py` only need `SPARTA` (a single stock binary — no AD
involved). `compare_modes.py` is different again: it drives
`ad_src/build/{shape_main,shape_main_ad}` by path (defaults relative to
this repo, override with `--stock-bin`/`--ad-bin`), not raw SPARTA
binaries.

```bash
export SPARTA_STOCK=/path/to/stock/build/src/spa_serial
export SPARTA_AD=/path/to/ad/build/src/spa_ad
```

## In CI

| Script | Workflow |
|---|---|
| `fmf_coeff.py` (self-test) | `.github/workflows/ad.yml` |
| `flatplate_value_match.py` | `.github/workflows/ad.yml` |
| `ctest_fallback_gate.py` | `.github/workflows/ad.yml` |
| `ad_stochastic_equivalence.py` | `.github/workflows/ad_distribution.yml` |
| `ad_convergence_sweep.py` | `.github/workflows/ad_distribution.yml` |
| `ad_distribution_test.py` | `.github/workflows/ad_distribution.yml` |

Everything else below is a manual/investigation tool, not wired into CI.

## Analytic reference (Phase A)

- **`fmf_coeff.py`** — closed-form free-molecular-flow (FMF) aerodynamic
  coefficients for a flat panel: a zero-noise analytic ground truth for
  both the drag *value* and, by differentiation, `d(drag)/d(design)` —
  something a noisy DSMC run can't provide directly. Self-verified, no
  SPARTA dependency. Everything below compares against this.

## Value-match gates (Phase B — no AD needed)

- **`flatplate_value_match.py`** — stock SPARTA reproduces the closed-form
  FMF panel pressure at every incidence angle; specular ⇒ shear ≈ 0. Pins
  units, `q_inf` normalization, speed-ratio and angle conventions before
  any derivative comparison is meaningful.
- **`thinplate_value_match.py`** — same idea on a small closed
  thin-rectangle plate, reading a single dumped surf element's press
  (not a reduce-sum across all four sides, which would mix physically
  different flow conditions).

## Derivative-match ladder (Phase C/D — needs AD)

- **`flatplate_derivative_match.py`** — seeds a Sacado direction into the
  freestream velocity (`d(vx,vy)/d(delta)`), so `d(pressure)/d(delta)`
  propagates through the full emission → free-flight → reflection →
  tally pipeline automatically; compares against `fmf_coeff.py`'s
  closed-form derivative.
- **`thinplate_derivative_match.py`** — same comparison, but seeds the
  *plate's own point coordinates* via `src/read_surf.cpp`'s
  `SPARTA_AD_SEED_ALPHA` rotation hook instead of the freestream —
  sidesteps `FINDING 1` (an emission-rate truncation bug, since fixed)
  entirely, since the freestream is never seeded.
- **`measure_theory_check.py`** — **the decisive self-test** for the
  flux-measure bias (`docs/AD_GRADIENTS.md`). No SPARTA needed, runs in
  seconds: an independent Monte Carlo sample from a drifting Maxwellian
  pins the exact-factor-of-2 ratio for specular pressure to four decimal
  places.
- **`hitcount_sensitivity.py`** — measures `d(hit count)/d(alpha)` on the
  stock binary (no AD) via `compute surf`'s raw `NUM` tally; supplies an
  independent validation target for the missing score-function term
  (the "capture-count" half of the derivative).
- **`synthetic_single_particle.py`** — deterministic, zero-RNG test (one
  hand-placed particle) isolating whether the AD *machinery* itself
  (reflection formula, normal computation, tally) is correct, independent
  of any statistical/grid-topology question. Passing here does NOT imply
  end-to-end correctness — it structurally cannot see the flux-measure
  bias, since fixing the sampling measure by construction is exactly what
  removes the bug it's blind to.
- **`multiseed_ad_test.py`** — averages `dP/dalpha` across independent
  seeds at one fixed angle, to distinguish a structural bias (mean
  converges tightly to the *same* offset from FD/analytic) from
  statistical noise (mean would converge toward FD/analytic given enough
  seeds).

## Stochastic-equivalence tests (AD build's RNG-branch-sensitive behavior)

A handful of gold-log regression tests are excluded from the AD build's
ctest suite because Sacado's expression-template arithmetic is not
guaranteed bit-identical to plain `double` on the same inputs (a real,
semantic difference — confirmed not a compiler/codegen artifact; see the
comment block in `cmake/common/set/sparta_cmake_defaults.cmake`). That
occasional 1-ULP divergence can flip a probabilistic accept/reject branch,
desyncing the RNG stream for the rest of a run. These three scripts
replace exact-match gold-log comparison with statistical equivalence for
exactly those excluded tests:

- **`ad_stochastic_equivalence.py`** — TOST-style mean-equivalence check
  at a fixed sample size: is the AD-vs-stock mean difference within noise?
- **`ad_convergence_sweep.py`** — sweeps sample size `N` and checks the
  AD-vs-stock mean difference shrinks like `1/sqrt(N)` (sampling noise)
  rather than converging to a nonzero constant (a systematic bias a
  single fixed-`N` check could miss).
- **`ad_distribution_test.py`** — permutation Kolmogorov-Smirnov test
  comparing the full empirical distributions, not just their means —
  catches a bug that leaves the mean roughly unchanged but distorts the
  distribution's shape (e.g. a heavier tail from an occasional bad
  cascade after an RNG desync).
- **`ctest_fallback_gate.py`** — the CI-level glue: runs ctest, and for
  any test in the known-divergent set that fails exact gold-log
  comparison, falls back to the statistical checks above instead of
  failing the build outright.

## Sacado API validation

- **`sacado_seed_selftest.cpp`** — standalone (not a CMake target), links
  only against Sacado. Validates the exact derivative-seeding recipe used
  in `src/mixture.cpp`/`src/compute_boundary.cpp`/`src/read_surf.cpp`
  against this codebase's actual Sacado storage policy
  (`Sacado::Fad::SFad<double,N>` / `StaticFixedStorage`), specifically
  checking whether a documented landmine ("a default-constructed `SFad`
  has `size()==0`, so poking `.fastAccessDx(j)` doesn't seed it")
  actually applies here. Empirically it does not — `StaticFixedStorage`'s
  derivative array is fixed-size at compile time — but this is verified,
  not assumed.

## ad_src-specific

- **`compare_modes.py`** — runs one `ad_src` shape/objective/alpha config
  three ways (FD via `shape_main`, AD-uncorrected and AD-corrected via
  `shape_main_ad`) and tabulates value, gradient, and the AD/FD ratio per
  component. The direct, scriptable answer to "run with and without the
  fix" — see `ad_src/README.md`'s "Flipping the AD correction switch"
  section.

## Supporting gas data

`N.species` (single-species gas data, reused by every script and by
`ad_src/`'s own drivers).
