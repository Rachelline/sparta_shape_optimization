# Known gaps

Things that are deliberately unfinished, untested, or deferred — recorded
here so they're discoverable without reading `git log`. None of these are
silent: each one also shows up as a code comment at the relevant spot, a
runtime warning, or both.

## AD gradient correction

- **`nfactor_inverse`-weighted tallies (force/torque) are not corrected.**
  `--score-correction` only reaches `fluxscale`-weighted tallies
  (press/shear/heat-flux — `HeatFluxObjective`). It's a documented no-op
  for `DragObjective` (`fx`, scaled by `nfactor_inverse`). See
  `docs/AD_GRADIENTS.md`'s coverage table. Extending the correction here
  would need a `src/compute_surf.cpp` change (a `src/` touch, not
  `ad_src/`) — deferred as separate, later work by explicit decision, not
  overlooked.
- **Diffuse-wall reflection resampling and multi-collision compounding are
  not analyzed.** The `u/spval(u)` correction's derivation and validation
  (`docs/AD_GRADIENTS.md`) cover specular reflection's flux-measure bias
  specifically. Diffuse walls (an `|u|`-independent re-emission term
  breaks the exact-half symmetry specular reflection has) and collisional
  (non-free-molecular) runs (repeated resampling across `collide vss`
  events) may carry additional, uncharacterized bias beyond what
  `--score-correction` addresses. Not verified either way.

## ad_src shape-optimization framework

- **`DragObjective` on a 3D shape (`--shape powerlaw --objective drag`)
  is untested.** `fx` is a dimension-agnostic compute-surf keyword and
  should work mechanically, but no verification run has exercised this
  combination — every 3D verification in this codebase used
  `--objective heatflux`.
- **`shape_main`/`shape_main_ad` don't offer `--shape powerlaw`.** Only
  `shape_opt`/`shape_opt_ad` do. A single-point 3D evaluation currently
  has to go through the optimizer driver with `--max-iter 0`.
- **`MinSizeConstraint`'s 3D volume branch has no SPARTA-level
  verification**, only a standalone finite-difference self-test (no
  SPARTA involved) confirming the analytic gradient is correct
  (`grad` vs FD: ~4e-13). It has never been exercised through an actual
  `shape_opt --shape powerlaw --min-area N` optimization run.
- **`svg_shape.cpp` is hardcoded to `alpha[4]` / a 2D Bezier outline.** A
  second 2D shape family with a different `ndesign()` would need its own
  renderer; the current `shapes.svg` gate (`shape->dim() == 2`) doesn't
  protect against that case, only against 3D shapes.
- **No `--pl-*` flags exist on `shape_main` for the power-law body's own
  geometry parameters** (`L`, `Rmax`, `nx`, `ntheta`) — those are only
  configurable on `shape_opt`/`shape_opt_ad`.

## SPARTA AD build (not ad_src-specific)

- **Sacado's expression-template arithmetic is not guaranteed
  bit-identical to plain `double`**, even on identical inputs (confirmed
  semantic — not a compiler/codegen artifact; see
  `cmake/common/set/sparta_cmake_defaults.cmake`'s comment block for the
  traced mechanism, and `tools/ad_verify/`'s stochastic-equivalence
  scripts for the CI-level workaround). A handful of gold-log regression
  tests are excluded from the AD build's ctest suite for this reason and
  checked statistically instead. This is inherent to Sacado's
  implementation, not something this codebase can fix.

## Documentation

- **`docs/PLAN.md` and `docs/ad_phase_c_investigation/FINDINGS.md` are
  gitignored working notes**, not committed documentation — by explicit
  decision, their content was distilled into `docs/AD_GRADIENTS.md`
  (committed) rather than committing the originals. Anyone without access
  to those local files is not missing anything `AD_GRADIENTS.md` doesn't
  already cover for the AD-gradient-bias topic specifically; they do
  contain additional historical investigation detail (disproven
  hypotheses, exact numeric tables) not reproduced here.
