# SPARTA Shape Optimization

An AD-linked DSMC shape-optimization framework built on top of
[SPARTA](https://sparta.github.io/): pick a parametrized body shape, an
objective (drag or heat flux), optionally a size constraint, and either
finite-difference or Sacado forward-mode AD gradients feed an IPOPT
optimizer loop.

- **[Overview](ad_src/)** — architecture, capability matrix, the full AD
  build recipe, how to run each target, how to flip the AD
  score-correction switch, and extension HOWTOs.
- **[AD Gradients](docs/AD_GRADIENTS/)** — why forward-mode AD is biased
  on DSMC surface tallies, the `u/spval(u)` correction mechanism, and a
  coverage table of what's fixed and what isn't.
- **[Verification Scripts](tools/ad_verify/)** — index of the AD
  correctness verification scripts, what each proves, and which run in
  CI.
- **[Build (CMake)](BUILD_CMAKE/)** — the top-level SPARTA CMake build,
  including the automatic-differentiation option reference.
- **[Known Gaps](KNOWN_GAPS/)** — every deliberately-unfinished or
  untested piece, in one place.
