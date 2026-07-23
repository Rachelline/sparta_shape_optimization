# ##############################################################################
# This file sets common default options that all sparta builds use. These
# options can be overridden at configure time via `cmake -DVAR=VAL` or `cmake -C
# /path/to/preset/presets.cmake`
# ##############################################################################
set(SPARTA_DEFAULT_CXX_COMPILE_FLAGS
    -DSPARTA_GZIP
    CACHE
      STRING
      "Compiler flags used when building object files for the \"spa_\" executable"
)

set(SPARTA_MACHINE
    ""
    CACHE
      STRING
      "Suffix to append to spa binary (WON'T enable any features automatically)"
)

if(SPARTA_ENABLE_TESTING)
  set(SPARTA_ENABLED_TEST_SUITES
      "ablation"
      "adapt"
      "vibrate"
      "surf_collide"
      "surf"
      "surf_react_adsorb"
      "step"
      "spiky"
      "sphere"
      "jagged"
      # FAILING."implicit"
      "free"
      "flowfile"
      "emit"
      "collide"
      "circle"
      "chem"
      "cylinder"
      "axi"
      "ambi"
      "relax_const"
      "relax_variable"
      "thermostat"
      "bfield"
      "adjust_temp"
      "shock_tube"
      "variable_timestep"
      "surf_react_heatflux"
      "chem_rates"
      "custom"
      "explicit2implicit"
      "mfp_mct"
      "torque")

  set(SPARTA_DISABLED_TESTS
      "in.ablation.3d.reactions" # Failing
      "in.axi" # Failing
      "in.collide" # Failing
      "in.ambi" # Failing
      "in.cylinder" # Long runtime
      "in.jagged.3d" # Long runtime
      "in.jagged.3d.distributed" # Long runtime
      "in.custom.cube.read.restart" # Failing
      "in.custom.cube.set.restart" # Failing
      "in.custom.step.read.restart" # Failing
      "in.custom.step.set.restart" # Failing
  )

  # When running the KOKKOS regression tests (SPARTA_KOKKOS_EXACT, run with
  # "-k on -sf kk"), skip the inputs that use features which are not yet
  # KOKKOS-enabled and would error out at run time. These tests still run in
  # the non-KOKKOS configurations.
  if(SPARTA_KOKKOS_EXACT)
    list(APPEND SPARTA_DISABLED_TESTS
        # fix ave/grid for grid/surf inputs not yet supported in KOKKOS
        "in.ablation.2d"
        "in.ablation.3d"
        # surf_collide adiabatic/cll/td/impulsive styles not KOKKOS-enabled
        "in.beam.adiabatic"
        "in.beam.cll"
        "in.beam.impulsive"
        "in.beam.td"
        "in.circle.adiabatic"
        "in.circle.cll"
        "in.circle.impulsive"
        "in.circle.td"
        # surf_react gs/ps styles use a non-KOKKOS-enabled surf_collide method
        "in.beam.face.gs"
        "in.beam.face.gs_ps"
        "in.beam.face.ps"
        "in.beam.surf.gs"
        "in.beam.surf.gs_ps"
        "in.beam.surf.ps"
        "in.circle.gs"
        "in.circle.gs_ps"
        "in.circle.ps"
        # external field fix not KOKKOS-enabled
        "in.bfield"
        "in.bfield.grid"
    )
  endif()

  # NOTE on -DSPARTA_ENABLE_AD and gold-log regression tests: no tests are
  # statically excluded here for AD. A handful of inputs are known to
  # sometimes diverge from a fixed-seed gold log under AD even though they
  # are not buggy -- see the CI-level explanation below, not a hardcoded
  # skip list here.
  #
  # Root cause (traced live via bit-exact in-binary shadow comparison on
  # in.bfield, 2026-07-22 -- see docs/PLAN.md for the full trace): Sacado's
  # Fad expression-template evaluation of even the simplest arithmetic
  # (e.g. update.cpp's ballistic move "xnew[d] = x[d] + dtremain*v[d]") is
  # NOT guaranteed bit-identical to the same formula evaluated in plain
  # double, even given identical input values. This is NOT rare: instrumenting
  # that one line alone showed ~940k one-ULP-or-greater disagreements over a
  # single 1000-step/10000-particle in.bfield run (roughly one in twenty
  # particle-dimension updates, every step). It is also NOT a compiler
  # codegen artifact -- confirmed two ways: (1) disabling FMA contraction
  # (-ffp-contract=off) on both stock and AD builds changed nothing about the
  # divergence; (2) recompiling the ballistic-move code at -O0 produced the
  # exact same disagreement count as the default optimized build. The
  # difference is semantic, coming from Sacado's own expression-template
  # evaluation path, not from instruction scheduling/fusion.
  #
  # Given that volume of 1-ULP noise, across 10000s of particles and 1000s
  # of steps it is essentially certain that eventually one such difference
  # lands close enough to a probabilistic accept/reject branch (e.g.
  # collide_vss.cpp's "attempt += random->uniform()" then truncated to int,
  # or test_collision's vre/vremax-vs-random->uniform() comparison) or a
  # truncation-based grid-cell-index boundary (update.cpp's
  # "static_cast<int>(spval((xnew[0]-boxlo[0])/dx))") to flip which side of
  # that boundary a value falls on, desyncing the RNG stream and the
  # per-particle cell trajectory for the rest of that run. This is expected
  # floating-point non-determinism from using a different arithmetic
  # evaluation engine for the same formula -- structurally the same
  # phenomenon that already makes plain-double SPARTA non-bit-reproducible
  # across compilers/optimization levels/architectures, just made visible
  # here because Sacado's evaluation path differs from plain double on
  # (almost) every operation rather than only rare corner cases. There is no
  # build flag or optimization-level change that removes it -- it cannot be
  # "fixed" without hand-recomputing every value in plain double alongside
  # Sacado's derivative tracking, which defeats the point of using Sacado.
  #
  # "Diverges from the fixed-seed gold log" is not the same claim as "wrong".
  # Rather than maintain a hand-curated, must-be-kept-in-sync list of which
  # specific inputs are known to hit this (the previous design here), CI
  # (.github/workflows/ad.yml, job mpi-stubs-ad-sacado) instead runs the
  # FULL ctest suite unconditionally and, for any test that fails, falls
  # back to tools/ad_verify/ad_stochastic_equivalence.py +
  # ad_convergence_sweep.py against that specific input: the stronger,
  # statistically meaningful claim that the AD build's DISTRIBUTION of
  # outcomes across many seeds is unbiased relative to stock's -- i.e. this
  # is an ordinary different draw of the same correct random process, not a
  # systematic error. Only if that fallback ALSO fails does the job report
  # a real failure. This is automatic and evidence-based per input, rather
  # than needing a human to notice a new divergence and add it to a list.
  # See tools/ad_verify/ctest_fallback_gate.py.

  list(APPEND __DEFAULT_MPI_RANKS "1")
  list(APPEND __DEFAULT_MPI_RANKS "4")
endif()
