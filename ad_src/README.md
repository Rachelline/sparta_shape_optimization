# ad_src: shape optimization on top of SPARTA

A standalone CMake project (not wired into `sparta/cmake` or
`sparta/src/CMakeLists.txt`) that runs DSMC shape optimization against SPARTA
through its public C library interface (`src/library.h`): pick a
parametrized body shape, an objective (drag or heat flux), optionally a
size constraint, and either finite-difference or Sacado forward-mode AD
gradients feed an IPOPT optimizer loop.

This file is the entry point for extending it. For the *theory* behind the
AD gradient bias and what `--score-correction` does and doesn't fix, see
[`docs/AD_GRADIENTS.md`](../docs/AD_GRADIENTS.md). For an index of the
verification scripts, see
[`tools/ad_verify/README.md`](../tools/ad_verify/README.md). For known gaps
and deferred work, see [`KNOWN_GAPS.md`](../KNOWN_GAPS.md).

## Architecture

Four interfaces, each with exactly one job:

```
Shape            geometry: alpha -> triangulated/polylined surface mesh
Objective        physics: SPARTA compute/fix commands -> a scalar QOI
Constraint       pure geometry: a size (area/volume) bound on alpha
shape_case.h     the runner: builds a Shape's mesh, drives SPARTA, calls
                 an Objective, does FD or reads back an AD gradient
```

```cpp
// shape.h -- one interface per shape family (2D or 3D)
class Shape {
  int dim() const;                                     // 2 or 3
  int ndesign() const;                                  // len(alpha)
  void bounds(double *lo, double *hi) const;
  void to_mesh(const double *alpha, SurfMesh &m) const;  // alpha -> mesh
  void jacobian(const double *alpha,
               std::vector<double> &jac) const;         // d(mesh)/d(alpha)
  bool validate(const SurfMesh &m, std::string *why) const;  // optional
};

// objective.h -- one interface per QOI
class Objective {
  void setup(void *spa, int navg) const;                // issue compute/fix
  double extract(void *spa, int ndesign, double *grad) const;  // read result
};

// constraint.h -- one interface per size bound
class Constraint {
  double eval(const Shape &shape, const double *alpha) const;
  void grad(const Shape &shape, const double *alpha, double *grad) const;
  double lower() const; double upper() const; const char *name() const;
};
```

`shape_case.h`'s `evaluate()`/`evaluate_avg()`/`grad_fd()` are the only
functions that ever touch SPARTA's library API directly — they're
dimension-aware (2D `Lines` vs 3D `Triangles` surf files, `dimension 2` vs
`3`, per-shape box/grid/boundary from `RunConfig`) but shape-family-agnostic:
adding a new `Shape` never requires touching `shape_case.cpp`.

`ShapeTNLP` (an IPOPT `TNLP`) wraps a `const Shape&`/`const Objective&` pair
plus an optional list of `const Constraint*`, caches value/gradient between
IPOPT's paired `eval_f`/`eval_grad_f` calls, and drives the actual
optimization loop.

### Files

| File | Role |
|---|---|
| `shape.h` | `Shape`/`SurfMesh` interface |
| `bezier_shape.h/.cpp`, `bezier_geom.h/.cpp` | 2D symmetric-Bezier-body `Shape` (adapter + pure geometry kernel) |
| `power_law_shape.h/.cpp`, `power_law_body.h/.cpp` | 3D power-law body-of-revolution `Shape` (adapter + pure geometry kernel) |
| `polygon2d.h/.cpp` | shared 2D polygon utilities (shoelace area + gradient, min edge length) |
| `objective.h`, `surf_sum_objective.h/.cpp` | `Objective` interface + shared "sum one compute-surf keyword" implementation |
| `drag_objective.h`, `heat_flux_objective.h` | thin named `SurfSumObjective` subclasses (`fx`/`nfactor_inverse`, `etot`/`fluxscale`) |
| `constraint.h`, `min_size_constraint.h/.cpp` | `Constraint` interface + generic area(2D)/volume(3D) floor |
| `ad_seed.h` | the AD *read* seam: `ad_extract()` (seeding lives in `src/read_surf.cpp`'s `SPARTA_AD_SEED_JACFILE` hook) |
| `shape_case.h/.cpp` | `RunConfig`, `AdCorrections`, `evaluate`/`evaluate_avg`/`grad_fd` |
| `sparta_util.h/.cpp` | thin SPARTA library wrappers (`cmd`, `open_sparta`/`close_sparta`, `extract_compute`, `die`) |
| `cli.h/.cpp` | shared strict CLI parsing (`parse_doubles`/`parse_ints`/`trim`) |
| `run_output.h/.cpp` | shared output-directory naming |
| `shape_main.cpp` | single-point evaluator CLI (`shape_main`/`shape_main_ad`) |
| `shape_opt_main.cpp` | IPOPT optimizer CLI (`shape_opt`/`shape_opt_ad`) |
| `shape_tnlp.h/.cpp` | the IPOPT `TNLP` implementation |
| `svg_shape.h/.cpp` | before/after SVG renderer (2D/Bezier-specific by design) |

## Capability matrix (honest, as of this writing)

| | 2D (`bezier`) | 3D (`powerlaw`) |
|---|---|---|
| Single-point eval (`shape_main`) | yes | no — only wired into `shape_opt` |
| IPOPT optimization (`shape_opt`) | yes | yes |
| `--objective drag` | yes, validated | untested (should work — `fx` is dimension-agnostic) |
| `--objective heatflux` | yes, validated (needs `--wall-accom > 0`) | yes, validated |
| `--min-area` size constraint | yes (shoelace area) | yes (divergence-theorem volume) |
| AD gradient (`_ad` binaries) | yes | yes |
| `--score-correction` effect | **no-op for drag** (fx routes through `nfactor_inverse`, not `fluxscale`); **live for heatflux** (`etot` routes through `fluxscale`) | same rule, same caveat |
| `shapes.svg` output | yes | no (renderer is Bezier-specific; gated on `dim()==2`) |

See [`docs/AD_GRADIENTS.md`](../docs/AD_GRADIENTS.md)'s coverage table for
exactly which tallies the AD correction reaches.

## Building

### Stock (finite-difference gradients only)

Needs a stock SPARTA static library, built automatically as a CMake custom
command (shells out to `make -C src mode=lib serial`, adding
`-DSPARTA_UNORDERED_MAP` — see the CMakeLists.txt comment for why this
fork's `Makefile.serial` needs that override). No SPARTA build step of your
own required.

```bash
cd ad_src
cmake -S . -B build
cmake --build build
./build/shape_main --alpha 1.3,1.0,2.7,0.8
```

IPOPT is optional (`brew install ipopt` on macOS; headers under
`coin-or/`). If not found, `shape_main`/`shape_main_ad` still build and
`shape_opt`/`shape_opt_ad` are skipped with a `message(WARNING ...)`.

### AD (Sacado forward-mode gradients)

This needs two things built *first*, outside `ad_src/`: a minimal
Trilinos/Sacado install, and a top-level SPARTA CMake build configured with
`-DSPARTA_ENABLE_AD=ON` linked against it. `ad_src/CMakeLists.txt` does not
build either of these itself — the top-level AD path has real complexity
(Kokkos ODR guard, C++17+ bump, Sacado discovery) not worth duplicating in a
nested custom command.

**1. Build a minimal Trilinos/Sacado.** No Homebrew formula exists for
Sacado alone; build Trilinos with everything but Sacado disabled:

```bash
git clone --depth 1 https://github.com/trilinos/Trilinos.git trilinos_src
cmake -S trilinos_src -B trilinos_build \
  -DCMAKE_INSTALL_PREFIX=/path/to/sacado_install \
  -DCMAKE_BUILD_TYPE=Release \
  -DTrilinos_ENABLE_Fortran=OFF \
  -DTrilinos_ENABLE_Sacado=ON \
  -DTrilinos_ENABLE_ALL_PACKAGES=OFF \
  -DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES=OFF \
  -DTPL_ENABLE_MPI=OFF \
  -DCMAKE_CXX_STANDARD=20            # current Trilinos main requires 20 or 23
cmake --build trilinos_build -j8 --target install
```

Troubleshooting: `-DTrilinos_ENABLE_Fortran=OFF` matters — without it,
CMake probes for a working Fortran compiler even though Sacado is pure
C++, and on a machine with a broken/absent Fortran toolchain
(`ld: library not found for -lSystem` from a stray `gfortran`) the whole
configure fails for a dependency Sacado never needed.

**2. Build SPARTA itself with AD enabled**, pointed at that Sacado install
(see [`BUILD_CMAKE.md`](../BUILD_CMAKE.md)'s AD section for the CMake
option reference and the one-Kokkos-version rule):

```bash
cmake -S cmake -B sparta_ad_build -C cmake/presets/mac.cmake \
  -DSPARTA_MACHINE=ad -DSPARTA_ENABLE_AD=ON \
  -DSACADO_ROOT=/path/to/sacado_install
cmake --build sparta_ad_build -j8
```

This produces `sparta_ad_build/src/libsparta_ad.a`.

**3. Configure `ad_src` against both**, and build:

```bash
cd ad_src
cmake -S . -B build \
  -DSPARTA_AD_BUILD_DIR=/path/to/sparta_ad_build \
  -DSACADO_ROOT=/path/to/sacado_install
cmake --build build
```

If both are found, this adds `shape_main_ad`/`shape_opt_ad` alongside the
stock targets. `SPARTA_AD_NDIR` (CMake cache var, default `4`) must match
whatever the `SPARTA_AD_BUILD_DIR` library was itself built with — it's an
ABI-affecting parameter (`sfloat`'s `Sacado::Fad::SFad<double,N>` size);
mismatched values on either side of the link boundary are undefined
behavior, not a performance question. `4` matches `BezierShape`'s
`ndesign()==4`, so a normal-default AD build gets a Bezier shape's whole
gradient from one SPARTA run.

If `SPARTA_AD_BUILD_DIR`/Sacado aren't found, the AD targets are skipped
with a `message(WARNING ...)` — the stock targets are unaffected either way.

## Running

All binaries look for `N.species` (and `N.vss`, unless `--nocoll`) in the
current directory — run from `ad_src/`.

```bash
# single-point evaluation, finite-difference gradient
./build/shape_main --alpha 1.3,1.0,2.7,0.8 --seeds 12345

# same, AD gradient (needs shape_main_ad)
./build/shape_main_ad --alpha 1.3,1.0,2.7,0.8 --seeds 12345

# IPOPT-driven optimization, 2D drag, with a minimum-area constraint
./build/shape_opt --alpha 1.3,1.0,2.7,0.8 --min-area 1.5 --experiment run1

# IPOPT-driven optimization, 3D heat flux (absorbs what used to be a
# separate standalone power_law_main.cpp driver)
./build/shape_opt_ad --shape powerlaw --objective heatflux \
  --alpha 1.4 --wall-accom 1.0 --experiment pl_run1
```

`shape_opt`/`shape_opt_ad` write `output/<fd|ad>_<date>_<experiment>/`:
`config.txt`, `ipopt.log`, `trajectory.csv`, `result.txt`,
`shape_initial.txt`/`shape_final.txt` (generic mesh dumps, both dims), and
`shapes.svg` (2D shapes only).

## Flipping the AD correction switch

Every driver takes `--score-correction` / `--no-score-correction` (default
off), with the effective state echoed in the startup banner and, for
`shape_opt(_ad)`, in `config.txt`/`result.txt` too. Passing it to a stock
(non-AD) build prints a warning and is otherwise ignored — it can't do
anything there.

```bash
./build/shape_main_ad --alpha 1.3,1.0,2.7,0.8 --score-correction
```

Raw escape hatch, for driving SPARTA directly outside these drivers:
`SPARTA_AD_SCORE_CORRECTION=1` env var, read by `src/compute_surf.cpp`'s
`ComputeSurf::init()`.

`tools/ad_verify/compare_modes.py` runs one config three ways (FD,
AD-uncorrected, AD-corrected) and tabulates value/gradient/AD-FD-ratio —
the fastest way to see what the switch actually changes for a given
objective:

```bash
python3 ../tools/ad_verify/compare_modes.py --alpha 1.3,1.0,2.7,0.8
```

**Important limitation**: the correction only reaches `fluxscale`-weighted
tallies (press/shear/heat-flux — i.e. `HeatFluxObjective`). It is a
documented no-op for `nfactor_inverse`-weighted tallies (force/torque —
i.e. `DragObjective`'s `fx`). See
[`docs/AD_GRADIENTS.md`](../docs/AD_GRADIENTS.md) for why, and
[`KNOWN_GAPS.md`](../KNOWN_GAPS.md) for the deferred fix.

## Extending

### Add a new `Shape`

Implement `shape.h`'s five methods (`dim`, `ndesign`, `bounds`, `to_mesh`,
`jacobian`; `validate` is optional). Look at `bezier_shape.h/.cpp` (2D) or
`power_law_shape.h/.cpp` (3D) for the pattern: a thin adapter wrapping a
pure geometry kernel (`bezier_geom.h`/`power_law_body.h`) that stays
SPARTA-free and array-in/array-out, so a different AD engine could seed
`jacobian()`'s output without this file changing. `jacobian()` takes
`alpha` — don't assume it's alpha-independent even if your first
implementation happens to be linear.

Wire it into `shape_main.cpp`/`shape_opt_main.cpp`'s `--shape` selection
(construct it, add an `if (cfg.shape_name == "...")` branch) and add its
CMakeLists.txt source files to every target that should offer it.

### Add a new `Objective`

If it's "sum one `compute surf ... <keyword>` value over the whole body"
(most QOIs are), just add a named subclass to `surf_sum_objective.h`
following `DragObjective`/`HeatFluxObjective`'s one-liner pattern —
`SurfSumObjective(keyword, tag)` handles the compute/fix/reduce chain and
the AD-vs-stock `extract()` split for you.

Otherwise, implement `objective.h` directly: `setup()` issues whatever
`compute`/`fix` commands your QOI needs (called once, right after
`read_surf`/`surf_modify`, before the settle/averaging run); `extract()`
reads the result back after the averaging window (and, in the AD build,
the gradient — go through `ad_extract()` in `ad_seed.h`, never poke
`.fastAccessDx()` directly, so the AD-engine seam stays in one place).

### Add a new `Constraint`

Implement `constraint.h`'s five methods. `min_size_constraint.cpp` is the
reference for both branches: the 2D shoelace-area gradient (linear in each
point, chained through `jacobian()`) and the 3D divergence-theorem volume
gradient (per-triangle scalar-triple-product cyclic identity, same
chaining) — copy whichever branch's structure matches your constraint's
own geometry, not necessarily both.

## Known gaps

See [`KNOWN_GAPS.md`](../KNOWN_GAPS.md) for the full list (the
`nfactor_inverse` AD correction, `DragObjective` on a 3D body being
untested, `shape_main` not offering `--shape powerlaw`, and others). None
of them are silent — each shows up as either a code comment at the exact
spot, a stderr warning at runtime, or an entry there.
