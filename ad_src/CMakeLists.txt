# ----------------------------------------------------------------------
# Standalone build for ad_src/'s shape-optimization driver.
# Milestone 1: Parametrization x Objective, FD gradients (shape_main).
# Milestone 2: ShapeTNLP + IPOPT optimizer loop (shape_opt).
# Milestone 3: AD-linked targets (shape_main_ad/shape_opt_ad), a general
#   Constraint mechanism (--min-area), and a standalone 3D power-law-body
#   heatflux example. See the AD-targets section below for how
#   SPARTA_AD_BUILD_DIR/Sacado discovery works.
# Milestone 4: cleanup, 2D/3D unification (Parametrization replaced by
#   shape.h's Shape interface, used by both BezierShape and
#   PowerLawShape; MinAreaConstraint generalized to MinSizeConstraint),
#   and driver consolidation: the standalone power_law_main.cpp /
#   power_law_opt(_ad) targets are gone -- shape_opt(_ad)'s
#   --shape powerlaw absorbs that capability via the same Shape/
#   ShapeTNLP/Constraint machinery the 2D case already used.
#
# This is a self-contained CMake project, deliberately NOT wired into
# sparta/cmake or sparta/src/CMakeLists.txt (no add_subdirectory anywhere
# references ad_src). Configuring/building the main SPARTA CMake tree
# (sparta/build, with its SPARTA_ENABLE_AD/SPARTA_AD_NDIR options) is
# completely unaffected by this file.
#
# The STOCK (FD-gradient) side of this project builds SPARTA as a side
# effect, via a CMake custom command that shells out to SPARTA's own
# classic src/Makefile in "mode=lib" (same as `cd src && make mode=lib
# serial` upstream): produces src/libsparta_serial.a and
# src/STUBS/libmpi_stubs.a. Those live in sparta/src/, alongside (but not
# overwriting) whatever per-machine .o dir (Obj_serial) the classic
# Makefile already uses.
#
# The AD side does NOT use that classic-Makefile path (this fork has no
# src/MAKE/Makefile.serial_ad) -- it instead links against an already-
# built top-level SPARTA CMake build configured with
# -DSPARTA_ENABLE_AD=ON, pointed at via the SPARTA_AD_BUILD_DIR cache
# variable (see the AD-targets section below for the exact recipe).
#
# Usage, from ad_src/:
#   cmake -S . -B build
#   cmake --build build
#   ./build/shape_main --alpha 1.3,1.0,2.7,0.8   # run from ad_src/ itself,
#                                                  # where N.species/N.vss live
# See ad_src/README.md for the full build recipe, including the AD side
# (Trilinos/Sacado + a top-level -DSPARTA_ENABLE_AD=ON SPARTA build).
# ----------------------------------------------------------------------

cmake_minimum_required(VERSION 3.16)
project(sparta_shape_opt NONE)
enable_language(CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(SPARTA_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/..)
set(SPARTA_SRC_DIR ${SPARTA_ROOT}/src)
find_program(MAKE_EXECUTABLE NAMES make gmake REQUIRED)

# ---- IPOPT discovery (for shape_opt) -----------------------------------
# IPOPT is optional: if not found, shape_main still builds and only
# shape_opt is skipped. Install with `brew install ipopt` (headers ->
# /opt/homebrew/include/coin-or, lib -> /opt/homebrew/lib/libipopt.dylib).
# No pkg-config needed. (Confirmed present this session: 3.14.19, which
# satisfies get_curr_iterate()'s >=3.14 requirement -- see shape_tnlp.cpp.)
find_path(IPOPT_INCLUDE_DIR NAMES IpTNLP.hpp
  PATHS /opt/homebrew/include/coin-or /usr/local/include/coin-or
        /opt/homebrew/include /usr/local/include)
find_library(IPOPT_LIB NAMES ipopt
  PATHS /opt/homebrew/lib /usr/local/lib)
if(IPOPT_INCLUDE_DIR AND IPOPT_LIB)
  set(HAVE_IPOPT ON)
  message(STATUS "IPOPT found: ${IPOPT_LIB} (include ${IPOPT_INCLUDE_DIR})")
else()
  set(HAVE_IPOPT OFF)
  message(WARNING
    "IPOPT not found -- skipping shape_opt.\n"
    "Install it with `brew install ipopt`, then re-run cmake.")
endif()

# ---- STUBS (MPI stand-in for serial builds) ---------------------------

set(SPARTA_MPI_STUBS_A ${SPARTA_SRC_DIR}/STUBS/libmpi_stubs.a)
add_custom_command(
  OUTPUT ${SPARTA_MPI_STUBS_A}
  COMMAND ${MAKE_EXECUTABLE} -C ${SPARTA_SRC_DIR}/STUBS
  COMMENT "Building SPARTA MPI stubs (src/STUBS)")
add_custom_target(sparta_build_mpi_stubs DEPENDS ${SPARTA_MPI_STUBS_A})
add_library(mpi_stubs STATIC IMPORTED)
set_target_properties(mpi_stubs PROPERTIES IMPORTED_LOCATION ${SPARTA_MPI_STUBS_A})
add_dependencies(mpi_stubs sparta_build_mpi_stubs)

# ---- stock SPARTA library: classic `make mode=lib serial` -------------

# NOTE: src/MAKE/Makefile.serial in THIS fork lacks -DSPARTA_UNORDERED_MAP
# (present in the reference fork's copy) -- without it, src/hash3.h falls
# back to <tr1/unordered_map>, which modern Clang/libc++ no longer ships,
# and the classic-Makefile build fails outright. Passed as a make
# command-line override here (Makefile.serial's SPARTA_INC uses plain
# `=`, no `override` directive, so this safely wins) rather than editing
# the checked-in machine file for what's really a pre-existing toolchain
# gap, unrelated to this feature.
set(SPARTA_SERIAL_A ${SPARTA_SRC_DIR}/libsparta_serial.a)
add_custom_command(
  OUTPUT ${SPARTA_SERIAL_A}
  COMMAND ${MAKE_EXECUTABLE} -C ${SPARTA_SRC_DIR} mode=lib serial
          SPARTA_INC=-DSPARTA_GZIP\ -DSPARTA_UNORDERED_MAP
  DEPENDS ${SPARTA_MPI_STUBS_A}
  COMMENT "Building SPARTA (serial, static lib)")
add_custom_target(sparta_build_serial DEPENDS ${SPARTA_SERIAL_A})
add_library(sparta_serial STATIC IMPORTED)
set_target_properties(sparta_serial PROPERTIES IMPORTED_LOCATION ${SPARTA_SERIAL_A})
add_dependencies(sparta_serial sparta_build_serial)

# ---- shape_main: the framework, single-point evaluator ----------------
# Links against SPARTA through its public C interface (src/library.h),
# exactly like the reference's drag_run.

add_executable(shape_main
  shape_main.cpp
  shape_case.cpp
  bezier_geom.cpp
  bezier_shape.cpp
  polygon2d.cpp
  surf_sum_objective.cpp
  sparta_util.cpp
  cli.cpp)
target_compile_definitions(shape_main PRIVATE SPARTA_UNORDERED_MAP)
target_include_directories(shape_main PRIVATE ${SPARTA_SRC_DIR} ${SPARTA_SRC_DIR}/STUBS)
target_link_libraries(shape_main PRIVATE sparta_serial mpi_stubs)

# ---- shape_opt: IPOPT-driven optimizer loop, both dims -----------------
# Same SPARTA dependency as shape_main, plus IPOPT's TNLP interface.
# --shape bezier|powerlaw selects the 2D or 3D Shape at runtime -- see
# shape_opt_main.cpp for the absorbed power-law-body demo.

if(HAVE_IPOPT)
  add_executable(shape_opt
    shape_opt_main.cpp
    shape_tnlp.cpp
    svg_shape.cpp
    shape_case.cpp
    bezier_geom.cpp
    bezier_shape.cpp
    polygon2d.cpp
    power_law_body.cpp
    power_law_shape.cpp
    surf_sum_objective.cpp
    min_size_constraint.cpp
    sparta_util.cpp
    cli.cpp
    run_output.cpp)
  target_compile_definitions(shape_opt PRIVATE SPARTA_UNORDERED_MAP)
  target_include_directories(shape_opt PRIVATE
    ${SPARTA_SRC_DIR} ${SPARTA_SRC_DIR}/STUBS ${IPOPT_INCLUDE_DIR})
  target_link_libraries(shape_opt PRIVATE sparta_serial mpi_stubs ${IPOPT_LIB})
endif()

# ---- AD-linked targets: shape_main_ad / shape_opt_ad (Milestone 3) -----
# AD gradients carry a known low bias on DSMC surface tallies (docs/
# ad_phase_c_investigation/FINDINGS.md). A runtime toggle
# (SPARTA_AD_SCORE_CORRECTION env var / --score-correction CLI flag)
# corrects it for fluxscale-weighted tallies (press/shear/heat-flux);
# off by default, so gradients are used exactly as the solver returns
# them unless explicitly requested otherwise. force/torque tallies
# (nfactor_inverse-scaled) are not yet covered -- see docs/AD_GRADIENTS.md.
# Optional: skipped entirely (stock targets above unaffected) unless both
# SPARTA_AD_BUILD_DIR and Sacado are found.
#
# SPARTA_AD_BUILD_DIR: path to an ALREADY-configured-and-built top-level
# SPARTA CMake build (cmake/CMakeLists.txt) with -DSPARTA_ENABLE_AD=ON,
# e.g.:
#   cmake -S cmake -B <dir> -C cmake/presets/mac.cmake \
#     -DSPARTA_MACHINE=ad -DSPARTA_ENABLE_AD=ON -DSACADO_ROOT=<sacado install>
#   cmake --build <dir>
# Not auto-built by this file: the top-level AD path has real complexity
# (Kokkos ODR guard, C++17+ bump, Sacado discovery) not worth duplicating
# in a nested custom command. See ad_src/README.md for the full recipe.
#
# SPARTA_AD_NDIR must match whatever SPARTA_AD_BUILD_DIR was itself
# configured with (default 4 on both sides -- matches BezierShape's
# ndesign()==4, letting a normal-default AD build get a shape's whole
# gradient from one SPARTA run). This is an ABI-affecting
# parameter (sfloat's Sacado::Fad::SFad<double,N> size); mismatched values
# on either side of the link boundary is undefined behavior, not just a
# performance question.
set(SPARTA_AD_BUILD_DIR "" CACHE PATH
  "Path to a pre-built -DSPARTA_ENABLE_AD=ON top-level SPARTA CMake build. Empty = skip shape_main_ad/shape_opt_ad.")
set(SPARTA_AD_NDIR "4" CACHE STRING
  "Must match the SPARTA_AD_NDIR the SPARTA_AD_BUILD_DIR library was built with.")

if(SPARTA_AD_BUILD_DIR)
  list(APPEND CMAKE_MODULE_PATH ${SPARTA_ROOT}/cmake/Modules)
  find_package(Sacado)

  find_library(SPARTA_AD_LIB NAMES sparta_ad sparta
    PATHS ${SPARTA_AD_BUILD_DIR}/src NO_DEFAULT_PATH)
  if(NOT SPARTA_AD_LIB)
    file(GLOB _sparta_ad_lib_candidates ${SPARTA_AD_BUILD_DIR}/src/libsparta_*.a)
    if(_sparta_ad_lib_candidates)
      list(GET _sparta_ad_lib_candidates 0 SPARTA_AD_LIB)
    endif()
  endif()

  if(Sacado_FOUND AND SPARTA_AD_LIB)
    message(STATUS "AD SPARTA found: ${SPARTA_AD_LIB} (Sacado: ${SACADO_INCLUDE_DIRS})")
    add_library(sparta_ad_lib STATIC IMPORTED)
    set_target_properties(sparta_ad_lib PROPERTIES IMPORTED_LOCATION ${SPARTA_AD_LIB})

    set(AD_COMPILE_DEFS SPARTA_AD SPARTA_AD_NDIR=${SPARTA_AD_NDIR} SPARTA_UNORDERED_MAP)

    add_executable(shape_main_ad
      shape_main.cpp shape_case.cpp bezier_geom.cpp bezier_shape.cpp polygon2d.cpp
      surf_sum_objective.cpp sparta_util.cpp cli.cpp)
    target_compile_definitions(shape_main_ad PRIVATE ${AD_COMPILE_DEFS})
    target_include_directories(shape_main_ad PRIVATE ${SPARTA_SRC_DIR} ${SPARTA_SRC_DIR}/STUBS)
    target_link_libraries(shape_main_ad PRIVATE sparta_ad_lib mpi_stubs Sacado::Sacado)
    # Sacado's headers need C++17+ (matches the top-level build's own bump,
    # cmake/common/process/sparta_build_options.cmake: "AD needs at least
    # C++17"); the project-wide CMAKE_CXX_STANDARD above (14) is fine for
    # the stock targets but not this one.
    set_target_properties(shape_main_ad PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON)

    if(HAVE_IPOPT)
      add_executable(shape_opt_ad
        shape_opt_main.cpp shape_tnlp.cpp svg_shape.cpp shape_case.cpp
        bezier_geom.cpp bezier_shape.cpp polygon2d.cpp
        power_law_body.cpp power_law_shape.cpp
        surf_sum_objective.cpp min_size_constraint.cpp
        sparta_util.cpp cli.cpp run_output.cpp)
      target_compile_definitions(shape_opt_ad PRIVATE ${AD_COMPILE_DEFS})
      target_include_directories(shape_opt_ad PRIVATE
        ${SPARTA_SRC_DIR} ${SPARTA_SRC_DIR}/STUBS ${IPOPT_INCLUDE_DIR})
      target_link_libraries(shape_opt_ad PRIVATE sparta_ad_lib mpi_stubs Sacado::Sacado ${IPOPT_LIB})
      set_target_properties(shape_opt_ad PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON)
    endif()
  else()
    message(WARNING
      "SPARTA_AD_BUILD_DIR set but Sacado and/or the AD library weren't "
      "found there -- skipping shape_main_ad/shape_opt_ad.")
  endif()
endif()
