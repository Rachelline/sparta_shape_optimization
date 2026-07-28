# ----------------------------------------------------------------------
# Standalone build for ad_src/'s shape-optimization driver (Milestone 1:
# Parametrization x Objective, FD gradients only -- no AD, no IPOPT yet).
#
# This is a self-contained CMake project, deliberately NOT wired into
# sparta/cmake or sparta/src/CMakeLists.txt (no add_subdirectory anywhere
# references ad_src). Configuring/building the main SPARTA CMake tree
# (sparta/build, with its SPARTA_ENABLE_AD/SPARTA_AD_NDIR options) is
# completely unaffected by this file -- structure ported from the
# reference implementation at
#   /Users/olenamelnikov/Desktop/sparta_AD/sparta/ad_src/
# which shells out the exact same way, for the exact same reason (that
# fork has no CMake AD support at all, only the classic Makefile path;
# this fork's own SPARTA_ENABLE_AD CMake option is a DIFFERENT mechanism
# that this file does not use -- see the deferred-AD note below).
#
# SPARTA itself is built as a side effect, via a CMake custom command
# that shells out to SPARTA's own classic src/Makefile in "mode=lib"
# (same as `cd src && make mode=lib serial` upstream): produces
# src/libsparta_serial.a and src/STUBS/libmpi_stubs.a. Those live in
# sparta/src/, alongside (but not overwriting) whatever per-machine .o
# dir (Obj_serial) the classic Makefile already uses.
#
# Deferred (see docs/PLAN.md "What this milestone deliberately does NOT
# include" and the plan this was built from): an AD-linked target. This
# fork has no src/MAKE/Makefile.serial_ad for a classic-Makefile shell-
# out to target (unlike the reference fork); the better long-term answer
# is almost certainly linking this fork's REAL SPARTA_ENABLE_AD CMake
# build directly (reusing FindSacado.cmake/SPARTA_AD_NDIR) instead of
# adding a second, parallel classic-Makefile AD path here -- left as
# follow-up work, not decided in this file.
#
# Usage, from ad_src/:
#   cmake -S . -B build
#   cmake --build build
#   ./build/shape_main --alpha 1.3,1.0,2.7,0.8   # run from ad_src/ itself,
#                                                  # where N.species/N.vss
#                                                  # live (matches the
#                                                  # reference's own
#                                                  # "run from a dir with
#                                                  # air.species/air.vss"
#                                                  # convention)
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

# ---- shape_main: the framework, FD gradients only ----------------------
# Links against SPARTA through its public C interface (src/library.h),
# exactly like the reference's drag_run.

add_executable(shape_main
  shape_main.cpp
  shape_case.cpp
  bezier_geom.cpp
  bezier_parametrization.cpp
  drag_objective.cpp
  heat_flux_objective.cpp)
target_compile_definitions(shape_main PRIVATE SPARTA_UNORDERED_MAP)
target_include_directories(shape_main PRIVATE ${SPARTA_SRC_DIR} ${SPARTA_SRC_DIR}/STUBS)
target_link_libraries(shape_main PRIVATE sparta_serial mpi_stubs)
