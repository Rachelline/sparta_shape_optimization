# Find Sacado (Trilinos forward-mode automatic differentiation).
#
# Hand-written Module-mode finder, matching the idiom of the other SPARTA
# dependency finders (FindFFTW3.cmake, FindMKL.cmake, and the inline IPOPT
# block in ad_src) -- deliberately NOT Trilinos's generated SacadoConfig.cmake,
# so Sacado is located the same way as every other TPL in this repo.
#
# Only the Sacado headers + libsacado are referenced. Per the one-Kokkos rule
# (see BUILD_CMAKE.md's "Automatic differentiation (AD)" section and
# cmake/common/process/sparta_build_options.cmake), the AD+Kokkos build must
# get Kokkos from SPARTA's single Kokkos, NOT from Sacado's install -- so this
# finder intentionally does NOT expose Sacado's bundled include/kokkos
# directory.
#
# Point at a minimal Trilinos install (Sacado [+ Kokkos]) via -DSACADO_ROOT=...
# or the SACADO_ROOT environment variable.
#
# Sets:
#   SACADO_INCLUDE_DIRS, SACADO_LIBRARIES, Sacado_FOUND
#   imported target Sacado::Sacado

find_path(SACADO_INCLUDE_DIR
  NAMES Sacado.hpp
  HINTS ${SACADO_ROOT} $ENV{SACADO_ROOT}
  PATH_SUFFIXES include)

find_library(SACADO_LIBRARY
  NAMES sacado
  HINTS ${SACADO_ROOT} $ENV{SACADO_ROOT}
  PATH_SUFFIXES lib lib64)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Sacado DEFAULT_MSG
  SACADO_LIBRARY SACADO_INCLUDE_DIR)

if(Sacado_FOUND)
  set(SACADO_INCLUDE_DIRS ${SACADO_INCLUDE_DIR})
  set(SACADO_LIBRARIES ${SACADO_LIBRARY})

  if(NOT TARGET Sacado::Sacado)
    add_library(Sacado::Sacado UNKNOWN IMPORTED)
    set_target_properties(Sacado::Sacado PROPERTIES
      IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
      IMPORTED_LOCATION "${SACADO_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${SACADO_INCLUDE_DIRS}")
  endif()
endif()

mark_as_advanced(SACADO_INCLUDE_DIR SACADO_LIBRARY)
