# ##############################################################################
# This file process boolean options set by
# cmake/common/set/sparta_build_options.cmake and sets up the correct SPARTA
# targets based on those options. This file is also responsible for handeling
# dependencies among boolean options.
# ##############################################################################

# ################### BEGIN PROCESS EXTRA OPTIONS ####################
if(SPARTA_ENABLE_ALL_PKGS AND SPARTA_DISABLE_ALL_PKGS)
  message(FATAL_ERROR "SPARTA_ENABLE_ALL_PKGS=SPARTA_DISABLE_ALL_PKGS=ON")
endif()

if(SPARTA_ENABLE_ALL_PKGS)
  foreach(pkg IN LISTS SPARTA_PKG_LIST)
    get_property(
      pkg_help
      CACHE ${pkg}
      PROPERTY HELPSTRING)
    set(${pkg}
        ON
        CACHE STRING "${pkg_help}" FORCE)
  endforeach()
endif()

if(SPARTA_DISABLE_ALL_PKGS)
  foreach(pkg IN LISTS SPARTA_PKG_LIST)
    get_property(
      pkg_help
      CACHE ${pkg}
      PROPERTY HELPSTRING)
    set(${pkg}
        OFF
        CACHE STRING "${pkg_help}" FORCE)
  endforeach()
endif()
# ################### END   PROCESS EXTRA OPTIONS ####################

# ################### BEGIN PROCESS MPI TPL/PKG ####################
if(BUILD_MPI AND PKG_MPI_STUBS)
  message(WARNING "Both BUILD_MPI: ${BUILD_MPI} and PKG_MPI_STUBS: "
                  "${PKG_MPI_STUBS} are selected. Defaulting to PKG_MPI_STUBS.")
endif()

if(BUILD_MPI AND NOT PKG_MPI_STUBS)
  find_package(MPI REQUIRED)
  # TODO: if NOT MPI_FOUND, handle finding mpi installs
  set(TARGET_SPARTA_BUILD_MPI MPI::MPI_CXX)
else()
  set(PKG_MPI_STUBS ON)
  set(BUILD_MPI OFF)
  set(TARGET_SPARTA_BUILD_MPI pkg_mpi_stubs)
endif()
list(APPEND TARGET_SPARTA_BUILD_TPLS ${TARGET_SPARTA_BUILD_MPI})
# ################### END PROCESS MPI TPL/PKG ####################

# ################### BEGIN PROCESS FFT TPL/PKG ####################

if((NOT PKG_FFT) AND (NOT (FFT_KOKKOS STREQUAL "OFF")))
  message(FATAL_ERROR  "Setting FFT_KOKKOS library requires PKG_FFT: ON.")
endif()

if((NOT PKG_FFT) AND (NOT (FFT STREQUAL "OFF")))
  message(FATAL_ERROR  "Setting FFT library requires PKG_FFT: ON.")
endif()

if(PKG_FFT)

  if(FFT STREQUAL "OFF")
    set(FFT "KISS" CACHE STRING "Select a FFT TPL for CPU from FFTW3, MKL, or KISS. Default: KISS." FORCE)
  endif()

  set(TARGET_SPARTA_PKG_FFT pkg_fft)
  list(APPEND TARGET_SPARTA_PKGS ${TARGET_SPARTA_PKG_FFT})

  string(TOUPPER ${FFT_KOKKOS} FFT_KOKKOS)

  set(SPARTA_DEFAULT_CXX_COMPILE_FLAGS -DFFT_${FFT}
                                       ${SPARTA_DEFAULT_CXX_COMPILE_FLAGS})

  if(NOT FFT STREQUAL "KISS")
    find_package(${FFT} REQUIRED)

    set(TARGET_SPARTA_BUILD_FFT ${FFT}::${FFT})
    list(APPEND TARGET_SPARTA_BUILD_TPLS ${TARGET_SPARTA_BUILD_FFT})

    set(TARGET_SPARTA_PKG_FFT pkg_fft)
    list(APPEND TARGET_SPARTA_PKGS ${TARGET_SPARTA_PKG_FFT})
  endif()

  if(SPARTA_ENABLE_TESTING)
    set(SPARTA_ENABLED_TEST_SUITES ${SPARTA_ENABLED_TEST_SUITES} "fft")
  endif()

  if(PKG_KOKKOS)
    if(FFT_KOKKOS STREQUAL "OFF")
      set(FFT_KOKKOS "KISS" CACHE STRING "Select a FFT TPL for Kokkos from CUFFT, HIPFFT, FFTW3, MKL, MKL_GPU, or KISS. Default: KISS." FORCE)
    endif()

    string(TOUPPER ${FFT_KOKKOS} FFT_KOKKOS)
    if(USE_EXTERNAL_KOKKOS)
      find_package(Kokkos REQUIRED CONFIG)
    endif()

    if(FFT_KOKKOS STREQUAL "CUFFT" AND NOT Kokkos_ENABLE_CUDA)
      message(FATAL_ERROR  "FFT_KOKKOS: ${FFT_KOKKOS} requires Kokkos_ENABLE_CUDA: ON.")
    endif()

    if(FFT_KOKKOS STREQUAL "HIPFFT" AND NOT Kokkos_ENABLE_HIP)
      message(FATAL_ERROR  "FFT_KOKKOS: ${FFT_KOKKOS} requires Kokkos_ENABLE_HIP: ON.")
    endif()

    if(FFT_KOKKOS STREQUAL "MKL_GPU" AND NOT Kokkos_ENABLE_SYCL)
      message(FATAL_ERROR  "FFT_KOKKOS: ${FFT_KOKKOS} requires Kokkos_ENABLE_SYCL: ON.")
    endif()

    if(FFT_KOKKOS STREQUAL "FFTW3" OR FFT_KOKKOS STREQUAL "MKL")
      if((NOT Kokkos_ENABLE_SERIAL) AND (NOT Kokkos_ENABLE_OPENMP))
        message(FATAL_ERROR  "FFT_KOKKOS: ${FFT_KOKKOS} requires either Kokkos_ENABLE_OPENMP or Kokkos_ENABLE_SERIAL")
      endif()

      if(Kokkos_ENABLE_CUDA OR Kokkos_ENABLE_HIP OR Kokkos_ENABLE_ROCM OR Kokkos_ENABLE_SYCL)
        message(FATAL_ERROR  "FFT_KOKKOS: ${FFT_KOKKOS} cannot run with a kokkos GPU backend.")
      endif()
    endif()

    if(Kokkos_ENABLE_CUDA)
      if(NOT ((FFT_KOKKOS STREQUAL "KISS") OR (FFT_KOKKOS STREQUAL "CUFFT")))
        message(FATAL_ERROR "The CUDA backend of Kokkos requires either KISS FFT or CUFFT.")
      elseif(FFT_KOKKOS STREQUAL "KISS")
        message(WARNING "Using KISS FFT with the CUDA backend of Kokkos may be sub-optimal.")
      elseif(FFT_KOKKOS STREQUAL "CUFFT")
        find_library(CUFFT_LIBRARY cufft)
        if (CUFFT_LIBRARY STREQUAL "CUFFT_LIBRARY-NOTFOUND")
          message(FATAL_ERROR "Required cuFFT library not found. Check your environment or set CUFFT_LIBRARY to its location")
        endif()
        link_libraries(${CUFFT_LIBRARY})
      endif()
    elseif(Kokkos_ENABLE_HIP)
      if(NOT ((FFT_KOKKOS STREQUAL "KISS") OR (FFT_KOKKOS STREQUAL "HIPFFT")))
        message(FATAL_ERROR "The HIP backend of Kokkos requires either KISS FFT or HIPFFT.")
      elseif(FFT_KOKKOS STREQUAL "KISS")
        message(WARNING "Using KISS FFT with the HIP backend of Kokkos may be sub-optimal.")
      elseif(FFT_KOKKOS STREQUAL "HIPFFT")
        include(DetectHIPInstallation)
        find_package(hipfft REQUIRED)
        link_libraries(hip::hipfft)
      endif()
    elseif(Kokkos_ENABLE_SYCL)
      if(NOT ((FFT_KOKKOS STREQUAL "KISS") OR (FFT_KOKKOS STREQUAL "MKL_GPU")))
        message(FATAL_ERROR "The SYCL backend of Kokkos requires either KISS FFT or MKL_GPU.")
      elseif(FFT_KOKKOS STREQUAL "KISS")
        message(WARNING "Using KISS FFT with the SYCL backend of Kokkos may be sub-optimal.")
      elseif(FFT_KOKKOS STREQUAL "MKL_GPU")
        find_package(MKL REQUIRED)
        link_libraries(mkl_sycl_dft mkl_intel_ilp64 mkl_tbb_thread mkl_core tbb)
      endif()
    endif()

    if(FFT_KOKKOS STREQUAL "FFTW3" OR FFT_KOKKOS STREQUAL "MKL")
      find_package(${FFT_KOKKOS} REQUIRED)
    endif()

    set(SPARTA_DEFAULT_CXX_COMPILE_FLAGS -DFFT_KOKKOS_${FFT_KOKKOS}
                                         ${SPARTA_DEFAULT_CXX_COMPILE_FLAGS})

  endif()

endif()

# ################### BEGIN PROCESS PKGS ####################

if(PKG_KOKKOS)

  # As of version 5.0.0 Kokkos requires C++20
  if(CMAKE_CXX_STANDARD LESS 20)
    message(FATAL_ERROR "The KOKKOS package requires the C++ standard to
  be set to at least C++20")
  endif()

########################################################################
# consistency checks and Kokkos options/settings required by SPARTA
if(Kokkos_ENABLE_HIP)
  option(Kokkos_ENABLE_HIP_MULTIPLE_KERNEL_INSTANTIATIONS "Enable multiple kernel instantiations with HIP" ON)
  mark_as_advanced(Kokkos_ENABLE_HIP_MULTIPLE_KERNEL_INSTANTIATIONS)
  option(Kokkos_ENABLE_ROCTHRUST "Use RoCThrust library" ON)
  mark_as_advanced(Kokkos_ENABLE_ROCTHRUST)
endif()

if(Kokkos_ENABLE_SERIAL)
  if(NOT (Kokkos_ENABLE_OPENMP OR Kokkos_ENABLE_THREADS OR
    Kokkos_ENABLE_CUDA OR Kokkos_ENABLE_HIP OR Kokkos_ENABLE_SYCL
    OR Kokkos_ENABLE_OPENMPTARGET))
  option(Kokkos_ENABLE_ATOMICS_BYPASS "Disable atomics for Kokkos Serial Backend" ON)
  mark_as_advanced(Kokkos_ENABLE_ATOMICS_BYPASS)
  endif()
endif()
########################################################################

  set(TARGET_SPARTA_PKG_KOKKOS pkg_kokkos)
  list(APPEND TARGET_SPARTA_PKGS ${TARGET_SPARTA_PKG_KOKKOS})
  set(SPARTA_DEFAULT_CXX_COMPILE_FLAGS -DSPARTA_KOKKOS
                                       ${SPARTA_DEFAULT_CXX_COMPILE_FLAGS})
  # SPARTA_KOKKOS_EXACT makes the KOKKOS package reproduce non-KOKKOS results
  # exactly, which allows the KOKKOS build to be regression tested against the
  # existing (non-KOKKOS) gold-standard log files.
  if(SPARTA_KOKKOS_EXACT)
    set(SPARTA_DEFAULT_CXX_COMPILE_FLAGS -DSPARTA_KOKKOS_EXACT
                                         ${SPARTA_DEFAULT_CXX_COMPILE_FLAGS})
  endif()
  # PKG_KOKKOS depends on BUILD_KOKKOS
  set(BUILD_KOKKOS ON)
endif()
# ################### END PROCESS PKGS ####################

# ################### BEGIN PROCESS TPLS ####################
if(BUILD_KOKKOS)
  # The 2 lines below find an install version of Kokkos. We are using and
  # in-tree copy of kokkos in sparta for now. find_package(KOKKOS REQUIRED)
  # set(TARGET_SPARTA_BUILD_KOKKOS Kokkos::kokkos)
  if(USE_EXTERNAL_KOKKOS)
    find_package(Kokkos REQUIRED CONFIG)
    set(TARGET_SPARTA_BUILD_KOKKOS Kokkos::kokkos)
  else()
    set(TARGET_SPARTA_BUILD_KOKKOS kokkos)
  endif()
  list(APPEND TARGET_SPARTA_BUILD_TPLS ${TARGET_SPARTA_BUILD_KOKKOS})
  # BUILD_KOKKOS does not depend on PKG_KOKKOS, do not attempt to resolve
  # dependency
endif()

if(BUILD_JPEG)
  find_package(JPEG REQUIRED)
  set(TARGET_SPARTA_BUILD_JPEG JPEG::JPEG)
  list(APPEND TARGET_SPARTA_BUILD_TPLS ${TARGET_SPARTA_BUILD_JPEG})
  set(SPARTA_DEFAULT_CXX_COMPILE_FLAGS -DSPARTA_JPEG
                                       ${SPARTA_DEFAULT_CXX_COMPILE_FLAGS})
endif()

if(BUILD_PNG)
  find_package(PNG REQUIRED)
  set(TARGET_SPARTA_BUILD_PNG PNG::PNG)
  list(APPEND TARGET_SPARTA_BUILD_TPLS ${TARGET_SPARTA_BUILD_PNG})
  set(SPARTA_DEFAULT_CXX_COMPILE_FLAGS -DSPARTA_PNG
                                       ${SPARTA_DEFAULT_CXX_COMPILE_FLAGS})
endif()

# AD (Sacado forward-mode automatic differentiation): sfloat = SFad.
# Header-only SFad, so Sacado::Sacado just adds an include path + libsacado.
if(SPARTA_ENABLE_AD)
  find_package(Sacado REQUIRED)
  list(APPEND TARGET_SPARTA_BUILD_TPLS Sacado::Sacado)
  set(SPARTA_DEFAULT_CXX_COMPILE_FLAGS
      -DSPARTA_AD -DSPARTA_AD_NDIR=${SPARTA_AD_NDIR}
      ${SPARTA_DEFAULT_CXX_COMPILE_FLAGS})

  # AD needs at least C++17. Kokkos already forces C++20 (cmake/CMakeLists.txt),
  # so this only lifts the AD-without-Kokkos build off the C++11 default.
  if(CMAKE_CXX_STANDARD LESS 17)
    set(CMAKE_CXX_STANDARD 17)
  endif()

  # -------- One-Kokkos guard (AD + Kokkos) --------
  # Sacado's Kokkos and SPARTA's Kokkos MUST be the exact same version, or C++
  # ODR is violated (two definitions of Kokkos::View etc.) -> silent wrong
  # results, especially on GPU. SFad is header-only, so it binds to whichever
  # Kokkos is on the include path; the AD+Kokkos build must therefore use ONE
  # Kokkos (USE_EXTERNAL_KOKKOS + Kokkos_ROOT = the Trilinos install Sacado was
  # built against), never the unmanaged in-tree lib/kokkos.
  if(PKG_KOKKOS)
    if(NOT USE_EXTERNAL_KOKKOS)
      message(FATAL_ERROR
        "SPARTA_ENABLE_AD + PKG_KOKKOS requires USE_EXTERNAL_KOKKOS=ON pointing "
        "(via Kokkos_ROOT) at the SAME Kokkos that Sacado was built against. "
        "The in-tree lib/kokkos is a separate, version-unmanaged Kokkos build; "
        "linking it alongside Sacado's Kokkos is an ODR violation.")
    endif()

    # SPARTA's expected Kokkos version, parsed from lib/kokkos.
    file(STRINGS "${SPARTA_TPL_DIR}/kokkos/CMakeLists.txt" _spa_kv
         REGEX "set\\(Kokkos_VERSION_(MAJOR|MINOR|PATCH)")
    string(REGEX MATCH "MAJOR ([0-9]+)" _m "${_spa_kv}")
    set(_spa_maj "${CMAKE_MATCH_1}")
    string(REGEX MATCH "MINOR ([0-9]+)" _m "${_spa_kv}")
    set(_spa_min "${CMAKE_MATCH_1}")
    string(REGEX MATCH "PATCH ([0-9]+)" _m "${_spa_kv}")
    set(_spa_pat "${CMAKE_MATCH_1}")
    set(SPARTA_KOKKOS_EXPECTED "${_spa_maj}.${_spa_min}.${_spa_pat}")

    # Kokkos_VERSION is set by the find_package(Kokkos CONFIG) that the
    # USE_EXTERNAL_KOKKOS path already ran above.
    if(DEFINED Kokkos_VERSION AND SPARTA_KOKKOS_EXPECTED MATCHES "^[0-9]+\\.[0-9]+\\.[0-9]+$")
      if(NOT Kokkos_VERSION VERSION_EQUAL SPARTA_KOKKOS_EXPECTED)
        message(FATAL_ERROR
          "AD+Kokkos version mismatch: external Kokkos (Sacado's) is "
          "${Kokkos_VERSION} but SPARTA lib/kokkos expects "
          "${SPARTA_KOKKOS_EXPECTED}. Point Sacado/Kokkos at the Trilinos "
          "release whose bundled Kokkos == ${SPARTA_KOKKOS_EXPECTED} "
          "(e.g. trilinos-release-17-0-0 ships Kokkos 5.0.2). "
          "Rule: match Sacado to Kokkos -- bump Sacado, never downgrade Kokkos.")
      endif()
    else()
      message(WARNING
        "AD+Kokkos: could not determine both Kokkos versions to verify they "
        "match (SPARTA expects '${SPARTA_KOKKOS_EXPECTED}', external Kokkos_VERSION "
        "='${Kokkos_VERSION}'). Ensure they are identical to avoid an ODR "
        "violation.")
    endif()
  endif()
endif()

if(PKG_PYTHON)

  set(TARGET_SPARTA_PKG_PYTHON pkg_python)
  list(APPEND TARGET_SPARTA_PKGS ${TARGET_SPARTA_PKG_PYTHON})
  set(SPARTA_DEFAULT_CXX_COMPILE_FLAGS -DSPARTA_PYTHON
                                       ${SPARTA_DEFAULT_CXX_COMPILE_FLAGS})

  if(NOT Python_INTERPRETER)
    find_package(Python COMPONENTS Interpreter)
  endif()
  find_package(Python REQUIRED COMPONENTS Interpreter Development)
  link_libraries(Python::Python)
  set(SPARTA_DEFAULT_CXX_COMPILE_FLAGS -DSPARTA_PYTHON
                                       ${SPARTA_DEFAULT_CXX_COMPILE_FLAGS})
endif()

# ################### END PROCESS TPLS ####################

# ################### BEGIN COMBINE CXX FLAGS ####################
set(SPARTA_DEFAULT_CXX_COMPILE_FLAGS ${SPARTA_CXX_COMPILE_FLAGS}
                                     ${SPARTA_DEFAULT_CXX_COMPILE_FLAGS})
# ################### BEGIN COMBINE CXX FLAGS ####################

if(BUILD_MPI)
  set_property(
    TARGET ${TARGET_SPARTA_BUILD_MPI}
    PROPERTY INTERFACE_COMPILE_OPTIONS ${SPARTA_DEFAULT_CXX_COMPILE_FLAGS})

  # -DSPARTA_AD above reaches every package target (FFT, PYTHON, KOKKOS, ...)
  # because they all link TARGET_SPARTA_BUILD_MPI to pick up shared compile
  # options. sfloat.h's AD branch needs Sacado's headers wherever that flag
  # is visible, so the include path has to ride along on the same target --
  # attaching Sacado::Sacado only to the main sparta lib (src/CMakeLists.txt)
  # is not enough; package targets never link the main lib itself.
  if(SPARTA_ENABLE_AD)
    set_property(
      TARGET ${TARGET_SPARTA_BUILD_MPI}
      APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${SACADO_INCLUDE_DIRS})
  endif()
endif()

if(SPARTA_CTEST_CONFIGS)
  string(REPLACE " " ";" SPARTA_CTEST_CONFIGS "${SPARTA_CTEST_CONFIGS}")
endif()
