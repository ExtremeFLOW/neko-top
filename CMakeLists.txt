# Copyright (c) 2023-2024, The Neko-TOP Authors
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
#   * Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer.
#
#   * Redistributions in binary form must reproduce the above
#     copyright notice, this list of conditions and the following
#     disclaimer in the documentation and/or other materials provided
#     with the distribution.
#
#   * Neither the name of the authors nor the names of its
#     contributors may be used to endorse or promote products derived
#     from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
# ============================================================================ #
# Neko-TOP: Neko Topology Optimization Package
# ============================================================================ #

cmake_minimum_required(VERSION 3.21)
project(Neko-TOP VERSION 0.0.1 LANGUAGES Fortran)

# Set language standards
set(CMAKE_Fortran_STANDARD 2008)

# ............................................................................ #
# Define cmake variables

# Define the list of available build types
set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build")
set(CMAKE_BUILD_TYPE_LIST
    Debug
    Release
    Testing
)
set(CMAKE_Fortran_PREPROCESS ON)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# ............................................................................ #
# Options

# Assign default values to options
option(BUILD_EXAMPLES "Build the examples" OFF)
option(BUILD_DOCS "Build documentation" OFF)
option(BUILD_TESTING "Build the tests" OFF)
option(ONLY_DOCS "Only build the documentation" OFF)

# ============================================================================ #
# Documentation
# ============================================================================ #

# Add the documentation target to ALL if we have "BUILD_DOCS" enabled
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/documentation)
if (ONLY_DOCS)
    return()
endif()

# ============================================================================ #
# Setup External Libraries
# ============================================================================ #
# Check for the required libraries.
find_package(BLAS QUIET)
find_package(LAPACK QUIET)

# Check for the optional libraries.
find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED COMPONENTS Fortran)

# We need PKG Config to find the JSON-Fortran and Neko libraries
find_package(PkgConfig REQUIRED)

# Check that the MPI Fortran 2008 module is available
if (NOT MPI_Fortran_HAVE_F08_MODULE)
    message(FATAL_ERROR "MPI Fortran 2008 module is required for Neko-TOP.")
endif()
include_directories(${MPI_Fortran_INCLUDE_DIRS})
include_directories(${MPI_Fortran_MODULE_DIR})

# Set the default external directory
if (DEFINED ENV{EXTERNAL_DIR})
    set(EXTERNAL_DIR $ENV{EXTERNAL_DIR})
else()
    set(EXTERNAL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external)
endif()

# ............................................................................ #
# Setup of JSON-Fortran

# Set the proper paths for the JSON_FORTRAN_DIR directory
if (DEFINED ENV{JSON_FORTRAN_DIR} AND IS_ABSOLUTE $ENV{JSON_FORTRAN_DIR})
    set(JSON_FORTRAN_DIR $ENV{JSON_FORTRAN_DIR} CACHE STRING
        "Path to the JSON-Fortran installation")
elseif(DEFINED ENV{JSON_FORTRAN_DIR})
    set(JSON_FORTRAN_DIR ${EXTERNAL_DIR}/$ENV{JSON_FORTRAN_DIR}
        CACHE STRING "Path to the JSON-Fortran installation")
else()
    set(JSON_FORTRAN_DIR ${EXTERNAL_DIR}/json-fortran
        CACHE STRING "Path to the JSON-Fortran installation")
endif()

set(CMAKE_PREFIX_PATH ${JSON_FORTRAN_DIR} ${CMAKE_PREFIX_PATH})

# Find the JSON-Fortran library
pkg_check_modules(json-fortran REQUIRED IMPORTED_TARGET json-fortran)

# ............................................................................ #
# Setup of the Neko library

# Set the proper paths for the NEKO_DIR
if (DEFINED ENV{NEKO_DIR} AND IS_ABSOLUTE $ENV{NEKO_DIR})
    set(NEKO_DIR $ENV{NEKO_DIR} CACHE STRING
        "Path to the Neko installation")
elseif(DEFINED ENV{NEKO_DIR})
    set(NEKO_DIR ${EXTERNAL_DIR}/$ENV{NEKO_DIR}
        CACHE STRING "Path to the Neko installation")
else()
    set(NEKO_DIR ${EXTERNAL_DIR}/neko
        CACHE STRING "Path to the Neko installation")
endif()

set(CMAKE_PREFIX_PATH ${NEKO_DIR} ${CMAKE_PREFIX_PATH})

# Find the Neko library
pkg_check_modules(neko REQUIRED IMPORTED_TARGET neko)
pkg_get_variable(DEVICE_TYPE neko backend)

# ............................................................................ #
# HDF5 Support
#
# Neko-TOP never links HDF5 itself. HAVE_HDF5 only gates which sources are
# compiled (see sources/IO/CMakeLists.txt) and which tests are registered;
# every HDF5 include and link flag already arrives transitively through
# PkgConfig::neko. Searching for HDF5 independently is therefore not merely
# redundant, it is harmful: Neko discovers HDF5 through pkg-config
# (m4/ax_hdf5.m4) or, on HPE Cray systems, through the compiler wrapper
# (m4/ax_cray.m4), and a CMake-side find_package() silently disagrees whenever
# those mechanisms disagree. On LUMI the cray-hdf5-parallel module ships a
# pkg-config file but no CMake package config, so find_package() ignored the
# loaded module and selected an unrelated system-wide installation instead,
# leaving Neko and Neko-TOP built against different HDF5 libraries.
#
# So we do not search. We ask the Neko library we are actually linking against
# what it uses, then confirm we can compile and link against the same thing.

# HDF5_DIR is retained for reporting only: it names the installation that
# scripts/dependencies.sh:find_hdf5 provisioned and that Neko was configured
# `--with-hdf5=` against. It is deliberately not added to CMAKE_PREFIX_PATH --
# HDF5_DIR is CMake's reserved config-mode variable, and find_package() used to
# overwrite it with the directory holding hdf5-config.cmake.
if (NOT DEFINED HDF5_DIR)
    if (DEFINED ENV{HDF5_DIR} AND IS_ABSOLUTE $ENV{HDF5_DIR})
        set(HDF5_DIR $ENV{HDF5_DIR} CACHE STRING
            "Path to the HDF5 installation")
    elseif(DEFINED ENV{HDF5_DIR})
        set(HDF5_DIR ${EXTERNAL_DIR}/$ENV{HDF5_DIR}
            CACHE STRING "Path to the HDF5 installation")
    endif()
endif()

# Determine which flavour of HDF5, if any, the Neko library links against.
#
# Neko compiles src/io/hdf5_session.F90, hdf5_file.F90 and vtkhdf_file.F90
# unconditionally and guards their bodies with `#ifdef HAVE_HDF5`, so libneko
# carries undefined HDF5 symbols if and only if Neko was built with HDF5. Neko
# calls the parallel-only h5pset_fapl_mpio/h5pset_dxpl_mpio routines itself, so
# one scan answers both "has HDF5" and "is it parallel". Symbols are matched on
# the bare routine name rather than on any compiler's name mangling, so this
# behaves identically for GNU and Cray-compiled archives.
#
# Sets in the caller's scope:
#   NEKO_HDF5_STATUS  PARALLEL, SERIAL, PRESENT, NONE or UNKNOWN
#   NEKO_HDF5_ORIGIN  how that conclusion was reached, for diagnostics
function(neko_top_probe_neko_hdf5)
    set(status "UNKNOWN")
    set(origin "no usable evidence")

    find_program(NEKO_TOP_NM NAMES nm llvm-nm)
    mark_as_advanced(NEKO_TOP_NM)

    # Static archive only, matching what scripts/dependencies.sh:find_neko
    # looks for. A shared libneko would need `nm -D` and could be stripped, so
    # rather than risk a wrong answer we fall through to the flag fallback.
    file(GLOB neko_archives "${NEKO_DIR}/lib*/libneko.a")

    if (NEKO_TOP_NM AND neko_archives)
        list(GET neko_archives 0 archive)
        execute_process(
            COMMAND ${NEKO_TOP_NM} --undefined-only "${archive}"
            OUTPUT_VARIABLE symbols
            ERROR_VARIABLE nm_stderr
            RESULT_VARIABLE nm_status
            OUTPUT_STRIP_TRAILING_WHITESPACE)

        if (nm_status EQUAL 0)
            string(TOLOWER "${symbols}" symbols)
            if (symbols MATCHES "h5pset_fapl_mpio")
                set(status "PARALLEL")
            elseif (symbols MATCHES "h5open_f")
                set(status "SERIAL")
            else()
                set(status "NONE")
            endif()
            set(origin "undefined symbols in ${archive}")
        endif()
    endif()

    # Fallback for the case where nm is unavailable or the archive could not be
    # read. The flags Neko publishes reveal whether it uses HDF5, but not which
    # flavour, hence PRESENT. Note that this cannot reach a NONE verdict: on
    # HPE Cray systems the ftn wrapper supplies HDF5 implicitly and neko.pc
    # mentions it nowhere, so an absence of flags proves nothing.
    if (status STREQUAL "UNKNOWN")
        string(TOLOWER "${neko_CFLAGS};${neko_LDFLAGS}" neko_flags)
        if (neko_flags MATCHES "hdf5")
            set(status "PRESENT")
            set(origin "HDF5 flags published in neko.pc")
        endif()
    endif()

    set(NEKO_HDF5_STATUS ${status} PARENT_SCOPE)
    set(NEKO_HDF5_ORIGIN ${origin} PARENT_SCOPE)
endfunction()

# Confirm that Neko-TOP can actually compile and link the parallel HDF5 calls
# its own I/O sources make, using nothing but the flags Neko publishes plus
# MPI. Mirrors the real usage in sources/IO/hdf5/design_hdf5_io.f90. Because
# check_fortran_source_compiles() builds an executable, this exercises exactly
# the link step that used to fail with `undefined reference to
# h5pset_fapl_mpio_f_` long after configure had reported success.
function(neko_top_check_hdf5_usable result)
    set(CMAKE_REQUIRED_QUIET TRUE)
    set(CMAKE_REQUIRED_INCLUDES
        ${neko_INCLUDE_DIRS}
        ${MPI_Fortran_INCLUDE_DIRS}
        ${MPI_Fortran_MODULE_DIR})

    if (neko_LINK_LIBRARIES)
        set(CMAKE_REQUIRED_LIBRARIES ${neko_LINK_LIBRARIES})
    else()
        set(CMAKE_REQUIRED_LIBRARIES ${neko_LDFLAGS})
    endif()
    list(APPEND CMAKE_REQUIRED_LIBRARIES ${MPI_Fortran_LIBRARIES})

    check_fortran_source_compiles([==[
program neko_top_hdf5_probe
  use hdf5
  use mpi_f08, only: MPI_COMM_WORLD, MPI_INFO_NULL
  implicit none
  integer :: ierr
  integer(hid_t) :: fapl_id, xfer_id
  call h5open_f(ierr)
  call h5pcreate_f(H5P_FILE_ACCESS_F, fapl_id, ierr)
  call h5pset_fapl_mpio_f(fapl_id, MPI_COMM_WORLD%mpi_val, &
       MPI_INFO_NULL%mpi_val, ierr)
  call h5pcreate_f(H5P_DATASET_XFER_F, xfer_id, ierr)
  call h5pset_dxpl_mpio_f(xfer_id, H5FD_MPIO_COLLECTIVE_F, ierr)
  call h5close_f(ierr)
end program neko_top_hdf5_probe
]==] ${result} SRC_EXT f90)

    set(${result} ${${result}} PARENT_SCOPE)
endfunction()

include(CheckFortranSourceCompiles)
neko_top_probe_neko_hdf5()

set(HAVE_HDF5 FALSE)
if (NEKO_HDF5_STATUS STREQUAL "NONE")

    # Neko has no HDF5, so neither may we. Compiling our HDF5 I/O against some
    # other HDF5 the compiler happens to find is what produces two HDF5
    # libraries in one executable, and with it the teardown crashes seen
    # previously. Disable quietly; a build without HDF5 is a valid build.
    message(STATUS "Neko links no HDF5, disabling Neko-TOP HDF5 support "
        "(determined from ${NEKO_HDF5_ORIGIN})")

elseif (NEKO_HDF5_STATUS STREQUAL "SERIAL")

    # Defensive: no Neko built today can reach here, because everything Neko
    # compiles under `#ifdef HAVE_HDF5` references h5pset_fapl_mpio_f, so an
    # HDF5-enabled Neko always classifies as PARALLEL. The branch exists so
    # that a future Neko growing a serial HDF5 path fails here, at configure
    # time, rather than at link time in every Neko-TOP executable.
    message(FATAL_ERROR
        "Neko links a serial HDF5, but Neko-TOP's HDF5 I/O uses the "
        "MPI-parallel interface (h5pset_fapl_mpio_f, h5pset_dxpl_mpio_f) "
        "unconditionally.\n"
        "  Determined from: ${NEKO_HDF5_ORIGIN}\n"
        "Rebuild Neko against an MPI-parallel HDF5, or rebuild it without "
        "HDF5 to build Neko-TOP with HDF5 support disabled.")

else()

    neko_top_check_hdf5_usable(NEKO_TOP_HDF5_USABLE)

    if (NEKO_TOP_HDF5_USABLE)
        set(HAVE_HDF5 TRUE)
        add_compile_definitions(HAVE_HDF5=1)
    elseif (NEKO_HDF5_STATUS STREQUAL "UNKNOWN")
        message(WARNING
            "Could not determine whether Neko links HDF5 "
            "(${NEKO_HDF5_ORIGIN}), and Neko-TOP cannot compile against a "
            "parallel HDF5 either. Continuing with HDF5 support disabled.")
    else()
        message(FATAL_ERROR
            "Neko links HDF5 but Neko-TOP cannot compile and link against "
            "it, so the two would disagree about HDF5 support.\n"
            "  Neko HDF5:       ${NEKO_HDF5_STATUS} "
            "(from ${NEKO_HDF5_ORIGIN})\n"
            "  Fortran compiler: ${CMAKE_Fortran_COMPILER}\n"
            "See the check_fortran_source_compiles entry in "
            "${CMAKE_BINARY_DIR}/CMakeFiles/CMakeConfigureLog.yaml for the "
            "compiler output. Build Neko and Neko-TOP in the same "
            "environment, with the same HDF5 module loaded.")
    endif()

endif()

# ............................................................................ #
# Setup PFUnit

if (DEFINED ENV{PFUNIT_DIR} AND IS_ABSOLUTE $ENV{PFUNIT_DIR})
    set(PFUNIT_DIR $ENV{PFUNIT_DIR} CACHE STRING
        "Path to the pFUnit installation")
elseif(DEFINED ENV{PFUNIT_DIR})
    set(PFUNIT_DIR ${EXTERNAL_DIR}/$ENV{PFUNIT_DIR}
        CACHE STRING "Path to the pFUnit installation")
else()
    set(PFUNIT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/pfunit
        CACHE STRING "Path to the pFUnit installation")
endif()

# ============================================================================ #
# Setup compiler specific options
# ============================================================================ #

# Set the compiler flags
include(${CMAKE_Fortran_COMPILER_ID} OPTIONAL RESULT_VARIABLE compiler_flags_found)
if (NOT compiler_flags_found)
    message(WARNING "No specific compiler flags found for "
        "${CMAKE_Fortran_COMPILER_ID}.")
endif()

# ============================================================================ #
# Setup the Neko-TOP library
# ============================================================================ #

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/sources)

# Add a static library target for use in examples and tests
# This ensures the complete Neko-TOP library is built and linked before any
# examples or tests are built.
add_library(neko-top-a STATIC IMPORTED)
set_target_properties(neko-top-a PROPERTIES
    IMPORTED_LOCATION "${CMAKE_BINARY_DIR}/libneko-top.a"
)
add_custom_target(neko-top_ready
    DEPENDS neko-top
    COMMENT "Neko-TOP library is ready.")
add_dependencies(neko-top-a neko-top_ready)

# ============================================================================ #
# Build the examples
# ============================================================================ #

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/examples)

# ============================================================================ #
# Build the tests
# ============================================================================ #

# Setup testing
include(CTest)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/tests)

# ============================================================================ #
# Write messages to the user about the options and paths if we are at the top
# level of the build.

if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)

    string(REPEAT "=" 80 BARRIER)

    message(STATUS ${BARRIER})
    message(STATUS "Neko-TOP: Neko Topology Optimization Package")
    message(STATUS "")
    message(STATUS "Current configuration:")
    message(STATUS "- CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
    message(STATUS "- DEVICE_TYPE:      ${DEVICE_TYPE}")
    message(STATUS "- BUILD_DOCS:       ${BUILD_DOCS}")
    message(STATUS "- BUILD_TESTING:    ${BUILD_TESTING}")
    message(STATUS "- BUILD_EXAMPLES:   ${BUILD_EXAMPLES}")
    message(STATUS "")
    message(STATUS "External dependencies:")
    message(STATUS "- JSON_FORTRAN_DIR: ${JSON_FORTRAN_DIR}")
    message(STATUS "- NEKO_DIR:         ${NEKO_DIR}")
    message(STATUS "- PFUNIT_DIR:       ${PFUNIT_DIR}")
    if (HAVE_HDF5)
        message(STATUS "- HDF5:             Enabled, following Neko "
            "(${NEKO_HDF5_STATUS})")
    else()
        message(STATUS "- HDF5:             Disabled, Neko links no HDF5")
    endif()
    if (DEFINED HDF5_DIR)
        message(STATUS "- HDF5_DIR:         ${HDF5_DIR}")
    endif()
    message(STATUS ${BARRIER})
    message(STATUS "")

endif()

# ============================================================================ #
# End of file
