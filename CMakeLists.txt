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

if (NOT DEFINED HDF5_DIR)
    if (DEFINED ENV{HDF5_DIR} AND IS_ABSOLUTE $ENV{HDF5_DIR})
        set(HDF5_DIR $ENV{HDF5_DIR} CACHE STRING
            "Path to the HDF5 installation")
    elseif(DEFINED ENV{HDF5_DIR})
        set(HDF5_DIR ${EXTERNAL_DIR}/$ENV{HDF5_DIR}
            CACHE STRING "Path to the HDF5 installation")
    endif()
endif()

if (DEFINED HDF5_DIR)
    set(CMAKE_PREFIX_PATH ${HDF5_DIR} ${CMAKE_PREFIX_PATH})
endif()

find_package(HDF5 QUIET COMPONENTS Fortran)
if (HDF5_FOUND AND NOT HDF5_IS_PARALLEL)
    message(STATUS "Ignoring HDF5 without MPI support at "
        "${HDF5_INCLUDE_DIRS}")
    set(HDF5_FOUND FALSE)
endif()

if (HDF5_FOUND)
    set(HAVE_HDF5 TRUE)
    add_compile_definitions(HAVE_HDF5=1)
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
    if (DEFINED HDF5_DIR AND IS_DIRECTORY "${HDF5_DIR}")
        message(STATUS "- HDF5:             ${HDF5_DIR}")
    elseif (HDF5_FOUND)
        message(STATUS "- HDF5:             Found (HDF5 support enabled)")
    else()
        message(STATUS "- HDF5:             Not found (HDF5 support disabled)")
    endif()
    message(STATUS ${BARRIER})
    message(STATUS "")

endif()

# ============================================================================ #
# End of file
