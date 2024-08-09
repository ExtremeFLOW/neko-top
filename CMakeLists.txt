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

cmake_minimum_required(VERSION 3.18)
project(Neko-TOP VERSION 0.0.1 LANGUAGES Fortran)

# Set language standards
set(CMAKE_Fortran_STANDARD 2008)

# Enable devices
option(DEVICE_TYPE "Enable device support" OFF)
if (DEVICE_TYPE STREQUAL "CUDA")
    include(CheckLanguage)
    check_language(CUDA)
    if (CMAKE_CUDA_COMPILER)
        enable_language(CUDA)
        add_compile_definitions(HAVE_CUDA)
    endif()
endif()

# ............................................................................ #
# Options

option(BUILD_DOCS "Build documentation" OFF)
option(BUILD_EXAMPLES "Build the examples" ON)
option(BUILD_TESTING "Build the tests" OFF)

# Advanced options
include(CMakeDependentOption)

option(ONLY_DOCS "Only build the documentation" OFF)
if (ONLY_DOCS)
    set(BUILD_DOCS ON)
    set(BUILD_EXAMPLES OFF)
    set(BUILD_TESTING OFF)
endif()

if (CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE "Release"
        CACHE STRING "Choose the type of build" FORCE)
endif()

# Set default paths for the external dependencies
set(JSON_FORTRAN_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/json-fortran
    CACHE STRING "Path to the JSON-Fortran installation")
set(NEKO_DIR         ${CMAKE_CURRENT_SOURCE_DIR}/external/neko
    CACHE STRING "Path to the Neko installation")
set(PFUNIT_DIR       ${CMAKE_CURRENT_SOURCE_DIR}/external/pFUnit
    CACHE STRING "Path to the pFUnit installation")

# ............................................................................ #
# Define compiler flags

set(CMAKE_Fortran_PREPROCESS ON)
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(
        $<$<COMPILE_LANGUAGE:Fortran>:-Og>
        $<$<COMPILE_LANGUAGE:Fortran>:-Wall>
        $<$<COMPILE_LANGUAGE:Fortran>:-Wextra>
        $<$<COMPILE_LANGUAGE:Fortran>:-Wpedantic>
        $<$<COMPILE_LANGUAGE:Fortran>:-Wno-unused-dummy-argument>
        $<$<COMPILE_LANGUAGE:Fortran>:-Werror>
    )
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    add_compile_options(
        $<$<COMPILE_LANGUAGE:Fortran>:-O3>
        $<$<COMPILE_LANGUAGE:Fortran>:-march=native>
        $<$<COMPILE_LANGUAGE:Fortran>:-mtune=native>
        $<$<COMPILE_LANGUAGE:Fortran>:-funroll-loops>
        $<$<COMPILE_LANGUAGE:Fortran>:-flto>
        $<$<COMPILE_LANGUAGE:Fortran>:-fwhole-program>
        $<$<COMPILE_LANGUAGE:Fortran>:-Werror>
    )
endif()

# ============================================================================ #
# Documentation
# ============================================================================ #

if (BUILD_DOCS)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/documentation)
endif()

if (ONLY_DOCS)
    return()
endif()

# ============================================================================ #
# Setup External Libraries
# ............................................................................ #
# Check for the required libraries.
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# Check for the optional libraries.
find_package(MPI REQUIRED)

# ............................................................................ #
# Setup of JSON-Fortran

# Set the proper paths for the JSON_FORTRAN_DIR directory
if (IS_DIRECTORY ${JSON_FORTRAN_DIR})
    set(CMAKE_PREFIX_PATH ${JSON_FORTRAN_DIR} ${CMAKE_PREFIX_PATH})
endif()

# Check that the JSON-Fortran library was found
find_package(jsonfortran-gnu REQUIRED)

# ............................................................................ #
# Setup of the Neko library

# Check that the NEKO_DIR directories exist
if (NOT IS_DIRECTORY ${NEKO_DIR})
    message(FATAL_ERROR
        "Neko directory not found (${NEKO_DIR})."
        "Please run `git submodule update --init` to download the external "
        "libraries."
    )
endif()

# Set the proper paths for the NEKO_DIR and JSON_FORTRAN_DIR directories
set(ENV{PKG_CONFIG_PATH} ${NEKO_DIR}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH})

find_package(PkgConfig REQUIRED)
pkg_check_modules(neko REQUIRED IMPORTED_TARGET neko)

# ============================================================================ #
# Setup the Neko-TOP library
# ============================================================================ #

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/sources)

# ============================================================================ #
# Build the examples
# ============================================================================ #

if (BUILD_EXAMPLES)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/examples)
endif()

# ============================================================================ #
# Build the tests
# ============================================================================ #

# Setup testing
include(CTest)
if (BUILD_TESTING)
    find_package(PFUNIT REQUIRED
        HINTS ${PFUNIT_DIR}
    )
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/tests)
endif()

# ============================================================================ #
# Write messages to the user about the options and paths if we are at the top
# level of the build.

if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)

    string(REPEAT "=" 80 BARRIER)

    message(STATUS "")
    message(STATUS ${BARRIER})
    message(STATUS "Neko-TOP: Neko Topology Optimization Package")
    message(STATUS "")
    message(STATUS "Current configuration:")
    message(STATUS "- CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
    message(STATUS "- BUILD_DOCS: ${BUILD_DOCS}")
    message(STATUS "- BUILD_EXAMPLES: ${BUILD_EXAMPLES}")
    message(STATUS "- BUILD_TESTING: ${BUILD_TESTING}")
    message(STATUS "- ONLY_DOCS: ${ONLY_DOCS}")
    message(STATUS "- DEVICE_TYPE: ${DEVICE_TYPE}")
    message(STATUS "")
    message(STATUS "External dependencies:")
    message(STATUS "- JSON_FORTRAN_DIR: ${JSON_FORTRAN_DIR}")
    message(STATUS "- NEKO_DIR: ${NEKO_DIR}")
    message(STATUS "- PFUNIT_DIR: ${PFUNIT_DIR}")
    message(STATUS ${BARRIER})
    message(STATUS "")

endif()

# ============================================================================ #
# End of file