cmake_minimum_required(VERSION 3.18)
project(Neko-TOP VERSION 0.0.1 LANGUAGES C CXX Fortran)

# ============================================================================ #
# Setup Compilers

# Set language standards
set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_Fortran_STANDARD 2008)

# ============================================================================ #
# Documentation
# ============================================================================ #

option(BUILD_DOCS "Build documentation" OFF)
option(ONLY_DOCS "Only build the documentation" OFF)
if (BUILD_DOCS OR ONLY_DOCS)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/documentation)

    if (ONLY_DOCS)
        return()
    endif()
endif()

# ============================================================================ #
# Setup External Libraries
# ............................................................................ #
# Check for the required libraries.
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# Check for the optional libraries.
find_package(MPI REQUIRED)

# If the user did not specify a build type, default to Release
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build.")
endif()

# Set default paths for the external dependencies
set(JSON_FORTRAN_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/json-fortran
    CACHE STRING "Path to the JSON-Fortran installation")
set(NEKO_DIR         ${CMAKE_CURRENT_SOURCE_DIR}/external/neko
    CACHE STRING "Path to the Neko installation")
set(PFUNIT_DIR       ${CMAKE_CURRENT_SOURCE_DIR}/external/pFUnit
    CACHE STRING "Path to the pFUnit installation")

# ............................................................................ #
# Setup of JSON-Fortran

# Check that the JSON_FORTRAN_DIR directories exist
if (NOT IS_DIRECTORY ${JSON_FORTRAN_DIR})
    message(FATAL_ERROR
        "JSON-Fortran directory not found (${JSON_FORTRAN_DIR})."
        "Please run `git submodule update --init` to download the external libraries."
    )
endif()

# Find the JSON-Fortran library and add it to the CMake search path
set(jsonfortran-gnu_DIR ${JSON_FORTRAN_DIR})
if (IS_DIRECTORY ${jsonfortran-gnu_DIR}/lib)
    set(jsonfortran-gnu_DIR ${jsonfortran-gnu_DIR}/lib)
elseif(IS_DIRECTORY ${jsonfortran-gnu_DIR}/lib64)
    set(jsonfortran-gnu_DIR ${jsonfortran-gnu_DIR}/lib64)
elseif(IS_DIRECTORY ${jsonfortran-gnu_DIR}/lib32)
    set(jsonfortran-gnu_DIR ${jsonfortran-gnu_DIR}/lib32)
endif()

set(jsonfortran-gnu_DIR ${jsonfortran-gnu_DIR}/cmake/jsonfortran-gnu-8.3.0)

# Check that the JSON-Fortran library was found
find_package(jsonfortran-gnu REQUIRED)

# ............................................................................ #
# Setup of the Neko library

# Check that the NEKO_DIR directories exist
if (NOT IS_DIRECTORY ${NEKO_DIR})
    message(FATAL_ERROR
        "Neko directory not found (${NEKO_DIR})."
        "Please run `git submodule update --init` to download the external libraries."
    )
endif()

# Set the proper paths for the NEKO_DIR and JSON_FORTRAN_DIR directories
set(ENV{PKG_CONFIG_PATH} ${NEKO_DIR}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH})

find_package(PkgConfig REQUIRED)
pkg_check_modules(neko REQUIRED IMPORTED_TARGET neko)

# ============================================================================ #
# Setup Build options
# ============================================================================ #
# Enable devices
include(CheckLanguage)

check_language(CUDA)
if (CMAKE_CUDA_COMPILER)
    enable_language(CUDA)
    add_compile_definitions(HAVE_CUDA)
endif()

# If we are building in Debug mode, enable -g and -Wall for all languages
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(
        $<$<NOT:$<COMPILE_LANGUAGE:CUDA>>:-Og>
        $<$<NOT:$<COMPILE_LANGUAGE:CUDA>>:-Wpedantic>
        -Wall -Wextra -Wno-unused-dummy-argument
    )
endif()

# Enforce compiler options
set(CMAKE_Fortran_PREPROCESS ON)

# ============================================================================ #
# Setup the Neko-TOP library
# ============================================================================ #

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/sources)

# ============================================================================ #
# Build the examples
# ============================================================================ #

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/examples)

# ============================================================================ #
# Build the tests
# ============================================================================ #

# Setup testing
option(BUILD_TESTING "Build the tests" OFF)

include(CTest)
if (BUILD_TESTING)
    find_package(PFUNIT REQUIRED)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/tests)
endif()

# ============================================================================ #