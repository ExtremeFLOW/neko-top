cmake_minimum_required(VERSION 3.18)
project(Neko-TOP VERSION 0.0.1 LANGUAGES C CXX Fortran)

# ============================================================================ #
# Setup Compilers

# Set language standards
set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_Fortran_STANDARD 2008)

# ============================================================================ #
# Setup External Libraries
# ............................................................................ #
# Check for the required libraries.
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# Check for the optional libraries.
find_package(MPI QUIET)

# If the user did not specify a build type, default to Release
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
endif()

# ============================================================================ #
# Enable devices
include(CheckLanguage)

# Enable CUDA if the compiler is found
check_language(CUDA)
if (CMAKE_CUDA_COMPILER)
    enable_language(CUDA)
    add_compile_definitions(HAVE_CUDA)
endif()

# ............................................................................ #
# Setup of JSON-Fortran

# If the user did not specify JSONFORTRAN_DIR, default to the
# external/json-fortran directory
set(JSONFORTRAN_DIR
    ${CMAKE_CURRENT_SOURCE_DIR}/external/json-fortran
    CACHE PATH
    "Path to the JSON-Fortran installation"
)

# Check that the JSONFORTRAN_DIR directories exist
if (NOT IS_DIRECTORY ${JSONFORTRAN_DIR})
    message(FATAL_ERROR
        "JSON-Fortran directory not found (${JSONFORTRAN_DIR})."
        "Please run `git submodule update --init` to download the external libraries."
    )
endif()

# Find the JSON-Fortran library and add it to the CMake search path
set(jsonfortran-gnu_DIR ${JSONFORTRAN_DIR})
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

# If the user did not specify NEKO_DIR, default to the external/neko directory
set(NEKO_DIR ${CMAKE_SOURCE_DIR}/external/neko CACHE PATH
    "Path to the Neko installation")

# Check that the NEKO_DIR directories exist
if (NOT IS_DIRECTORY ${NEKO_DIR})
    message(FATAL_ERROR
        "Neko directory not found (${NEKO_DIR})."
        "Please run `git submodule update --init` to download the external libraries."
    )
endif()

# Set the proper paths for the NEKO_DIR and JSONFORTRAN_DIR directories
set(ENV{PKG_CONFIG_PATH} ${NEKO_DIR})

find_package(PkgConfig REQUIRED)
pkg_check_modules(neko REQUIRED IMPORTED_TARGET neko)

# ............................................................................ #
# Setup of the pFUnit library

find_package(pFUnit QUIET)
if (NOT PFUNIT_FOUND)
    include(FetchContent)
    message(STATUS "pFUnit not found. Fetching pFUnit from GitHub")
    FetchContent_Declare(
        pFUnit
        GIT_REPOSITORY https://github.com/Goddard-Fortran-Ecosystem/pFUnit.git
        GIT_TAG        v4.8.0
    )
    FetchContent_MakeAvailable(pFUnit)
endif()

# ============================================================================ #
# Setup Build options
# ============================================================================ #

# If we are building in Debug mode, enable -g and -Wall for all languages
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(
        $<$<NOT:$<COMPILE_LANGUAGE:CUDA>>:-Og>
        $<$<NOT:$<COMPILE_LANGUAGE:CUDA>>:-Wpedantic>
        -Wall -Wextra -Wno-unused-dummy-argument
        )

    set(BUILD_TESTING ON CACHE BOOL "Build the tests" FORCE)
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    set(BUILD_TESTING OFF CACHE BOOL "Build the tests" FORCE)
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

include(CTest)
if (BUILD_TESTING)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/tests)
endif()

# ============================================================================ #
# Documentation
# ============================================================================ #

option(BUILD_DOCS "Build documentation" ON)
if (BUILD_DOCS)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/documentation)
endif()

# ============================================================================ #