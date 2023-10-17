cmake_minimum_required(VERSION 3.0.0)
project(Neko-TOP VERSION 0.1.0 LANGUAGES C CXX Fortran)

# ============================================================================ #
# Setup CMake

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

# If the user did not specify a build type, default to Release
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# ============================================================================ #
# Setup Compilers

# Set language standards
set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_Fortran_STANDARD 2008)

# ============================================================================ #
# Setup External Libraries
# ............................................................................ #
# Setup of JSON-Fortran

# If the user did not specify JSONFORTRAN_DIR, default to the
# External/json-fortran directory
set(JSONFORTRAN_DIR ${CMAKE_CURRENT_SOURCE_DIR}/External/json-fortran CACHE PATH "Path to the JSON-Fortran installation")

# Check that the JSONFORTRAN_DIR directories exist
if (NOT IS_DIRECTORY ${JSONFORTRAN_DIR})
    message(FATAL_ERROR "JSON-Fortran directory not found (${JSONFORTRAN_DIR}). Please run `git submodule update --init` to download the external libraries.")
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

# If the user did not specify NEKO_DIR, default to the External/neko directory
set(NEKO_DIR ${CMAKE_CURRENT_SOURCE_DIR}/External/neko CACHE PATH "Path to the Neko installation")

# Check that the NEKO_DIR directories exist
if (NOT IS_DIRECTORY ${NEKO_DIR})
    message(FATAL_ERROR "Neko directory not found (${JSONFORTRAN_DIR}). Please run `git submodule update --init` to download the external libraries.")
endif()

# Set the proper paths for the NEKO_DIR and JSONFORTRAN_DIR directories
set(ENV{PKG_CONFIG_PATH} ${NEKO_DIR})

find_package(PkgConfig REQUIRED)
pkg_check_modules(neko REQUIRED IMPORTED_TARGET neko)

# ============================================================================ #
# Add examples

add_subdirectory(Examples)