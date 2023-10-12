cmake_minimum_required(VERSION 3.0.0)
project(Neko-TOP VERSION 0.1.0 LANGUAGES C CXX Fortran)

# ============================================================================ #
# Setup CMake

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

# ============================================================================ #
# Setup Compilers

# Set language standards
set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_Fortran_STANDARD 2008)

# ============================================================================ #
# Setup External Libraries

set(EXTERNAL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/External)
set(NEKO_DIR ${EXTERNAL_DIR}/neko)

# Extend the search path for pkg-config files to include the neko installation.
cmake_path(SET pkgconfig_path ${NEKO_DIR} NORMALIZE)
set(ENV{PKG_CONFIG_PATH} ${pkgconfig_path})

# Find the JSON-Fortran library
set(jsonfortran-gnu_DIR ${EXTERNAL_DIR}/json-fortran)
if (IS_DIRECTORY ${jsonfortran-gnu_DIR}/lib)
    set(jsonfortran-gnu_DIR ${jsonfortran-gnu_DIR}/lib)
elseif(IS_DIRECTORY ${jsonfortran-gnu_DIR}/lib64)
    set(jsonfortran-gnu_DIR ${jsonfortran-gnu_DIR}/lib64)
elseif(IS_DIRECTORY ${jsonfortran-gnu_DIR}/lib32)
    set(jsonfortran-gnu_DIR ${jsonfortran-gnu_DIR}/lib32)
endif()
set(jsonfortran-gnu_DIR ${jsonfortran-gnu_DIR}/cmake/jsonfortran-gnu-8.3.0)

# ============================================================================ #
# Add examples

add_subdirectory(Examples)
