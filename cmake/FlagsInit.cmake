# ============================================================================ #
# Per-build-type default Fortran flags, keyed by CMAKE_Fortran_COMPILER_ID.
#
# This file is loaded through CMAKE_USER_MAKE_RULES_OVERRIDE_Fortran (set in
# CMakeLists.txt, before `project()`). CMake includes it after the compiler
# information files have set their own CMAKE_Fortran_FLAGS_<CONFIG>_INIT
# defaults, but before those are cached into CMAKE_Fortran_FLAGS_<CONFIG>, so
# setting the *_INIT variables here is the documented way to replace CMake's
# per-config defaults instead of patching the cached flags afterwards. See:
# https://cmake.org/cmake/help/latest/variable/CMAKE_USER_MAKE_RULES_OVERRIDE.html
#
# Flags that are not a function of build type (warnings, the
# NEKO_TOP_RUNTIME_CHECKS option) stay out of this file; they are applied in
# cmake/<COMPILER_ID>.cmake through add_compile_options()/generator
# expressions after `project()`, where ordinary (non-_INIT) variables and
# generator expressions are available.

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_Fortran_FLAGS_DEBUG_INIT   "-g -Og")
    set(CMAKE_Fortran_FLAGS_TESTING_INIT "-g -O2")
    set(CMAKE_Fortran_FLAGS_RELEASE_INIT "-O3 -funroll-loops")

elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
    set(CMAKE_Fortran_FLAGS_DEBUG_INIT   "-g -O0 -m1 -RABCDS")
    set(CMAKE_Fortran_FLAGS_TESTING_INIT "-g -O2")
    set(CMAKE_Fortran_FLAGS_RELEASE_INIT "-O2 -m4 -e0")

else()
    # Unknown compiler: only seed Testing, which otherwise would not exist at
    # all (it is not one of CMake's four built-in configurations). Seed it
    # from RelWithDebInfo -- itself an optimised-plus-debug-info config, the
    # closest built-in match to what Testing is for -- rather than leaving it
    # empty. Leave Debug/Release at whatever CMake's own compiler-info file
    # already set for this compiler.
    set(CMAKE_Fortran_FLAGS_TESTING_INIT
        "${CMAKE_Fortran_FLAGS_RELWITHDEBINFO_INIT}")
endif()
