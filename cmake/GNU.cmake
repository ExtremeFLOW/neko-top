# ============================================================================ #
# This file contains the compiler flags for the project.
#
# This allow us to set the flags for different compilers and different build
# types. The flags are set using the `add_compile_options` command.
#
# ============================================================================ #
# Set the compiler flags for GNU compilers
# ============================================================================ #

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(
        # Debugging flags
        $<$<COMPILE_LANGUAGE:Fortran>:-g>
        $<$<COMPILE_LANGUAGE:Fortran>:-Og>

        # Extra warnings for debug builds, but not errors
        $<$<COMPILE_LANGUAGE:Fortran>:-Wall>
        $<$<COMPILE_LANGUAGE:Fortran>:-Wextra>
        $<$<COMPILE_LANGUAGE:Fortran>:-Wno-unused-dummy-argument>
        $<$<COMPILE_LANGUAGE:Fortran>:-Wno-function-elimination>
        $<$<COMPILE_LANGUAGE:Fortran>:-Wno-missing-include-dirs>
        $<$<COMPILE_LANGUAGE:Fortran>:-Wcompare-reals>
        $<$<COMPILE_LANGUAGE:Fortran>:-Wconversion>
        $<$<COMPILE_LANGUAGE:Fortran>:-pedantic-errors>
    )
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    add_compile_options(
        $<$<COMPILE_LANGUAGE:Fortran>:-w>
        $<$<COMPILE_LANGUAGE:Fortran>:-O3>
        $<$<COMPILE_LANGUAGE:Fortran>:-funroll-loops>
    )
elseif(CMAKE_BUILD_TYPE STREQUAL "Testing")
    add_compile_options(
        # Debugging flags
        $<$<COMPILE_LANGUAGE:Fortran>:-g>
        $<$<COMPILE_LANGUAGE:Fortran>:-O2>

        # Extra warnings for testing builds, but not errors
        $<$<COMPILE_LANGUAGE:Fortran>:-Wall>
        $<$<COMPILE_LANGUAGE:Fortran>:-Wextra>
        $<$<COMPILE_LANGUAGE:Fortran>:-Wno-unused-dummy-argument>
        $<$<COMPILE_LANGUAGE:Fortran>:-Wno-function-elimination>
        $<$<COMPILE_LANGUAGE:Fortran>:-Wno-missing-include-dirs>
        $<$<COMPILE_LANGUAGE:Fortran>:-Wcompare-reals>
        $<$<COMPILE_LANGUAGE:Fortran>:-Wconversion>
        $<$<COMPILE_LANGUAGE:Fortran>:-pedantic-errors>
    )
endif()

# ============================================================================ #
# Define the runtime checks for GNU compilers
# ============================================================================ #
#
# Runtime checks:
# - bounds/do/pointer/recursion checking
# - FPE traps, and
# - poisoning values for uninitialised variables.

if (NEKO_TOP_RUNTIME_CHECKS STREQUAL "ON")
    add_compile_options(
        $<$<COMPILE_LANGUAGE:Fortran>:-fcheck=all>
        $<$<COMPILE_LANGUAGE:Fortran>:-ffpe-trap=invalid,zero,overflow>
        $<$<COMPILE_LANGUAGE:Fortran>:-finit-real=snan>
        $<$<COMPILE_LANGUAGE:Fortran>:-finit-integer=-99999999>
        $<$<COMPILE_LANGUAGE:Fortran>:-finit-character=1>
        $<$<COMPILE_LANGUAGE:Fortran>:-finit-logical=true>
        $<$<COMPILE_LANGUAGE:Fortran>:-finit-derived>
    )
elseif (NEKO_TOP_RUNTIME_CHECKS STREQUAL "AUTO"
        AND CMAKE_BUILD_TYPE STREQUAL "Debug")
     add_compile_options(
        $<$<COMPILE_LANGUAGE:Fortran>:-fcheck=all,no-array-temps,no-bounds>
        $<$<COMPILE_LANGUAGE:Fortran>:-ffpe-trap=invalid,zero,overflow>
        $<$<COMPILE_LANGUAGE:Fortran>:-finit-real=snan>
        $<$<COMPILE_LANGUAGE:Fortran>:-finit-integer=-99999999>
        $<$<COMPILE_LANGUAGE:Fortran>:-finit-character=1>
        $<$<COMPILE_LANGUAGE:Fortran>:-finit-logical=true>
        $<$<COMPILE_LANGUAGE:Fortran>:-finit-derived>
    )
elseif (NEKO_TOP_RUNTIME_CHECKS STREQUAL "AUTO"
        AND NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    # No runtime checks for Release and Testing builds
elseif (NOT NEKO_TOP_RUNTIME_CHECKS STREQUAL "OFF")
    message(FATAL_ERROR
        "Unknown NEKO_TOP_RUNTIME_CHECKS: '${NEKO_TOP_RUNTIME_CHECKS}'. "
        "Valid options are AUTO, ON or OFF.")
endif()
