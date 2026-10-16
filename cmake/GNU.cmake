# ============================================================================ #
# This file contains the compiler flags for the project.
#
# This allow us to set the flags for different compilers and different build
# types. The flags are set using the `add_compile_options` command.

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

# ............................................................................ #
# Runtime checks: bounds/do/pointer/recursion checking, FPE traps, and
# poisoning values for uninitialised variables. Orthogonal to build type --
# NEKO_TOP_RUNTIME_CHECKS (declared in CMakeLists.txt) is AUTO, ON or OFF, and
# defaults to OFF.
#
# There is deliberately no -finit-local-zero: it sets the integer, real,
# logical and character initialisers all to zero, so placed after the
# -finit-* options below it would silently replace their poison values with
# zero (gfortran's options.cc), defeating the point of poisoning. There is
# also no -fbacktrace: it is gfortran's default already.
if (NEKO_TOP_RUNTIME_CHECKS STREQUAL "ON")
    set(_neko_top_apply_runtime_checks TRUE)
elseif (NEKO_TOP_RUNTIME_CHECKS STREQUAL "AUTO")
    # Reserved to mean "Debug builds only". Until the Nek-style collapsed
    # indexing in the sources (x(i,1,1,1) with i running over the whole
    # array) is rewritten, -fcheck=bounds rejects it, so AUTO behaves as OFF
    # for now and is safe to leave at its default.
    set(_neko_top_apply_runtime_checks FALSE)
elseif (NEKO_TOP_RUNTIME_CHECKS STREQUAL "OFF")
    set(_neko_top_apply_runtime_checks FALSE)
else()
    message(FATAL_ERROR
        "Unknown NEKO_TOP_RUNTIME_CHECKS: '${NEKO_TOP_RUNTIME_CHECKS}'. "
        "Valid options are AUTO, ON or OFF.")
endif()

if (_neko_top_apply_runtime_checks)
    add_compile_options(
        $<$<COMPILE_LANGUAGE:Fortran>:-fcheck=all,no-array-temps>
        $<$<COMPILE_LANGUAGE:Fortran>:-ffpe-trap=invalid,zero,overflow>
        $<$<COMPILE_LANGUAGE:Fortran>:-finit-real=snan>
        $<$<COMPILE_LANGUAGE:Fortran>:-finit-integer=-99999999>
        $<$<COMPILE_LANGUAGE:Fortran>:-finit-character=1>
        $<$<COMPILE_LANGUAGE:Fortran>:-finit-logical=true>
        $<$<COMPILE_LANGUAGE:Fortran>:-finit-derived>
    )
endif()
unset(_neko_top_apply_runtime_checks)
