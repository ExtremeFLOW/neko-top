# ============================================================================ #
# This file contains the compiler flags for the project that are not purely a
# function of CMAKE_BUILD_TYPE.
#
# Per-build-type optimisation/debug-info flags (CMAKE_Fortran_FLAGS_DEBUG,
# _TESTING and _RELEASE) are set in cmake/FlagsInit.cmake, through CMake's
# CMAKE_USER_MAKE_RULES_OVERRIDE_Fortran hook, so that they are the same kind
# of overridable per-config cache variable as any of CMake's own build types.
# What is left here is flags that cut across build types instead of selecting
# one: warnings, shared by Debug and Testing, and the NEKO_TOP_RUNTIME_CHECKS
# option, which is orthogonal to build type entirely.

# ............................................................................ #
# Warnings for Debug and Testing.
#
# Warnings are not errors by default. CI's Debug lane (check_compile.yml)
# turns them into errors with -DCMAKE_COMPILE_WARNING_AS_ERROR=ON, and a local
# build can do the same. Release gets no warning flags.
#
# -Warray-temporaries and -Wimplicit-interface are deliberately not in this
# set: array temporaries are a performance concern rather than a correctness
# one, and the existing sources trigger the warning at dozens of sites. Add
# either through CMAKE_Fortran_FLAGS to audit for them.
set(_neko_top_fortran_debug_testing
    $<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug,Testing>>)
add_compile_options(
    $<${_neko_top_fortran_debug_testing}:-Wall>
    $<${_neko_top_fortran_debug_testing}:-Wextra>
    $<${_neko_top_fortran_debug_testing}:-Wno-unused-dummy-argument>
    $<${_neko_top_fortran_debug_testing}:-Wno-function-elimination>
    $<${_neko_top_fortran_debug_testing}:-Wno-missing-include-dirs>
    $<${_neko_top_fortran_debug_testing}:-Wcompare-reals>
    $<${_neko_top_fortran_debug_testing}:-Wconversion>
    $<${_neko_top_fortran_debug_testing}:-pedantic-errors>
)
unset(_neko_top_fortran_debug_testing)

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
