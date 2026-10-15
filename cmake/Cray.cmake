# ============================================================================ #
# This file contains the compiler flags for the project that are not purely a
# function of CMAKE_BUILD_TYPE.
#
# Per-build-type Fortran flags (CMAKE_Fortran_FLAGS_DEBUG, _TESTING and
# _RELEASE) are set in cmake/FlagsInit.cmake, through CMake's
# CMAKE_USER_MAKE_RULES_OVERRIDE_Fortran hook. What is left here is the HIP
# device-compile options, which apply identically to Release and Testing.
# Debug carries no HIP flags.
#
# NOTE: not verified in this environment -- there is no Cray compiler here.
set(_neko_top_hip_release_testing
    $<AND:$<COMPILE_LANGUAGE:HIP>,$<CONFIG:Release,Testing>>)
add_compile_options(
    $<${_neko_top_hip_release_testing}:-w>  # Suppress all warnings
    $<${_neko_top_hip_release_testing}:-O3> # Optimize to level 3
    $<${_neko_top_hip_release_testing}:--offload-arch=gfx90a> # HIP target arch
)
unset(_neko_top_hip_release_testing)
