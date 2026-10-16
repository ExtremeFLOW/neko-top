# ============================================================================ #
# This file contains the Cray compiler flags for the project.
#
# This allow us to set the flags for different compilers and different build
# types. The flags are set using the `add_compile_options` command.
#
# ============================================================================ #
# Set the compiler flags for Cray compilers
# ============================================================================ #

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(
        # Debugging flags
        $<$<COMPILE_LANGUAGE:Fortran>:-g>  # Enable debugging
        $<$<COMPILE_LANGUAGE:Fortran>:-O0> # Disable optimization
        $<$<COMPILE_LANGUAGE:Fortran>:-m1> # Set message level to verbose
        $<$<COMPILE_LANGUAGE:Fortran>:-ei> # Set uninitialized values to NaN

        $<$<COMPILE_LANGUAGE:HIP>:-g>      # Enable debugging
        $<$<COMPILE_LANGUAGE:HIP>:-O0>     # Optimize to level 0
        $<$<COMPILE_LANGUAGE:HIP>:-Wall>   # Suppress all warnings
        $<$<COMPILE_LANGUAGE:HIP>:--offload-arch=gfx90a> # Set the target architecture for HIP
    )
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    add_compile_options(
        $<$<COMPILE_LANGUAGE:Fortran>:-e0> # Initialize to 0
        $<$<COMPILE_LANGUAGE:Fortran>:-O2> # Optimize to level 2
        $<$<COMPILE_LANGUAGE:Fortran>:-m4> # Set message level to Error only

        $<$<COMPILE_LANGUAGE:HIP>:-w>      # Suppress all warnings
        $<$<COMPILE_LANGUAGE:HIP>:-O3>     # Optimize to level 3
        $<$<COMPILE_LANGUAGE:HIP>:--offload-arch=gfx90a> # Set the target architecture for HIP
    )
elseif(CMAKE_BUILD_TYPE STREQUAL "Testing")
    add_compile_options(
        $<$<COMPILE_LANGUAGE:Fortran>:-g>  # Enable debugging
        $<$<COMPILE_LANGUAGE:Fortran>:-O2> # Optimize to level 2
        $<$<COMPILE_LANGUAGE:Fortran>:-m4> # Set message level to Error only

        $<$<COMPILE_LANGUAGE:HIP>:-g>      # Enable debugging
        $<$<COMPILE_LANGUAGE:HIP>:-O3>     # Optimize to level 3
        $<$<COMPILE_LANGUAGE:HIP>:--offload-arch=gfx90a> # Set the target architecture for HIP
    )
endif()

# ============================================================================ #
# Define the runtime checks for Cray compilers
# ============================================================================ #
#
# Runtime checks (Capitalized throw error on first occurrence):
# - a:    Enables argument present check on entry.
# - b:    Enables checking of array bounds.
# - c:    Enables conformance checking of array operands in array expressions.
# - p:    Generates run time code to check the association or allocation status
#         of referenced POINTER variables, ALLOCATABLE arrays, or
#         assumed-shape arrays.
# - s:    Enables checking of character substring bounds.

if (NEKO_TOP_RUNTIME_CHECKS STREQUAL "ON")
    add_compile_options(
        $<$<COMPILE_LANGUAGE:Fortran>:-ei>      # Initialised values to NaN
        $<$<COMPILE_LANGUAGE:Fortran>:-RABCDPS> # Run time checks
    )
elseif (NEKO_TOP_RUNTIME_CHECKS STREQUAL "AUTO"
        AND CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(
        $<$<COMPILE_LANGUAGE:Fortran>:-ei>      # Initialised values to NaN
        $<$<COMPILE_LANGUAGE:Fortran>:-RAbCDpS> # Run time checks (Bounds warning)
    )
elseif (NOT NEKO_TOP_RUNTIME_CHECKS STREQUAL "OFF")
    message(FATAL_ERROR
        "Unknown NEKO_TOP_RUNTIME_CHECKS: '${NEKO_TOP_RUNTIME_CHECKS}'. "
        "Valid options are AUTO, ON or OFF.")
endif()
