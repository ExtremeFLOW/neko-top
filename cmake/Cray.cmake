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
        $<$<COMPILE_LANGUAGE:Fortran>:-g>  # Enable debugging
        $<$<COMPILE_LANGUAGE:Fortran>:-m1> # Set message level to verbose
        $<$<COMPILE_LANGUAGE:Fortran>:-O0> # Disable optimization

        $<$<COMPILE_LANGUAGE:Fortran>:-RABCDS> # Run time checks
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
endif()
