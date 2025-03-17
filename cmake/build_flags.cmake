# ============================================================================ #
# This file contains the compiler flags for the project.
#
# This allow us to set the flags for different compilers and different build
# types. The flags are set using the `add_compile_options` command.

# ============================================================================ #
# Set the compiler flags for GNU compilers
# ============================================================================ #
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        add_compile_options(
            # Debugging flags
            $<$<COMPILE_LANGUAGE:Fortran>:-g>
            $<$<COMPILE_LANGUAGE:Fortran>:-Og>

            # Warnings which for debug builds are errors
            $<$<COMPILE_LANGUAGE:Fortran>:-Wall>
            $<$<COMPILE_LANGUAGE:Fortran>:-Wextra>
            $<$<COMPILE_LANGUAGE:Fortran>:-Werror>
            $<$<COMPILE_LANGUAGE:Fortran>:-Wno-unused-dummy-argument>
            $<$<COMPILE_LANGUAGE:Fortran>:-Wno-function-elimination>
            $<$<COMPILE_LANGUAGE:Fortran>:-Wno-missing-include-dirs>
            $<$<COMPILE_LANGUAGE:Fortran>:-Wcompare-reals>
            $<$<COMPILE_LANGUAGE:Fortran>:-pedantic-errors>
        )
    elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
        add_compile_options(
            $<$<COMPILE_LANGUAGE:Fortran>:-O3>
            $<$<COMPILE_LANGUAGE:Fortran>:-funroll-loops>
            $<$<COMPILE_LANGUAGE:Fortran>:-flto>
            $<$<COMPILE_LANGUAGE:Fortran>:-fwhole-program>
            $<$<COMPILE_LANGUAGE:Fortran>:-Werror>
        )
    elseif(CMAKE_BUILD_TYPE STREQUAL "Testing")
        add_compile_options(
            # Debugging flags
            $<$<COMPILE_LANGUAGE:Fortran>:-g>
            $<$<COMPILE_LANGUAGE:Fortran>:-Og>

            # Warnings which for debug builds are errors
            $<$<COMPILE_LANGUAGE:Fortran>:-Wall>
            $<$<COMPILE_LANGUAGE:Fortran>:-Wextra>
            $<$<COMPILE_LANGUAGE:Fortran>:-Werror>
            $<$<COMPILE_LANGUAGE:Fortran>:-Wsurprising>
            $<$<COMPILE_LANGUAGE:Fortran>:-Warray-temporaries>
            $<$<COMPILE_LANGUAGE:Fortran>:-Wimplicit-interface>
            $<$<COMPILE_LANGUAGE:Fortran>:-Wno-unused-dummy-argument>
            $<$<COMPILE_LANGUAGE:Fortran>:-Wno-function-elimination>
            $<$<COMPILE_LANGUAGE:Fortran>:-Wno-missing-include-dirs>
            $<$<COMPILE_LANGUAGE:Fortran>:-Wcompare-reals>

            # Debugging checks
            $<$<COMPILE_LANGUAGE:Fortran>:-fcheck=all>
            $<$<COMPILE_LANGUAGE:Fortran>:-fbacktrace>
            $<$<COMPILE_LANGUAGE:Fortran>:-ffpe-trap=invalid,zero,overflow>
            $<$<COMPILE_LANGUAGE:Fortran>:-finit-real=snan>
            $<$<COMPILE_LANGUAGE:Fortran>:-finit-integer=-99999999>
            $<$<COMPILE_LANGUAGE:Fortran>:-finit-character=1>
            $<$<COMPILE_LANGUAGE:Fortran>:-finit-logical=true>
            $<$<COMPILE_LANGUAGE:Fortran>:-finit-local-zero>
        )
    endif()
endif()
