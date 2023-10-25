# ============================================================================ #
# Build definitions for all examples.
#
# Each subdirectory should contain a CMakeLists.txt file that defines the type
# of driver to use and the extra source files to compile. See
# CMakeListsTemplate.txt for an example.
#
# The following driver types are currently supported:
#  - "default" : Unmodified neko driver.
#  - "user"    : Driver initializing the neko library and calling the user code.
#                This is equivalent to calling makeneko from the original neko.
#  - "custom"  : Custom driver defined by a f90 file. Either set through the
#                DRIVER cmake vatiable in the example CMakeLists.txt file or
#                stored as a file named "driver.f90" in the example directory.


# ============================================================================ #
# Setup of the external libraries.

# Check for the required libraries.
find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)

# Check for the optional libraries.
find_package(CUDAToolkit QUIET)
find_package(MPI QUIET)

# ============================================================================ #
# Include the build_example function
include(build_example.cmake)

# Set the examples directory
set(EXAMPLES_DIR ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "Examples directory")

# ============================================================================ #
# Add the example directories.

add_subdirectory(neko_examples/tgv)
add_subdirectory(neko_examples/turb_channel)
add_subdirectory(neko_examples/2d_cylinder)
add_subdirectory(neko_examples/rayleigh-benard)
