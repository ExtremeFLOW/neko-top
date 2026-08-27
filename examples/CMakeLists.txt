# ============================================================================ #
# Build definitions for all examples.
#
# Each subdirectory should contain a CMakeLists.txt file that defines the type
# of driver to use and the extra source files to compile. See
# CMakeListsTemplate.txt for an example.
#
# The following driver types are currently supported:
# - "default" : Unmodified neko driver.
# - "user"    : Driver initializing the neko library and calling the user code.
#                This is equivalent to calling makeneko from the original neko.
# - "custom"  : Custom driver defined by a f90 file. Either set through the
#                DRIVER cmake variable in the example CMakeLists.txt file or
#                stored as a file named "driver.f90" in the example directory.
# - "topopt"  : Topology optimization driver.

# ============================================================================ #
# Include the build_example function
include(build_example)

# Set the examples directory
set(EXAMPLES_DIR ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "Examples directory")

# ============================================================================ #
# Define our custom target.

# Add a custom target for the examples
add_custom_target(Examples
    COMMENT "Build all the examples."
    DEPENDS neko-top
)
set_target_properties(Examples PROPERTIES FOLDER "Examples")

# Add the documentation target to ALL if we have "BUILD_EXAMPLES" enabled
if (BUILD_EXAMPLES)
    set_target_properties(Examples PROPERTIES EXCLUDE_FROM_ALL FALSE)
else()
    set_target_properties(Examples PROPERTIES EXCLUDE_FROM_ALL TRUE)
endif()

# ============================================================================ #
# Add the example directories.

add_subdirectory(${EXAMPLES_DIR}/easy-E)
add_subdirectory(${EXAMPLES_DIR}/mma)
add_subdirectory(${EXAMPLES_DIR}/passive_scalar)
add_subdirectory(${EXAMPLES_DIR}/rugby_ball)
add_subdirectory(${EXAMPLES_DIR}/1d_beam)
add_subdirectory(${EXAMPLES_DIR}/static_mixer)
add_subdirectory(${EXAMPLES_DIR}/unsteady_mixer)
add_subdirectory(${EXAMPLES_DIR}/unsteady_mixer_small)
add_subdirectory(${EXAMPLES_DIR}/unsteady_rugby_ball)
add_subdirectory(${EXAMPLES_DIR}/wiggles)

# If the benchmark_mixer folder exists, add it as well
if (EXISTS ${EXAMPLES_DIR}/benchmark_mixer/CMakeLists.txt)
    add_subdirectory(${EXAMPLES_DIR}/benchmark_mixer)
endif()

# Add the benchmark directory if it exists.
if (EXISTS ${EXAMPLES_DIR}/benchmark/CMakeLists.txt)
    add_subdirectory(${EXAMPLES_DIR}/benchmark)
endif()