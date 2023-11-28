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
# Determine if we wish to build the examples.
option(BUILD_EXAMPLES "Build the examples" ON)

if(NOT BUILD_EXAMPLES)
    return()
else()
    message(STATUS "Building examples")
endif()

# ============================================================================ #
# Include the build_example function
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
include(build_example)

# Set the examples directory
set(EXAMPLES_DIR ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "Examples directory")

# ============================================================================ #
# Add the example directories.

# Find all subdirectories in the examples directory
file(GLOB EXAMPLES_DIRS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*)

# Add each subdirectory, if it contains a CMakelists.txt file
foreach(EXAMPLE_DIR ${EXAMPLES_DIRS})
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${EXAMPLE_DIR}/CMakeLists.txt)
        add_subdirectory(${EXAMPLE_DIR})
    endif()
endforeach()
