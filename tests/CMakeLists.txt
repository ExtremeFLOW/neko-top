# CMakeLists.txt for the tests directory.
# This file is used to define the tests for the Neko-TOP library.

if (NOT BUILD_TESTING)
    return()
endif()

add_custom_target(Tests
    COMMENT "Build all tests."
    DEPENDS neko-top
)
set_target_properties(Tests PROPERTIES
    EXCLUDE_FROM_ALL $<IF:$<BOOL:BUILD_TESTING>,FALSE,TRUE>
    FOLDER "Tests"
)

# Set the shared test sources directory. This is where modules and other shared
# test files are located.
set(SHARED_TEST_SOURCES ${CMAKE_SOURCE_DIR}/tests/shared)

# ============================================================================ #
# Add subdirectories for unit and regression tests.

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/unit)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/regression)
