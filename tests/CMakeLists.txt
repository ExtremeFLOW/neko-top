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

# ============================================================================ #
# Lint the registered tests: every DEPENDS entry must name an existing test.
# CTest silently ignores any other entry, e.g. a build target. The lint is
# also labelled "unit" so that it runs in CI's `ctest -L unit` selection; it
# takes well under a second and needs neither MPI nor a build.

add_test(NAME lint_test_depends
    COMMAND ${CMAKE_COMMAND}
        -DCTEST_COMMAND=${CMAKE_CTEST_COMMAND}
        -DTEST_DIR=${CMAKE_BINARY_DIR}
        -P ${CMAKE_CURRENT_SOURCE_DIR}/check_test_depends.cmake
)
set_tests_properties(lint_test_depends PROPERTIES
    LABELS "lint;unit"
    TIMEOUT 60
)
