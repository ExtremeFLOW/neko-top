# CMakeLists.txt for the tests directory.
# This file is used to define the tests for the Neko-TOP library.

add_custom_target(Tests
    COMMENT "Build all tests."
    DEPENDS neko-top
)
set_target_properties(Tests PROPERTIES
    EXCLUDE_FROM_ALL $<IF:$<BOOL:BUILD_TESTING>,FALSE,TRUE>
    FOLDER "Tests"
)

# ============================================================================ #
# Define the top target for unit tests.

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/unit)

# ============================================================================ #
# Define the top target for regression tests.

add_custom_target(Regression EXCLUDE_FROM_ALL
    COMMENT "Build all regression tests."
    DEPENDS Tests
)
set_target_properties(Regression PROPERTIES FOLDER "Tests/Regression")
