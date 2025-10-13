# ============================================================================ #
# Define the top target for regression tests.

add_custom_target(Regression
    COMMENT "Build all regression tests."
)
set_target_properties(Regression PROPERTIES
    FOLDER "Tests/Regression"
)
add_dependencies(Tests Regression)

# Add each regression test subdirectory
