# ============================================================================ #
# Define the top target for unit tests.

add_custom_target(Unit
    COMMENT "Build all unit tests."
)
set_target_properties(Unit PROPERTIES
    FOLDER "Tests/Unit"
)
add_dependencies(Tests Unit)

# Add each unit test subdirectory
add_subdirectory(checkpointing)