# ============================================================================ #
# Define the top target for regression tests.

add_custom_target(Regression
    COMMENT "Build all regression tests."
    DEPENDS neko-top
)
set_target_properties(Regression PROPERTIES FOLDER "Regression Tests")

if (BUILD_TESTING)
    set_target_properties(Regression PROPERTIES EXCLUDE_FROM_ALL FALSE)
else()
    set_target_properties(Regression PROPERTIES EXCLUDE_FROM_ALL TRUE)
endif()
