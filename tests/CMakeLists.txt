# This file will define the tests for the Neko-TOP library.
# Here we define the custom targets for the tests. Which allow us to build all
# the tests by targetting the specific target.

# Add the custom target for the tests
add_custom_target(Tests DEPENDS neko-top)

# Add the unit tests
add_subdirectory(unit)

# Add the regression tests
add_subdirectory(regression)
