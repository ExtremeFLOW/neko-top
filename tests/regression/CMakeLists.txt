# This file is used to define the regression tests for the Neko-TOP library.

# Add the regression tests
add_subdirectory(rugby_ball/)
add_dependencies(Tests Tests-regression-rugby_ball)
