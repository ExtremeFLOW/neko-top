# Simple setup to compile the current example.
set(DRIVER_TYPE "custom")
set(DRIVER run.sh)

# These extra sources are picked up by build_example()
set(EXTRA_SOURCES
    insitu_cylinder_POD.f90
    usr_driver.f90
)

# This macro (from cmake/build_example.cmake) creates the "turb_pipe" target
# and sets OUTPUT_NAME to "neko" for you.
build_example(${DRIVER_TYPE})
