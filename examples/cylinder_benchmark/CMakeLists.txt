# Simple setup to compile the current example.
set(DRIVER_TYPE "user")

set(EXTRA_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/user_setup.f90")

build_example(${DRIVER_TYPE})
