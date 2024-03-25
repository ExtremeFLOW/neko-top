# Simple setup to compile the current example.
set(DRIVER_TYPE "user")

set(EXTRA_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/heat_source.f90")

build_example(${DRIVER_TYPE})
