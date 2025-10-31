# Simple setup to compile the current example.
set(DRIVER_TYPE "topopt-user")
set(EXTRA_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/mixer.f90)

build_example(${DRIVER_TYPE})
