# Simple setup to compile the current example.
# set(DRIVER_TYPE "custom")
# set(EXTRA_SOURCES
#     ${CMAKE_CURRENT_SOURCE_DIR}/static_mixer.f90
# )

set(DRIVER_TYPE "default")
set(EXTRA_SOURCES )

build_example(${DRIVER_TYPE})
