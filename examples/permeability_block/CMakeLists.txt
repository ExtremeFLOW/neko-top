# Simple setup to compile the current example.
set(DRIVER_TYPE "user")
set(EXTRA_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/permeability_user.f90
)

build_example(${DRIVER_TYPE})
