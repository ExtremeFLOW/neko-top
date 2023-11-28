# Simple setup to compile the current example.
set(DRIVER_TYPE "user")
set(EXTRA_SOURCES
    ${CMAKE_SOURCE_DIR}/sources/user_modules/permeability_user.f90
)

build_example(${DRIVER_TYPE})
