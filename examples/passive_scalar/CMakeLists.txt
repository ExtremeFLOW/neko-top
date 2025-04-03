# Simple setup to compile the current example.
set(DRIVER_TYPE "custom")
set(DRIVER_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/driver.f90)
set(EXTRA_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/passive_scalar_user.f90
)

build_example(${DRIVER_TYPE})
