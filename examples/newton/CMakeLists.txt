# Driver and source definitions for build_example()
set(DRIVER_TYPE "custom")
set(DRIVER ${CMAKE_CURRENT_SOURCE_DIR}/driver.f90)
build_example(${DRIVER_TYPE})