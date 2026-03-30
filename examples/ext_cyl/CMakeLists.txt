# ==============================================================================
# LightKrylov Example CMake File
# ==============================================================================

# Driver and source definitions for build_example()
set(DRIVER_TYPE "custom")
set(DRIVER ${CMAKE_CURRENT_SOURCE_DIR}/driver.f90)
set(EXTRA_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/cylinder.f90
)

build_example(${DRIVER_TYPE})
