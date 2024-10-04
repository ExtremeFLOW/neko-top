# Simple setup to compile the current example.
set(DRIVER_TYPE "custom")
set(DRIVER_SOURCE driver.f90)
set(EXTRA_SOURCES
    user.f90
    simulation_adjoint.f90
)

build_example(${DRIVER_TYPE})
