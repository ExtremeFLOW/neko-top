# Simple setup to compile the current example.
set(DRIVER_TYPE "custom")
set(DRIVER_SOURCE driver.f90)
set(EXTRA_SOURCES
    # adjoint.f90
    user.f90
    steady_simcomp.f90
    # topopt.f90
    # solve_adjoint.f90
)

build_example(${DRIVER_TYPE})
