# Simple setup to compile the current example.
set(DRIVER_TYPE "custom")
set(DRIVER_SOURCE driver.f90)
set(EXTRA_SOURCES 
    adjoint_mixing_scalar_source_term.f90 
    scalar_mixing_objective_function.f90
    user.f90 
    )

build_example(${DRIVER_TYPE})
