# Simple setup to compile the current example.
set(DRIVER_TYPE "custom")
set(DRIVER_SOURCE driver.f90)
set(EXTRA_SOURCES
    user.f90
    # topopt.f90
    simulation_adjoint.f90
    simple_brinkman_source_term.f90
    new_design.f90
    adjoint_minimum_dissipation_source_term.f90
    objective_function.f90
    minimum_dissipation_objective_function.f90
)

build_example(${DRIVER_TYPE})
