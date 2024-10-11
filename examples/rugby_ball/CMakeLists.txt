# Simple setup to compile the current example.
set(DRIVER_TYPE "custom")
set(DRIVER_SOURCE driver.f90)
set(EXTRA_SOURCES
    user.f90
    simulation_adjoint.f90
    source_terms/simple_brinkman_source_term.f90
    source_terms/adjoint_lube_source_term.f90
    source_terms/adjoint_minimum_dissipation_source_term.f90
    objectives/objective_function.f90
    objectives/minimum_dissipation_objective_function.f90
    objectives/volume_constraint.f90
    designs/topopt_design.f90
    problems/problem.f90
    problems/steady_state_problem.f90
    mapping_functions/mapping.f90
    mapping_functions/PDE_filter_mapping.f90
    mapping_functions/linear_mapping.f90
    mapping_functions/RAMP_mapping.f90
    mapping_functions/Casper_mapping.f90
)

build_example(${DRIVER_TYPE})
