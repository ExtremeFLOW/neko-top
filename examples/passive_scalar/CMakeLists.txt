# Simple setup to compile the current example.
set(DRIVER_TYPE "custom")
set(DRIVER_SOURCE driver.f90)
set(EXTRA_SOURCES 
    adjoint_mixing_scalar_source_term.f90 
    scalar_mixing_objective_function.f90
    user.f90 
    ../rugby_ball/simulation_adjoint.f90
    ../rugby_ball/source_terms/simple_brinkman_source_term.f90
    ../rugby_ball/source_terms/adjoint_lube_source_term.f90
    ../rugby_ball/source_terms/adjoint_minimum_dissipation_source_term.f90
    ../rugby_ball/objectives/objective_function.f90
    ../rugby_ball/objectives/minimum_dissipation_objective_function.f90
    ../rugby_ball/objectives/volume_constraint.f90
    ../rugby_ball/designs/topopt_design.f90
    ../rugby_ball/problems/problem.f90
    ../rugby_ball/problems/steady_state_problem.f90
    ../rugby_ball/mapping_functions/mapping.f90
    ../rugby_ball/mapping_functions/PDE_filter_mapping.f90
    ../rugby_ball/mapping_functions/linear_mapping.f90
    ../rugby_ball/mapping_functions/RAMP_mapping.f90
    ../rugby_ball/mapping_functions/Casper_mapping.f90
    )

build_example(${DRIVER_TYPE})
