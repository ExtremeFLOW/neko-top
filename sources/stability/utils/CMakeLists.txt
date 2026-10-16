target_sources(
    neko-top
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/neko_stability_propagator.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/floquet_output_utils.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/stability_case_utils.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/stability_solver_tolerances.f90
)
