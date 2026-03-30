# ============================================================================ #
# CMakeLists.txt for the Neko-TOP library
# ============================================================================ #
# Extend the neko-top library with the Problem module

target_sources(neko-top
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/objective_factory.f90

        # List of specific objectives
        ${CMAKE_CURRENT_SOURCE_DIR}/brinkman_dissipation_objective.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/viscous_dissipation_objective.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/scalar_mixing_objective_function.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/augmented_lagrangian_objective.f90
)