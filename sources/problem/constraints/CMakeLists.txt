# ============================================================================ #
# CMakeLists.txt for the Neko-TOP library
# ============================================================================ #
# Extend the neko-top library with the Problem module

target_sources(neko-top
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/constraint_factory.f90

        # List of specific Constraints
        ${CMAKE_CURRENT_SOURCE_DIR}/volume_constraint.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/dummy_constraint.f90
)