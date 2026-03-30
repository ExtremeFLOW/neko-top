# ============================================================================ #
# CMakeLists.txt for the Neko-TOP library
# ============================================================================ #
# Extend the neko-top library with the Problem module

target_sources(neko-top
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/problem.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/base_functional.f90

        ${CMAKE_CURRENT_SOURCE_DIR}/objective.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/constraint.f90
)

# Add subdirectories for objectives and constraints
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/objectives)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/constraints)
