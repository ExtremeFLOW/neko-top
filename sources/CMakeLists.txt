# ============================================================================ #
# CMakeLists.txt for the Neko-TOP library
# ============================================================================ #
# Set CMake Variables.

# Setup the global inclusion for fortran
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules
    CACHE PATH "Fortran module directory")
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})

add_library(neko-top STATIC)

# Set the output directory of the library.
set_target_properties(neko-top
    PROPERTIES
    OUTPUT_NAME "neko-top"
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
    INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_SOURCE_DIR}/include"
)

# Add the local neko extension library
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/neko_ext)

# Add the mapping extension library
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/mapping_functions)

# Add the adjoint library
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/adjoint)
# ============================================================================ #
# Define the library sources and targets
target_sources(neko-top
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/neko_top.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/user_modules/initial_conditions.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/topology_optimization/design.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/topology_optimization/user_module.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/topology_optimization/sensitivity.f90

        ${CMAKE_CURRENT_SOURCE_DIR}/source_terms/simple_brinkman_source_term.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/source_terms/adjoint_lube_source_term.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/source_terms/adjoint_minimum_dissipation_source_term.f90

        ${CMAKE_CURRENT_SOURCE_DIR}/designs/design.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/designs/topopt_design.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/simulation/simulation.f90

        # Problem setup
        # Components for the objectives and constraints
        ${CMAKE_CURRENT_SOURCE_DIR}/problems/problem.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/problems/base_functional.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/problems/objective.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/problems/objective_factory.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/problems/constraint.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/problems/constraint_factory.f90

        # List of specific problems, objectives and constraints
        ${CMAKE_CURRENT_SOURCE_DIR}/problems/objectives/lube_term_objective.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/problems/objectives/minimum_dissipation_objective.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/problems/constraints/volume_constraint.f90

        ${CMAKE_CURRENT_SOURCE_DIR}/simulation_components/steady_simcomp.f90

        ${CMAKE_CURRENT_SOURCE_DIR}/optimizer/optimizer.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/optimizer/optimizer_factory.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/optimizer/mma.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/optimizer/mma_optimizer.f90
        # ${CMAKE_CURRENT_SOURCE_DIR}/mma/bcknd/mma_vector.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/optimizer/bcknd/cpu/mma_cpu.f90
)

# ============================================================================ #
# Specify compile options

target_include_directories(neko-top
    PUBLIC
        ${CMAKE_Fortran_MODULE_DIRECTORY}
)

target_link_libraries(neko-top
    BLAS::BLAS
    LAPACK::LAPACK
    jsonfortran-gnu::jsonfortran
    PkgConfig::neko
)

if(MPI_FOUND)
    target_link_libraries(neko-top MPI::MPI_Fortran)
endif()

# ============================================================================ #
# Specify installation instructions

install(TARGETS neko-top
    DESTINATION ${CMAKE_SOURCE_DIR}/lib/
)
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
    DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}
)
