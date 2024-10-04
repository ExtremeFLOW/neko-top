# ============================================================================ #
# CMakeLists.txt for the Neko-TOP library
# ============================================================================ #
# Set CMake Variables.

# Setup the global inclusion for fortran
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})

add_library(neko-top STATIC)

# Add the local neko extension library
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/neko_ext)

# ============================================================================ #
# Define the library sources and targets
target_sources(neko-top
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/neko_top.f90

        ${CMAKE_CURRENT_SOURCE_DIR}/user_modules/initial_conditions.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/topology_optimization/design.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/topology_optimization/user_module.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/topology_optimization/sensitivity.f90

        ${CMAKE_CURRENT_SOURCE_DIR}/adjoint/adjoint_case.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/adjoint/adjoint_output.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/adjoint/adjoint_scheme.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/adjoint/adjoint_fctry.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/adjoint/adjoint_pnpn.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/adjoint/adjoint_ic.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/adjoint/adjoint_source_term.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/adjoint/adv_adjoint_no_dealias.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/adjoint/adv_adjoint_dealias.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/adjoint/advection_adjoint.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/adjoint/advection_adjoint_fctry.f90

        ${CMAKE_CURRENT_SOURCE_DIR}/simulation_components/steady_simcomp.f90

        ${CMAKE_CURRENT_SOURCE_DIR}/mma/mma.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/mma/bcknd/mma_vector.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/mma/bcknd/cpu/mma_cpu.f90
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
    DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/lib/
)
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
    DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}
)
