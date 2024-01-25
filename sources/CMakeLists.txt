# ============================================================================ #
# CMakeLists.txt for the Neko-TOP library
# ============================================================================ #
# Set CMake Variables.

# Setup the global inclusion for fortran
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})

add_library(Neko-TOP STATIC)

# Add the local neko extension library
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/neko_ext)

# ============================================================================ #
# Define the library sources and targets
target_sources(Neko-TOP
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/neko_top.f90

        ${CMAKE_CURRENT_SOURCE_DIR}/user_modules/initial_conditions.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/topology_optimization/design.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/topology_optimization/user_module.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/topology_optimization/sensitivity.f90
)

# ============================================================================ #
# Specify compile options

target_include_directories(Neko-TOP
    PUBLIC
        ${CMAKE_Fortran_MODULE_DIRECTORY}
)

target_link_libraries(Neko-TOP
    BLAS::BLAS
    LAPACK::LAPACK
    jsonfortran-gnu::jsonfortran
    PkgConfig::neko
)

if(MPI_FOUND)
    target_link_libraries(Neko-TOP MPI::MPI_Fortran)
endif()

if(CUDAToolkit_FOUND)
    target_link_libraries(Neko-TOP CUDA::cudart)
    set(ENV{HAVE_CUDA} TRUE)
endif()

# ============================================================================ #
# Specify installation instructions

install(TARGETS ${PROJECT_NAME} DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/Lib/ )
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY} DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})
