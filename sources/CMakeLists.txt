# ============================================================================ #
# CMakeLists.txt for the Neko-TOP library
# ============================================================================ #
# Set CMake Variables.
set( CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules )

# ============================================================================ #
# Define the library sources and targets
set(neko-top_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/neko_top.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/assorted/develop.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/assorted/topology.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/user_modules/initial_conditions.f90
)

add_library(Neko-TOP STATIC ${neko-top_SOURCES})

# ============================================================================ #
# Specify compile options

target_include_directories(Neko-TOP
    PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY}
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
endif()

# ============================================================================ #
# Specify installation instructions

install(TARGETS ${PROJECT_NAME} DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/Lib/ )
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY} DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})
