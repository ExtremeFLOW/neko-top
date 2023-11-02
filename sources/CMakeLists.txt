# ============================================================================ #
# CMakeLists.txt for the Neko-TOP library
# ============================================================================ #
# Set CMake Variables.
set( CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules )

# ============================================================================ #
# Define the library
file( GLOB test_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.f90 )

add_library( Neko-TOP STATIC ${test_SOURCES} )

# ============================================================================ #
# Find external libraries

find_package( PkgConfig REQUIRED )
pkg_check_modules( neko REQUIRED IMPORTED_TARGET neko )
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

find_package(MPI QUIET )
find_package(CUDAToolkit QUIET)

# ============================================================================ #
# Specify compile options

target_include_directories(Neko-TOP
    PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY}
)

target_link_libraries(Neko-TOP
    jsonfortran-gnu::jsonfortran LAPACK::LAPACK BLAS::BLAS PkgConfig::neko
)

if( MPI_FOUND )
    target_link_libraries(Neko-TOP MPI::MPI_Fortran)
endif()

if( CUDAToolkit_FOUND )
    target_link_libraries(Neko-TOP CUDA::cudart)
endif()

# ============================================================================ #
# Specify installation instructions

install( TARGETS ${PROJECT_NAME} DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/Lib/ )
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY} DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})