# Simple setup to compile the current example.
set(EXAMPLE_NAME "neko_examples_tgv")
set(EXAMPLE_TYPE "user")

# ============================================================================ #
# Setup of the external libraries.

# Ensure that the compiler can find Neko.
find_package(PkgConfig REQUIRED)
pkg_check_modules(neko REQUIRED IMPORTED_TARGET neko)

# Check for the required libraries.
find_package(jsonfortran-gnu REQUIRED)
find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)

# Check for the optional libraries.
find_package(CUDAToolkit QUIET)
find_package(MPI QUIET)

# ============================================================================ #
# Define the executable.

if (${EXAMPLE_TYPE} STREQUAL "user")
    set(SOURCES ${EXAMPLES_DIR}/usr_driver.f90 ${NEKO_DIR}/examples/tgv/tgv.f90)
else()
    set(SOURCES ${EXAMPLES_DIR}/driver.f90)
endif()

add_executable(${EXAMPLE_NAME} ${SOURCES})

# Set the output directory of the executable.
set_target_properties(${EXAMPLE_NAME}
    PROPERTIES
    OUTPUT_NAME "neko"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
)

# ============================================================================ #
# Link the executable to the required libraries.
target_link_libraries(${EXAMPLE_NAME}
    jsonfortran-gnu::jsonfortran LAPACK::LAPACK BLAS::BLAS PkgConfig::neko)

# If CUDA is available, link the executable to the CUDA runtime library.
if(CUDAToolkit_FOUND)
    target_link_libraries(${EXAMPLE_NAME} CUDA::cudart)
endif()

# If MPI is available, link the executable to the MPI library.
if(MPI_FOUND)
    target_link_libraries(${EXAMPLE_NAME} MPI::MPI_Fortran)
endif()

