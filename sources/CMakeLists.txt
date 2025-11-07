# ============================================================================ #
# CMakeLists.txt for the Neko-TOP library
# ============================================================================ #
# Set CMake Variables.

# Setup the global inclusion for fortran
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules
    CACHE PATH "Fortran module directory")

# ............................................................................ #
# Enable devices

if (${DEVICE_TYPE} STREQUAL "CUDA")
    enable_language(CUDA)
    add_compile_definitions(HAVE_CUDA=1)
elseif(${DEVICE_TYPE} STREQUAL "HIP")
    enable_language(HIP)
    add_compile_definitions(HAVE_HIP=1)

    add_compile_options($<$<COMPILE_LANGUAGE:HIP>:-fPIE>)

elseif(NOT ${DEVICE_TYPE} STREQUAL "CPU")
    message(FATAL_ERROR "Unsupported device type: ${DEVICE_TYPE}")
endif()

# ============================================================================ #
# Define the library sources and targets

add_library(neko-top)

# Add our partial libraries to the build.
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/neko_ext)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/mapping_functions)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/adjoint)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/source_terms)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/simulation_components)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/checkpoint)

target_sources(neko-top
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/neko_top.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/user_modules/initial_conditions.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/topology_optimization/user_module.f90

        ${CMAKE_CURRENT_SOURCE_DIR}/designs/design.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/designs/design_3dto1d.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/designs/design_factory.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/designs/design_brinkman.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/designs/design_simple.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/designs/design_simplefield.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/designs/set_optimization_ic.f90

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
        ${CMAKE_CURRENT_SOURCE_DIR}/problems/objectives/scalar_mixing_objective_function.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/problems/constraints/volume_constraint.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/problems/constraints/dummy_constraint.f90

        ${CMAKE_CURRENT_SOURCE_DIR}/optimizer/optimizer.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/optimizer/optimizer_factory.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/optimizer/mma.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/optimizer/mma_optimizer.f90
        # ${CMAKE_CURRENT_SOURCE_DIR}/mma/bcknd/mma_vector.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/optimizer/optimizer.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/optimizer/bcknd/cpu/mma_cpu.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/optimizer/bcknd/device/mma_device.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/optimizer/bcknd/device/device_mma_math.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/optimizer/bcknd/device/cuda_mma_math.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/optimizer/bcknd/device/cuda/mma.cu
        ${CMAKE_CURRENT_SOURCE_DIR}/optimizer/bcknd/device/hip_mma_math.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/optimizer/bcknd/device/hip/mma.hip
)

# ============================================================================ #
# Specify compile options

target_include_directories(neko-top PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
    $<INSTALL_INTERFACE:include/neko-top>
)

# ............................................................................ #
# Link to the required libraries.

target_link_libraries(neko-top
    PkgConfig::neko
    PkgConfig::json-fortran
    MPI::MPI_Fortran
    $<$<BOOL:${BLAS_FOUND}>:BLAS::BLAS>
    $<$<BOOL:${LAPACK_FOUND}>:LAPACK::LAPACK>
)

# ============================================================================ #
# Specify installation instructions

add_library(${PROJECT_NAME}::neko-top ALIAS neko-top)

# Install the neko-top library and its modules.
# These should adhere to the CMake installation conventions.

install(TARGETS neko-top
    EXPORT neko-top-targets
    ARCHIVE DESTINATION lib/neko-top
    LIBRARY DESTINATION lib/neko-top
    RUNTIME DESTINATION bin/neko-top
    INCLUDES DESTINATION include/neko-top
)
install(EXPORT neko-top-targets
    NAMESPACE ${PROJECT_NAME}::
    DESTINATION lib/cmake/${PROJECT_NAME}
)
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/
    DESTINATION include/neko-top
    FILES_MATCHING PATTERN "*.mod"
)
