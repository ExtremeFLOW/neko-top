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
    find_package(CUDAToolkit REQUIRED)

    if(DEFINED ENV{CUDA_ARCH})
        set(CMAKE_CUDA_ARCHITECTURES "$ENV{CUDA_ARCH}")
    endif()

    message(STATUS "CUDA architecture(s) set to: ${CMAKE_CUDA_ARCHITECTURES}")

elseif(${DEVICE_TYPE} STREQUAL "HIP")
    enable_language(HIP)
    enable_language(C)
    add_compile_definitions(HAVE_HIP=1)
    add_compile_options($<$<COMPILE_LANGUAGE:HIP>:-fPIE>)

    # Keep HIP compiler flags consistent with external Neko when provided.
    set(NEKO_TOP_HIP_FLAGS "")
    if(DEFINED ENV{NEKO_HIPCC_FLAGS} AND
       NOT "$ENV{NEKO_HIPCC_FLAGS}" STREQUAL "")
        set(NEKO_TOP_HIP_FLAGS "$ENV{NEKO_HIPCC_FLAGS}")
    elseif(DEFINED ENV{HIP_HIPCC_FLAGS} AND
           NOT "$ENV{HIP_HIPCC_FLAGS}" STREQUAL "")
        set(NEKO_TOP_HIP_FLAGS "$ENV{HIP_HIPCC_FLAGS}")
    endif()

    if(NOT "${NEKO_TOP_HIP_FLAGS}" STREQUAL "")
        separate_arguments(NEKO_TOP_HIP_FLAGS_LIST NATIVE_COMMAND
                           "${NEKO_TOP_HIP_FLAGS}")
        foreach(flag IN LISTS NEKO_TOP_HIP_FLAGS_LIST)
            add_compile_options($<$<COMPILE_LANGUAGE:HIP>:${flag}>)
        endforeach()
    endif()

    find_package(hipblas REQUIRED)
    find_package(hipsolver REQUIRED)

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
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/problem)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/mma)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/design)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/IO)

target_sources(neko-top
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/neko_top.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/user_modules/initial_conditions.f90

        ${CMAKE_CURRENT_SOURCE_DIR}/simulation/simulation.f90

        ${CMAKE_CURRENT_SOURCE_DIR}/optimizer/optimizer.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/optimizer/optimizer_factory.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/optimizer/mma_optimizer.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/continuation/continuation_scheduler.f90
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

target_link_libraries(neko-top PUBLIC
    PkgConfig::neko
    PkgConfig::json-fortran
    MPI::MPI_Fortran
    OpenMP::OpenMP_Fortran
    $<$<BOOL:${BLAS_FOUND}>:BLAS::BLAS>
    $<$<BOOL:${LAPACK_FOUND}>:LAPACK::LAPACK>
    $<$<BOOL:${CUDAToolkit_FOUND}>:CUDA::cusolver>
    $<$<BOOL:${CUDAToolkit_FOUND}>:CUDA::cudart>
    $<$<BOOL:${hipblas_FOUND}>:roc::hipblas>
    $<$<BOOL:${hipsolver_FOUND}>:roc::hipsolver>
    neko_cxx_support  # ADIOS2/MPI C++/stdc++ must come after libneko
    $<$<BOOL:${HDF5_PKG_FOUND}>:PkgConfig::HDF5>
    $<$<BOOL:${HDF5_FOUND}>:HDF5::HDF5>
)

# ============================================================================ #
# Specify installation instructions

add_library(${PROJECT_NAME}::neko-top ALIAS neko-top)

# Install the neko-top library and its modules.
# These should adhere to the CMake installation conventions.

install(TARGETS neko_cxx_support
    EXPORT neko-top-targets
)
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
