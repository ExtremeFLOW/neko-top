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
        ${CMAKE_CURRENT_SOURCE_DIR}/problems/objectives/augmented_lagrangian_objective.f90
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
    $<$<BOOL:${CUDAToolkit_FOUND}>:CUDA::cusolver>
    $<$<BOOL:${CUDAToolkit_FOUND}>:CUDA::cudart>
)

if(${DEVICE_TYPE} STREQUAL "HIP")
    list(APPEND CMAKE_PREFIX_PATH
        /opt/rocm
        /opt/rocm-6.3.3
    )

    # --- FIRST TRY: Dardel/legacy ROCm ---
    find_library(HIPSOLVER_LIB hipsolver)

    if (HIPSOLVER_LIB)
        message(STATUS "Found legacy hipSOLVER: ${HIPSOLVER_LIB}")

        find_library(HIPBLAS_LIB hipblas REQUIRED)
        find_library(AMDHIP64_LIB amdhip64 REQUIRED)

        target_link_libraries(neko-top
            ${HIPSOLVER_LIB}
            ${HIPBLAS_LIB}
            ${AMDHIP64_LIB}
        )

    else()
        message(STATUS "hipsolver NOT found in legacy locations → checking rocm-libraries")

        # --- Check what ROCm packages are actually installed ---
        execute_process(
            COMMAND dpkg -l | grep -E "(rocm|hip)" | grep -E "(lib|dev)"
            OUTPUT_VARIABLE ROCM_PACKAGES
        )
        message(STATUS "Installed ROCm packages:\n${ROCM_PACKAGES}")

        # --- Try to find hipsolver via rocm-libraries location ---
        message(STATUS "=== Searching for hipsolver in rocm-libraries locations ===")
        
        # Common locations in ROCm 7.x with rocm-libraries
        set(ROCM_LIBS_PATHS
            /opt/rocm-7.1.0
            /opt/rocm
            /usr/lib/x86_64-linux-gnu
            /usr/local/cuda/lib64  # Sometimes hipSOLVER puts things here for compatibility
        )
        
        # Try different library naming patterns
        set(HIPSOLVER_NAMES hipsolver libhipsolver.so.1 libhipsolver.so)
        
        foreach(roc_path ${ROCM_LIBS_PATHS})
            foreach(lib_name ${HIPSOLVER_NAMES})
                if(EXISTS "${roc_path}/lib/${lib_name}" OR EXISTS "${roc_path}/lib64/${lib_name}")
                    if(EXISTS "${roc_path}/lib/${lib_name}")
                        set(HIPSOLVER_LIB "${roc_path}/lib/${lib_name}")
                    else()
                        set(HIPSOLVER_LIB "${roc_path}/lib64/${lib_name}")
                    endif()
                    message(STATUS "Found hipsolver: ${HIPSOLVER_LIB}")
                    break()
                endif()
            endforeach()
            if(HIPSOLVER_LIB)
                break()
            endif()
        endforeach()

        # --- If still not found, try to install the correct package ---
        if(NOT HIPSOLVER_LIB)
            message(STATUS "hipsolver not found - attempting to identify correct package...")
            
            # Check available ROCm packages
            execute_process(
                COMMAND apt-cache search "hipsolver"
                OUTPUT_VARIABLE AVAILABLE_HIPSOLVER_PACKAGES
            )
            message(STATUS "Available hipsolver packages:\n${AVAILABLE_HIPSOLVER_PACKAGES}")
            
            # Try installing hipsolver from rocm-libraries
            message(STATUS "Installing hipsolver via rocm-libraries package...")
            execute_process(
                COMMAND apt-get update && apt-get install -y libhipsolver-dev libhipsolver1
                RESULT_VARIABLE INSTALL_RESULT
                OUTPUT_VARIABLE INSTALL_OUTPUT
                ERROR_VARIABLE INSTALL_ERROR
            )
            
            if(INSTALL_RESULT EQUAL 0)
                message(STATUS "Successfully installed hipsolver")
                find_library(HIPSOLVER_LIB hipsolver)
            else()
                message(STATUS "Install failed, output: ${INSTALL_OUTPUT}")
                message(STATUS "Install error: ${INSTALL_ERROR}")
                
                # Try alternative package names used in rocm-libraries
                execute_process(
                    COMMAND apt-get install -y rocm-hipsolver
                    RESULT_VARIABLE INSTALL_RESULT2
                )
                if(INSTALL_RESULT2 EQUAL 0)
                    find_library(HIPSOLVER_LIB hipsolver)
                endif()
            endif()
        endif()

        # --- If hipsolver is truly unavailable, use rocsolver directly ---
        if(NOT HIPSOLVER_LIB)
            message(STATUS "hipsolver unavailable - using rocsolver directly")
            find_library(ROCSOLVER_LIB 
                NAMES rocsolver
                PATHS /opt/rocm-7.1.0/lib /opt/rocm-7.1.0/lib64
                NO_DEFAULT_PATH
            )
            
            if(ROCSOLVER_LIB)
                message(STATUS "Using rocsolver directly: ${ROCSOLVER_LIB}")
                set(HIPSOLVER_LIB ${ROCSOLVER_LIB})
            else()
                message(WARNING "Neither hipsolver nor rocsolver found - LAPACK functionality disabled")
                # Continue without solver libraries
            endif()
        endif()

        # --- Find other libraries ---
        find_library(HIPBLAS_LIB hipblas REQUIRED)
        find_library(AMDHIP64_LIB amdhip64 REQUIRED)

        message(STATUS "hipSOLVER: ${HIPSOLVER_LIB}")
        message(STATUS "hipBLAS:   ${HIPBLAS_LIB}")
        message(STATUS "amdhip64:  ${AMDHIP64_LIB}")

        # --- FINAL LINK ---
        target_link_libraries(neko-top
            ${HIPBLAS_LIB}
            ${AMDHIP64_LIB}
        )
        
        if(HIPSOLVER_LIB)
            target_link_libraries(neko-top ${HIPSOLVER_LIB})
        endif()

    endif()
endif()


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
