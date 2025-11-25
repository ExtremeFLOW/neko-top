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
        message(STATUS "hipsolver NOT found in legacy locations → trying hip/hip-lang package method")

        # --- DEBUG: Find ALL ROCm installations and search EVERYWHERE ---
        message(STATUS "=== DEBUG: Searching for ROCm installations ===")
        
        # Find all possible ROCm directories
        file(GLOB ROCM_DIRS "/opt/rocm*")
        foreach(rocm_dir ${ROCM_DIRS})
            if(IS_DIRECTORY ${rocm_dir})
                message(STATUS "Found ROCm directory: ${rocm_dir}")
            endif()
        endforeach()

        # --- Try to locate HIP installation path ---
        set(HIP_ROOT "")
        message(STATUS "checking 1")
        find_package(hip-lang CONFIG QUIET)
        if(hip-lang_FOUND)
            message(STATUS "Found hip-lang package at: ${hip-lang_DIR}")
            get_filename_component(HIP_ROOT "${hip-lang_DIR}/../../.." ABSOLUTE)
            message(STATUS "Extracted HIP_ROOT: ${HIP_ROOT}")
        endif()
        message(STATUS "checking 2")

        if (NOT HIP_ROOT)
            # Fallback: use the most recent ROCm installation
            if(EXISTS "/opt/rocm-7.1.0")
                set(HIP_ROOT "/opt/rocm-7.1.0")
            elseif(EXISTS "/opt/rocm")
                set(HIP_ROOT "/opt/rocm")
            endif()
            message(STATUS "Using fallback HIP_ROOT: ${HIP_ROOT}")
        endif()

        if (HIP_ROOT)
            message(STATUS "Detected HIP root: ${HIP_ROOT}")
        else()
            message(FATAL_ERROR "Could not locate HIP installation via hip or hip-lang packages.")
        endif()

        # --- DEBUG: Search EVERYWHERE in HIP_ROOT for hipsolver ---
        message(STATUS "=== DEBUG: Searching for hipsolver in ${HIP_ROOT} ===")
        
        # Recursive search for any file containing "hipsolver"
        execute_process(
            COMMAND find "${HIP_ROOT}" -name "*hipsolver*" -type f
            OUTPUT_VARIABLE HIPSOLVER_FILES
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        
        if(HIPSOLVER_FILES)
            message(STATUS "Found hipsolver-related files:")
            string(REPLACE "\n" "\n  - " HIPSOLVER_FILES_FORMATTED "  - ${HIPSOLVER_FILES}")
            message(STATUS "${HIPSOLVER_FILES_FORMATTED}")
        else()
            message(STATUS "No hipsolver-related files found in ${HIP_ROOT}")
        endif()

        # Search for libhipsolver.so specifically in ALL subdirectories
        execute_process(
            COMMAND find "${HIP_ROOT}" -name "libhipsolver.so" -type f
            OUTPUT_VARIABLE HIPSOLVER_LIB_PATH
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )

        if(HIPSOLVER_LIB_PATH)
            message(STATUS "Found libhipsolver.so at: ${HIPSOLVER_LIB_PATH}")
            set(HIPSOLVER_LIB "${HIPSOLVER_LIB_PATH}")
        else()
            # Try alternative library names and locations
            message(STATUS "=== DEBUG: Trying alternative search methods ===")
            
            # Method 1: Search in common library directories
            set(SEARCH_DIRS
                "${HIP_ROOT}/lib"
                "${HIP_ROOT}/lib64"
                "${HIP_ROOT}/lib/x86_64-linux-gnu"
                "${HIP_ROOT}/hip/lib"
                "${HIP_ROOT}/hip/lib64"
                "${HIP_ROOT}/hipsolver/lib"
                "${HIP_ROOT}/hipsolver/lib64"
                "${HIP_ROOT}/rocm/lib"
                "${HIP_ROOT}/rocm/lib64"
            )
            
            foreach(search_dir ${SEARCH_DIRS})
                if(EXISTS "${search_dir}")
                    message(STATUS "Searching in: ${search_dir}")
                    file(GLOB LIB_FILES "${search_dir}/libhipsolver*")
                    foreach(lib_file ${LIB_FILES})
                        message(STATUS "  Found: ${lib_file}")
                    endforeach()
                endif()
            endforeach()

            # Method 2: Try find_library with all possible paths
            find_library(
                HIPSOLVER_LIB
                NAMES hipsolver
                PATHS ${SEARCH_DIRS}
                NO_DEFAULT_PATH
            )

            if(NOT HIPSOLVER_LIB)
                # Method 3: Last resort - search system wide
                message(STATUS "=== DEBUG: Trying system-wide search ===")
                find_library(HIPSOLVER_LIB hipsolver)
            endif()
        endif()

        if (HIPSOLVER_LIB)
            message(STATUS "Found hipSOLVER: ${HIPSOLVER_LIB}")
        else()
            message(FATAL_ERROR "Could not find libhipsolver.so anywhere! Searched in ${HIP_ROOT} and system-wide")
        endif()

        # --- Find other libraries with similar exhaustive search ---
        message(STATUS "=== DEBUG: Searching for hipblas and amdhip64 ===")
        
        # Search for hipblas
        execute_process(
            COMMAND find "${HIP_ROOT}" -name "libhipblas.so" -type f
            OUTPUT_VARIABLE HIPBLAS_LIB_PATH
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        
        if(HIPBLAS_LIB_PATH)
            set(HIPBLAS_LIB "${HIPBLAS_LIB_PATH}")
        else()
            find_library(HIPBLAS_LIB hipblas
                PATHS ${SEARCH_DIRS}
                NO_DEFAULT_PATH
            )
            if(NOT HIPBLAS_LIB)
                find_library(HIPBLAS_LIB hipblas)  # System-wide
            endif()
        endif()

        # Search for amdhip64
        execute_process(
            COMMAND find "${HIP_ROOT}" -name "libamdhip64.so" -type f
            OUTPUT_VARIABLE AMDHIP64_LIB_PATH
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        
        if(AMDHIP64_LIB_PATH)
            set(AMDHIP64_LIB "${AMDHIP64_LIB_PATH}")
        else()
            find_library(AMDHIP64_LIB amdhip64
                PATHS ${SEARCH_DIRS}
                NO_DEFAULT_PATH
            )
            if(NOT AMDHIP64_LIB)
                find_library(AMDHIP64_LIB amdhip64)  # System-wide
            endif()
        endif()

        if (NOT HIPBLAS_LIB)
            message(FATAL_ERROR "Could not find hipblas")
        endif()
        if (NOT AMDHIP64_LIB)
            message(FATAL_ERROR "Could not find amdhip64")
        endif()

        message(STATUS "hipSOLVER: ${HIPSOLVER_LIB}")
        message(STATUS "hipBLAS:   ${HIPBLAS_LIB}")
        message(STATUS "amdhip64:  ${AMDHIP64_LIB}")

        # --- FINAL LINK ---
        target_link_libraries(neko-top
            ${HIPSOLVER_LIB}
            ${HIPBLAS_LIB}
            ${AMDHIP64_LIB}
        )

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
