# Extend the Neko-TOP target with the optimizers
#
# Currently we have
# - MMA
# - Steepest descent with move limits

target_sources(neko-top
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/optimizer.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/optimizer_factory.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/mma_optimizer.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/steepest_descent_optimizer.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/steepest_descent_math.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/device_steepest_descent.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/cpu/steepest_descent_cpu.f90
)

if(DEVICE_TYPE STREQUAL "CUDA")
    target_sources(neko-top
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/device/cuda/steepest_descent.cu
    )

    target_include_directories(neko-top
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/device/cuda
    )

elseif(DEVICE_TYPE STREQUAL "HIP")
    target_sources(neko-top
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/device/hip/steepest_descent.hip
    )

    target_include_directories(neko-top
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/device/hip
    )
endif()
