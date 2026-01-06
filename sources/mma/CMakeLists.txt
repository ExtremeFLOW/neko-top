# ============================================================================ #
# CMakeLists.txt for the Neko-TOP library
# ============================================================================ #
# Extend the neko-top library with the MMA module
target_sources(neko-top
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/mma.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/mma_io_hdf5.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/cpu/mma_cpu.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/device/mma_device.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/device/device_mma_math.f90

        ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/device/cuda_mma_math.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/device/cuda/mma.cu

        ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/device/hip_mma_math.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/device/hip/mma.hip
)
