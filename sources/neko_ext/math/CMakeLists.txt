# Extend the Neko-TOP target with the math extension module
#
# This module provides additional math functions for the Neko-TOP target.
#
# The math extension module is implemented in Fortran and CUDA. The CUDA
# implementation is only compiled if the CUDAToolkit is available.

target_sources(neko-top
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/math_ext.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/device_math_ext.f90
)

target_include_directories(neko-top
    PRIVATE
        ${NEKO_DIR}/src
)

if(CMAKE_CUDA_COMPILER)
    target_sources(neko-top
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/device/cuda/math_ext.cu
            ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/device/cuda/math_ext_kernel.h
    )

    target_include_directories(neko-top
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/device/cuda
    )
    if(CMAKE_VERSION VERSION_GREATER "3.24.0")
        set_target_properties(neko-top PROPERTIES CUDA_ARCHITECTURES native)
    endif()
endif()
