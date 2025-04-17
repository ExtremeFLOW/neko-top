# Extend the Neko-TOP target with the mapping functions
# as well as their "adjoint" (chain rule)
#
# Currently we only have
# - linear
# - RAMP (both convexities)
# - PDE filter
#
# Possible extensions include
# - exponential
# - SIMP
# - sigmoid/heavyside step

target_sources(neko-top
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/mapping.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/mapping_fctry.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/mapping_handler.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/PDE_filter_mapping.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/linear_mapping.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/RAMP_mapping.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/device_RAMP_mapping.f90
)


if(DEVICE_TYPE STREQUAL "CUDA")
    target_sources(neko-top
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/device/cuda/RAMP_mapping.cu
            ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/device/cuda/RAMP_mapping_kernel.h
    )

    target_include_directories(neko-top
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/device/cuda
    )

elseif(DEVICE_TYPE STREQUAL "HIP")
    # HIP kernel compiled in separate library to ensure hipcc is used
    add_library(hip_kernels STATIC
        ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/device/hip/RAMP_mapping.hip
    )

    set_target_properties(hip_kernels PROPERTIES
        LINKER_LANGUAGE HIP
    )

    target_include_directories(hip_kernels
        PUBLIC
            ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/device/hip
    )
    target_include_directories(hip_kernels
        PRIVATE
            ${NEKO_DIR}/src
    )

    target_link_libraries(neko-top hip_kernels)

    # Add the kernel headers for inclusion (but NOT as sources!)
    target_include_directories(neko-top
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/device/hip
    )
endif()

