# Extend the Neko-TOP target with the mapping functions
# as well as their "adjoint" (chain rule)
#
# Currently we only have
# - linear
# - RAMP (both convexities)
# - PDE filter
# - SIMP
#
# Possible extensions include
# - exponential
# - sigmoid/heavyside step

target_sources(neko-top
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/mapping.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/mapping_fctry.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/mapping_handler.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/PDE_filter_mapping.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/linear_mapping.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/RAMP_mapping.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/SIMP_mapping.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/device_RAMP_mapping.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/device_SIMP_mapping.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/cpu/SIMP_mapping_cpu.f90
)

if(DEVICE_TYPE STREQUAL "CUDA")
    target_sources(neko-top
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/device/cuda/RAMP_mapping.cu
            ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/device/cuda/SIMP_mapping.cu
    )

    target_include_directories(neko-top
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/device/cuda
    )

elseif(DEVICE_TYPE STREQUAL "HIP")
    target_sources(neko-top
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/device/hip/RAMP_mapping.hip
            ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/device/hip/SIMP_mapping.hip
    )

    target_include_directories(neko-top
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/device/hip
    )
endif()
