# Extend the Neko-TOP target with the mapping functions
# as well as their "adjoint" (chain rule)
#
# Currently we only have
# - linear
# - RAMP
# - Borrvall & Petersson
# - PDE filter
# - SIMP
# - Heaviside mapping
#
# Possible extensions include
# - exponential

target_sources(neko-top
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/mapping.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/mapping_fctry.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/mapping_handler.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/PDE_filter_mapping.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/linear_mapping.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/RAMP_mapping.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/device_RAMP_mapping.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/Borrvall_Petersson_mapping.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/device_Borrvall_Petersson_mapping.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/SIMP_mapping.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/device_SIMP_mapping.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/cpu/RAMP_mapping_cpu.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/cpu/Borrvall_Petersson_mapping_cpu.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/cpu/SIMP_mapping_cpu.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/heaviside_mapping.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/device_heaviside_mapping.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/cpu/heaviside_mapping_cpu.f90
)

if(DEVICE_TYPE STREQUAL "CUDA")
    target_sources(neko-top
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/device/cuda/RAMP_mapping.cu
            ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/device/cuda/Borrvall_Petersson_mapping.cu
            ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/device/cuda/SIMP_mapping.cu
            ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/device/cuda/heaviside_mapping.cu
    )

    target_include_directories(neko-top
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/device/cuda
    )

elseif(DEVICE_TYPE STREQUAL "HIP")
    target_sources(neko-top
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/device/hip/RAMP_mapping.hip
            ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/device/hip/Borrvall_Petersson_mapping.hip
            ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/device/hip/SIMP_mapping.hip
            ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/device/hip/heaviside_mapping.hip
    )

    target_include_directories(neko-top
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/device/hip
    )
endif()
