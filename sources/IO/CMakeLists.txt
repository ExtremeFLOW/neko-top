# ============================================================================ #
# CMakeLists.txt for the Neko-TOP library
# ============================================================================ #
# Extend the neko-top library with the MMA module

target_sources(neko-top
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/hdf5/design_hdf5_checkpoint.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/hdf5/mma_hdf5_checkpoint.f90
)
