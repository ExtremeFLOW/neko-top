# ============================================================================ #
# CMakeLists.txt for the Neko-TOP library
# ============================================================================ #
# HDF5 support for Neko-TOP
# All files in this directory are only compiled if HDF5 is found. Each of them
# should be a submodule with implementations of HDF5 specific IO routines for
# different Neko-TOP objects.

if (NOT hdf5_FOUND)
    return()
endif()

target_sources(neko-top
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/design_hdf5_io.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/mma_hdf5_io.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/optimizer_hdf5_io.f90
)
