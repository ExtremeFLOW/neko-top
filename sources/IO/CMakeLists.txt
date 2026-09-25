# ============================================================================ #
# CMakeLists.txt for the Neko-TOP library
# ============================================================================ #
# HDF5 support for Neko-TOP
# All files in this directory are only compiled if HDF5 is found. Each of them
# should be a submodule with implementations of HDF5 specific IO routines for
# different Neko-TOP objects.

target_sources(neko-top
    PRIVATE
        $<$<BOOL:${HAVE_HDF5}>:${CMAKE_CURRENT_SOURCE_DIR}/hdf5/design_hdf5_io.f90>
        $<$<BOOL:${HAVE_HDF5}>:${CMAKE_CURRENT_SOURCE_DIR}/hdf5/mma_hdf5_io.f90>
        $<$<BOOL:${HAVE_HDF5}>:${CMAKE_CURRENT_SOURCE_DIR}/hdf5/optimizer_hdf5_io.f90>
)
