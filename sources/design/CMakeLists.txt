# ============================================================================ #
# CMakeLists.txt for the Neko-TOP library
# ============================================================================ #
# Extend the neko-top library with the Design module

target_sources(neko-top
    PRIVATE
        # Main modules for design objects
        ${CMAKE_CURRENT_SOURCE_DIR}/design.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/design_factory.f90

        # Utility modules for the design objects
        ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/design_io_hdf5.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/bcknd/set_optimization_ic.f90

        # Specific design types
        ${CMAKE_CURRENT_SOURCE_DIR}/design_types/design_brinkman.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/design_types/design_simple.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/design_types/design_simplefield.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/design_types/design_3dto1d.f90
)
