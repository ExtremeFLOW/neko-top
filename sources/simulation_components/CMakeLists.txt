# CMakeLists.txt for simulation components in Neko-TOP
#
# Extend the neko-top library with the simulation component implementations.
# Each simulation component is implemented in its own Fortran source file, while
# the neko_top_simcomps.f90 file aggregates them and extends the Neko simulation
# component system to include these new simulation components.

target_sources(neko-top
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/neko_top_simcomps.f90

        # Individual simulation components
        ${CMAKE_CURRENT_SOURCE_DIR}/steady_simcomp.f90
)
