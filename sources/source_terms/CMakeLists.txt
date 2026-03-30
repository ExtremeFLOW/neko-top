# CMakeLists.txt for source terms in Neko-TOP
#
# Extend the neko-top library with the source term implementations.
# Each source term is implemented in its own Fortran source file, while the
# neko_top_source_terms.f90 file aggregates them and extends the Neko source
# term system to include these new source terms.

target_sources(neko-top PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/neko_top_source_terms.f90

        # Individual source terms
        ${CMAKE_CURRENT_SOURCE_DIR}/simple_brinkman_source_term.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/adjoint_brinkman_dissipation_source_term.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/adjoint_viscous_dissipation_source_term.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/adjoint_mixing_scalar_source_term.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/adjoint_scalar_convection_source_term.f90
)
