# Extensions to the neko library.
#
# This is a place new functions can be added which are intended to be moved over
# to neko proper at some point. We should make sure the code structure is
# similar to the neko library so that it is easy to move over.
#
# The functions in this file are not guaranteed to be stable. They may be
# removed or changed at any time.

# Add the neko-ext library which collects all the extensions.
target_sources(neko-top
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/neko_top_source_terms.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/simple_brinkman_source_term.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/adjoint_lube_source_term.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/adjoint_minimum_dissipation_source_term.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/adjoint_mixing_scalar_source_term.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/adjoint_scalar_convection_source_term.f90
)
