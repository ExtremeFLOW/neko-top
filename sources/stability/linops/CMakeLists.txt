target_sources(
    neko-top
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/linop_base.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/linop_fixed_baseflow.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/linop_transient_baseflow.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/neko_resolvent.f90
)
