target_sources(
    neko-top
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/system_base.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/fixed_point_system.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/periodic_orbit_system.f90
)
