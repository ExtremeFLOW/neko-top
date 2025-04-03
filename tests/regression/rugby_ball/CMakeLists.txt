# Define the Rugby Ball regression test


# Import the neko-top static library
add_library(neko-top-a STATIC IMPORTED)
set_target_properties(neko-top-a PROPERTIES
    IMPORTED_LOCATION "${CMAKE_BINARY_DIR}/libneko-top.a"
)


add_executable(Tests-regression-rugby_ball EXCLUDE_FROM_ALL
    ${CMAKE_CURRENT_SOURCE_DIR}/driver.f90
)

# Set the output directory for the tests
set_target_properties(Tests-regression-rugby_ball
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(Tests-regression-rugby_ball
    PkgConfig::neko
    PkgConfig::json-fortran
    neko-top-a
    MPI::MPI_Fortran
)

target_include_directories(Tests-regression-rugby_ball
    PRIVATE
        ${CMAKE_BINARY_DIR}/modules
        ${CMAKE_Fortran_MODULE_DIRECTORY}
)
