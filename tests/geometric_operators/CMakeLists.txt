# Add the neko tests
add_pfunit_ctest(distance
    TEST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/distance.pf
    LINK_LIBRARIES PkgConfig::neko MPI::MPI_Fortran Neko-TOP
)
