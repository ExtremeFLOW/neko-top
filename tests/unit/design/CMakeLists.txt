# CMakeLists.txt for the design unit tests.
#
# design_3dto1d_write: the file written by design_3dto1d_t%write must hold
# the values of every rank, as held on the device on a device backend. It
# needs two ranks on every backend, so it is registered here rather than in
# the shared list in the parent CMakeLists.txt, which runs device builds with
# a single rank.

set(CMAKE_FOLDER "Tests/Unit/Design")

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/run/design_3dto1d_write)

add_pfunit_ctest(design_3dto1d_write
    TEST_SOURCES design_3dto1d_write.pf
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/run/design_3dto1d_write
    LINK_LIBRARIES Neko-TOP::neko-top
    MAX_PES 2
)
add_dependencies(Unit design_3dto1d_write)
set_tests_properties(design_3dto1d_write PROPERTIES
    LABELS "unit;design"
    ENVIRONMENT "NEKO_LOG_FILE=design_3dto1d_write.log"
    TIMEOUT 120
)

# Disable treating unused-variable as error for pfunit tests
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(design_3dto1d_write PRIVATE
        $<$<COMPILE_LANGUAGE:Fortran>:-Wno-unused-variable>
    )
endif()
