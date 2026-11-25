# CMakeLists.txt for the design unit tests.
#
# design_3dto1d_write: the file written by design_3dto1d_t%write must hold
# the values of every rank, as held on the device on a device backend. It
# needs two ranks on every backend, so it is registered here rather than in
# the shared list in the parent CMakeLists.txt, which runs device builds with
# a single rank.
#
# brinkman_design_no_mapping: a Brinkman design without a "mapping" must map
# a continuous design to the Brinkman amplitude, and a continuous sensitivity
# back to the design, by the identity, whether it is initialised from JSON or
# from components. Freeing a design must also drop its mask, and on a
# discontinuous design the empty cascade must be a projection whose backward
# map is its transpose.

set(CMAKE_FOLDER "Tests/Unit/Design")

# Ensure the Neko binaries can be added to the PATH and the library directories
set(neko_PATH "${neko_PREFIX}/bin:$ENV{PATH}")
string(REPLACE ";" ":"
    neko_LD_LIBRARY_PATH "${neko_LIBRARY_DIRS}:$ENV{LD_LIBRARY_PATH}")

# ============================================================================ #
# design_3dto1d_write

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

# ============================================================================ #
# brinkman_design_no_mapping: prepare the case file and the mesh

set(NO_MAPPING_DIR ${CMAKE_CURRENT_BINARY_DIR}/run/brinkman_design_no_mapping)
file(MAKE_DIRECTORY ${NO_MAPPING_DIR})

add_custom_command(
    OUTPUT ${NO_MAPPING_DIR}/brinkman_design_no_mapping.case
    COMMENT "Copying case file to the test directory"
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/brinkman_design_no_mapping.case
    COMMAND
        ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_CURRENT_SOURCE_DIR}/brinkman_design_no_mapping.case
        ${NO_MAPPING_DIR}/brinkman_design_no_mapping.case
)

add_custom_command(
    OUTPUT ${NO_MAPPING_DIR}/box.nmsh
    COMMENT "Generating the mesh for brinkman_design_no_mapping"
    COMMAND ${CMAKE_COMMAND} -E env
        LD_LIBRARY_PATH="${neko_LD_LIBRARY_PATH}" PATH="${neko_PATH}"
        genmeshbox 0 1 0 1 0 0.125 8 8 1 .false. .true. .true.
        > ${NO_MAPPING_DIR}/genmeshbox.log
    WORKING_DIRECTORY ${NO_MAPPING_DIR}
)

add_custom_target(brinkman_design_no_mapping_prepare
    DEPENDS
        ${NO_MAPPING_DIR}/brinkman_design_no_mapping.case
        ${NO_MAPPING_DIR}/box.nmsh
)

# ============================================================================ #
# brinkman_design_no_mapping: define the test itself

add_pfunit_ctest(brinkman_design_no_mapping
    TEST_SOURCES brinkman_design_no_mapping.pf
    WORKING_DIRECTORY ${NO_MAPPING_DIR}
    LINK_LIBRARIES Neko-TOP::neko-top
    MAX_PES 1
)

add_dependencies(brinkman_design_no_mapping brinkman_design_no_mapping_prepare)
add_dependencies(Unit brinkman_design_no_mapping)

set_tests_properties(brinkman_design_no_mapping PROPERTIES
    LABELS "unit;design"
    ENVIRONMENT "NEKO_LOG_FILE=brinkman_design_no_mapping.log"
    TIMEOUT 120
)

# Disable treating unused-variable as error for pfunit tests
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    foreach(test design_3dto1d_write brinkman_design_no_mapping)
        target_compile_options(${test} PRIVATE
            $<$<COMPILE_LANGUAGE:Fortran>:-Wno-unused-variable>
        )
    endforeach()
endif()
