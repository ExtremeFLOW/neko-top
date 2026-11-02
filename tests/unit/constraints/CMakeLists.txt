# CMakeLists.txt for the constraint unit tests.
#
# volume_constraint_free: a volume constraint initialised from JSON with a
# PDE_filter mapping must drop that mapping when it is freed or re-initialised
# from components, and then match a constraint that never had one.

set(CMAKE_FOLDER "Tests/Unit/Constraints")

# Ensure the Neko binaries can be added to the PATH and the library directories
set(neko_PATH "${neko_PREFIX}/bin:$ENV{PATH}")
string(REPLACE ";" ":"
    neko_LD_LIBRARY_PATH "${neko_LIBRARY_DIRS}:$ENV{LD_LIBRARY_PATH}")

set(VOLUME_FREE_DIR ${CMAKE_CURRENT_BINARY_DIR}/run/volume_constraint_free)
file(MAKE_DIRECTORY ${VOLUME_FREE_DIR})

# ============================================================================ #
# Prepare the case file and the mesh

add_custom_command(
    OUTPUT ${VOLUME_FREE_DIR}/volume_constraint_free.case
    COMMENT "Copying case file to the test directory"
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/volume_constraint_free.case
    COMMAND
        ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_CURRENT_SOURCE_DIR}/volume_constraint_free.case
        ${VOLUME_FREE_DIR}/volume_constraint_free.case
)

add_custom_command(
    OUTPUT ${VOLUME_FREE_DIR}/box.nmsh
    COMMENT "Generating the mesh for volume_constraint_free"
    COMMAND ${CMAKE_COMMAND} -E env
        LD_LIBRARY_PATH="${neko_LD_LIBRARY_PATH}" PATH="${neko_PATH}"
        genmeshbox 0 1 0 1 0 0.125 8 8 1 .false. .true. .true.
        > ${VOLUME_FREE_DIR}/genmeshbox.log
    WORKING_DIRECTORY ${VOLUME_FREE_DIR}
)

add_custom_target(volume_constraint_free_prepare
    DEPENDS
        ${VOLUME_FREE_DIR}/volume_constraint_free.case
        ${VOLUME_FREE_DIR}/box.nmsh
)

# ============================================================================ #
# Define the test itself

add_pfunit_ctest(volume_constraint_free
    TEST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/volume_constraint_free.pf
    WORKING_DIRECTORY ${VOLUME_FREE_DIR}
    LINK_LIBRARIES Neko-TOP::neko-top
    MAX_PES 1
)

add_dependencies(volume_constraint_free volume_constraint_free_prepare)
add_dependencies(Unit volume_constraint_free)

set_tests_properties(volume_constraint_free PROPERTIES
    LABELS "unit;constraints"
    ENVIRONMENT "NEKO_LOG_FILE=volume_constraint_free.log"
    TIMEOUT 120
)

# Disable treating unused-variable as error for pfunit tests
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(volume_constraint_free PRIVATE
        $<$<COMPILE_LANGUAGE:Fortran>:-Wno-unused-variable>
    )
endif()
