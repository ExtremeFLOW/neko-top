# CMakeLists.txt for optimizer_checkpoint unit test
# This file requires HDF5 support to be enabled.

if (NOT HAVE_HDF5)
    message(STATUS "Skipping optimizer_checkpoint test since HDF5 is not found")
    return()
endif()

# Define the CTest for optimizer_checkpoint
set(CMAKE_FOLDER "Tests/Unit/OptimizerCheckpoint")

# Ensure the Neko binaries can be added to the PATH and the library directories
set(neko_PATH "${neko_PREFIX}/bin:$ENV{PATH}")
string(REPLACE ";" ":"
    neko_LD_LIBRARY_PATH "${neko_LIBRARY_DIRS}:$ENV{LD_LIBRARY_PATH}")

# ============================================================================ #
# Define the test preparation to be run before the test
# This test requires the neko library to be built and installed

# Prepare case file
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/checkpoint.case
    COMMENT "Copying case file to the test directory"
    COMMAND
        ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_CURRENT_SOURCE_DIR}/checkpoint.case
        ${CMAKE_CURRENT_BINARY_DIR}/checkpoint.case
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_custom_target(optimizer_checkpoint_case
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/checkpoint.case
)

# Prepare script to setup mesh
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/prepare.sh
    COMMENT "Copying and execute prepare script"
    COMMAND
        ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_CURRENT_SOURCE_DIR}/prepare.sh
        ${CMAKE_CURRENT_BINARY_DIR}/prepare.sh
    COMMAND ${CMAKE_COMMAND} -E env
        LD_LIBRARY_PATH="${neko_LD_LIBRARY_PATH}" PATH="${neko_PATH}"
        ${CMAKE_CURRENT_BINARY_DIR}/prepare.sh > ${CMAKE_CURRENT_BINARY_DIR}/prepare.log
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_custom_target(optimizer_checkpoint_prepare
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/prepare.sh
)

# Clear old log and temporary files before running the test
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/cleanup.done
    COMMENT "Cleaning up old log and temporary files"
    COMMAND ${CMAKE_COMMAND} -E rm -f
        ${CMAKE_CURRENT_BINARY_DIR}/optimizer_checkpoint.h5
        ${CMAKE_CURRENT_BINARY_DIR}/optimizer_checkpoint.log
    COMMAND ${CMAKE_COMMAND} -E touch
        ${CMAKE_CURRENT_BINARY_DIR}/cleanup.done
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_custom_target(optimizer_checkpoint_cleanup
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/cleanup.done
)

# ============================================================================ #
# Define the test itself

add_pfunit_ctest(optimizer_checkpoint
    TEST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/optimizer_checkpoint.pf
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    LINK_LIBRARIES Neko-TOP::neko-top
    MAX_PES 2
)

# Set test properties
set_tests_properties(optimizer_checkpoint PROPERTIES
    LABELS "unit;optimizer;hdf5"
    ENVIRONMENT "NEKO_LOG_FILE=optimizer_checkpoint.log"
    TIMEOUT 300
)

# Add dependencies to ensure preparation is done before the test runs
add_dependencies(optimizer_checkpoint
    optimizer_checkpoint_cleanup
    optimizer_checkpoint_prepare
    optimizer_checkpoint_case
)

# Disable treating unused-variable as error for pfunit tests
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(optimizer_checkpoint PRIVATE
        $<$<COMPILE_LANGUAGE:Fortran>:-Wno-unused-variable>
    )
endif()
