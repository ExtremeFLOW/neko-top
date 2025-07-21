# Define the CTest for sensitivity

# Set the CMake folder for organization
set(CMAKE_FOLDER "Tests/Regression/Sensitivity")

# ============================================================================ #
# Ensure the neko library is found correctly

# Get the Neko-TOP module directory
get_target_property(NEKO_TOP_MODDIR neko-top Fortran_MODULE_DIRECTORY)
include_directories(${NEKO_TOP_MODDIR})
link_libraries(
    neko-top-a
    PkgConfig::neko
    PkgConfig::json-fortran
    MPI::MPI_Fortran
    $<$<BOOL:${BLAS_FOUND}>:BLAS::BLAS>
    $<$<BOOL:${LAPACK_FOUND}>:LAPACK::LAPACK>
)

# ============================================================================ #
# Define the executable used for the sensitivity tests

# ---------------------------------------------------------------------------- #
# Add the volume test executable

add_executable(volume_test_bin EXCLUDE_FROM_ALL
    ${CMAKE_CURRENT_SOURCE_DIR}/volume.f90)

# ============================================================================ #
# Define the test preparation to be run before the test
# This test requires the neko library to be built and installed

# Ensure everything is built before running the test
add_custom_target(sensitivity_test_prepare
    COMMENT "Preparing sensitivity test files")

# Locate the case file
add_custom_target(
    sensitivity_case_file
    COMMENT "Copying case file to the test directory"
    SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/sensitivity.case
    BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/sensitivity.case
    COMMAND
        ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_CURRENT_SOURCE_DIR}/sensitivity.case
        ${CMAKE_CURRENT_BINARY_DIR}/sensitivity.case
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

# Ensure the Neko binaries can be added to the PATH and the library directories
set(neko_PATH "${neko_PREFIX}/bin:$ENV{PATH}")
string(REPLACE ";" ":"
    neko_LD_LIBRARY_PATH "${neko_LIBRARY_DIRS}:$ENV{LD_LIBRARY_PATH}")

# If a prepare script exists, copy it to the test directory and run it
add_custom_target(
    sensitivity_prepare_script
    COMMENT "Copying and execute prepare script"
    SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/prepare.sh
    BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/prepare.sh
    COMMAND
        ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_CURRENT_SOURCE_DIR}/prepare.sh
        ${CMAKE_CURRENT_BINARY_DIR}/prepare.sh
    COMMAND ${CMAKE_COMMAND} -E env
        LD_LIBRARY_PATH="${neko_LD_LIBRARY_PATH}" PATH="${neko_PATH}"
        ./prepare.sh -N8 > ${CMAKE_CURRENT_BINARY_DIR}/prepare.log
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

# ============================================================================ #
# Define the dependency list

add_dependencies(Regression volume_test_bin)

# Ensure the test preparation is run before the test
add_dependencies(Regression sensitivity_test_prepare)
add_dependencies(sensitivity_test_prepare sensitivity_case_file)
add_dependencies(sensitivity_test_prepare sensitivity_prepare_script)

# ============================================================================ #
# Define the test itself

# Add the test to CTest
add_test(NAME volume_test
    COMMAND $<TARGET_FILE:volume_test_bin> sensitivity.case
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(volume_test
    PROPERTIES
        DEPENDS "volume_test_bin;sensitivity_test_prepare"
        LABELS "sensitivity;regression"
        ENVIRONMENT "NEKO_LOG_FILE=neko.log"
        TIMEOUT 300
)
