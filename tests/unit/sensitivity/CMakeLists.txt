# Define the CTest for sensitivity

# Set the CMake folder for organization
set(CMAKE_FOLDER "Tests/Unit/Sensitivity")

# Define the list of sensitivity tests
# each of these should be a custom driver.
set(test_list
    "minimum_dissipation.f90"
)

# ============================================================================ #
# Ensure the neko library is found correctly

# Get the Neko-TOP module directory
get_target_property(NEKO_TOP_MODDIR neko-top Fortran_MODULE_DIRECTORY)

# ============================================================================ #
# Define the test preparation to be run before the test
# This test requires the neko library to be built and installed

# Ensure everything is built before running the test
add_custom_target(sensitivity_test_prepare
    COMMENT "Preparing sensitivity test files")

# Ensure the sensitivity utils library is built
add_library(sensitivity_utils EXCLUDE_FROM_ALL
    ${CMAKE_CURRENT_SOURCE_DIR}/sensitivity.f90)
target_include_directories(sensitivity_utils PUBLIC ${NEKO_TOP_MODDIR})
target_link_libraries(sensitivity_utils
    PUBLIC
        neko-top-a
        PkgConfig::neko
        PkgConfig::json-fortran
        MPI::MPI_Fortran
        $<$<BOOL:${BLAS_FOUND}>:BLAS::BLAS>
        $<$<BOOL:${LAPACK_FOUND}>:LAPACK::LAPACK>
)
add_dependencies(sensitivity_test_prepare sensitivity_utils)

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

# Ensure the test preparation is run before the test
add_dependencies(Unit sensitivity_test_prepare)
add_dependencies(sensitivity_test_prepare sensitivity_case_file)
add_dependencies(sensitivity_test_prepare sensitivity_prepare_script)

# ============================================================================ #
# Define the test itself

foreach(test ${test_list})
    # Extract the test name from the file name
    get_filename_component(test_name ${test} NAME_WE)
    add_executable(${test_name}_sensitivity_bin EXCLUDE_FROM_ALL
        ${CMAKE_CURRENT_SOURCE_DIR}/${test}
    )
    target_link_libraries(${test_name}_sensitivity_bin
        sensitivity_utils
    )
    add_dependencies(Unit ${test_name}_sensitivity_bin)

    # Add the test to CTest
    add_test(NAME ${test_name}_sensitivity
        COMMAND $<TARGET_FILE:${test_name}_sensitivity_bin> sensitivity.case
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    set_tests_properties(${test_name}_sensitivity
        PROPERTIES
            DEPENDS "${test_name}_sensitivity_bin;sensitivity_test_prepare"
            LABELS "sensitivity;unit"
            ENVIRONMENT "NEKO_LOG_FILE=neko.log"
            TIMEOUT 5000
    )
    set_target_properties(${test_name}_sensitivity_bin PROPERTIES
        FOLDER "Tests/Unit/Sensitivity"
    )
endforeach()
