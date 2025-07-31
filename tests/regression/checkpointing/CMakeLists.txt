# Define the CTest for checkpointing

# ============================================================================ #
# Define the executable used for the checkpointing test

set(BASE_EXAMPLE "rugby_ball")
set(BASE_CASE_NAME "rugby_ball.case")
set(BASE_EXAMPLE_DIR "${CMAKE_SOURCE_DIR}/examples/${BASE_EXAMPLE}")

add_executable(checkpointing_test
    ${CMAKE_CURRENT_SOURCE_DIR}/driver.f90
)
target_link_libraries(checkpointing_test
    PRIVATE
        PkgConfig::neko
        PkgConfig::json-fortran
        Neko-TOP::neko-top
        MPI::MPI_Fortran
        $<$<BOOL:${BLAS_FOUND}>:BLAS::BLAS>
        $<$<BOOL:${LAPACK_FOUND}>:LAPACK::LAPACK>
)

# Ensure the Neko binaries can be added to the PATH and the library directories
set(neko_PATH "${neko_PREFIX}/bin:$ENV{PATH}")
string(REPLACE ";" ":"
    neko_LD_LIBRARY_PATH "${neko_LIBRARY_DIRS}:$ENV{LD_LIBRARY_PATH}")


# ============================================================================ #
# Define the test preparation to be run before the test
# This test requires the neko library to be built and installed

# Ensure everything is built before running the test
add_custom_target(prepare_checkpointing_test
    COMMENT "Preparing checkpointing test files"
)
add_dependencies(checkpointing_test
    prepare_checkpointing_test
)

# Locate the case file
if (NOT EXISTS "${BASE_EXAMPLE_DIR}/${BASE_CASE_NAME}")
    message(FATAL_ERROR "The case file ${BASE_CASE_NAME} does not exist in "
        "${BASE_EXAMPLE_DIR}.")
endif()
add_custom_target(
    checkpoint_case_file
    COMMENT "Copying case file to the test directory"
    SOURCES ${BASE_EXAMPLE_DIR}/${BASE_CASE_NAME}
    BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/checkpoint.case
    COMMAND
        ${CMAKE_COMMAND} -E copy_if_different
        ${BASE_EXAMPLE_DIR}/${BASE_CASE_NAME}
        ${CMAKE_CURRENT_BINARY_DIR}/checkpoint.case
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_dependencies(prepare_checkpointing_test
    checkpoint_case_file
)

# If a prepare script exists, copy it to the test directory and run it
if (EXISTS "${BASE_EXAMPLE_DIR}/prepare.sh")
    add_custom_target(
        checkpoint_prepare_script
        COMMENT "Copying and execute prepare script"
        SOURCES ${BASE_EXAMPLE_DIR}/prepare.sh
        COMMAND
            ${CMAKE_COMMAND} -E copy_if_different
            ${BASE_EXAMPLE_DIR}/prepare.sh
            ${CMAKE_CURRENT_BINARY_DIR}/prepare.sh
        COMMAND ${CMAKE_COMMAND} -E env
            LD_LIBRARY_PATH="${neko_LD_LIBRARY_PATH}" PATH="${neko_PATH}"
            ./prepare.sh > ${CMAKE_CURRENT_BINARY_DIR}/prepare.log
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    add_dependencies(prepare_checkpointing_test
        checkpoint_prepare_script
    )
else()
    message(WARNING "No prepare script found in ${BASE_EXAMPLE_DIR}.")
endif()

# ============================================================================ #
# Define the test itself

# Add the test to CTest
add_test(NAME checkpointing_test
    COMMAND $<TARGET_FILE:checkpointing_test> checkpoint.case
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
# Set the environment variables for the test
set_tests_properties(checkpointing_test PROPERTIES
    ENVIRONMENT "NEKO_LOG_FILE=neko.log"
)
