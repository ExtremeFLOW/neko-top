# Define the CTest for checkpointing
set(CMAKE_FOLDER "Tests/Unit/Checkpointing")

# ============================================================================ #
# Define reference case information
set(BASE_EXAMPLE "rugby_ball")
set(BASE_CASE_NAME "rugby_ball.case")
set(BASE_EXAMPLE_DIR "${CMAKE_SOURCE_DIR}/examples/${BASE_EXAMPLE}")

# Ensure the Neko binaries can be added to the PATH and the library directories
set(neko_PATH "${neko_PREFIX}/bin:$ENV{PATH}")
string(REPLACE ";" ":"
    neko_LD_LIBRARY_PATH "${neko_LIBRARY_DIRS}:$ENV{LD_LIBRARY_PATH}")

# ============================================================================ #
# Define the test preparation to be run before the test
# This test requires the neko library to be built and installed

# Locate the case file
if (NOT EXISTS "${BASE_EXAMPLE_DIR}/${BASE_CASE_NAME}")
    message(FATAL_ERROR "The case file ${BASE_CASE_NAME} does not exist in "
        "${BASE_EXAMPLE_DIR}.")
endif()
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/checkpoint.case
    COMMENT "Copying case file to the test directory"
    COMMAND
        ${CMAKE_COMMAND} -E copy_if_different
        ${BASE_EXAMPLE_DIR}/${BASE_CASE_NAME}
        ${CMAKE_CURRENT_BINARY_DIR}/checkpoint.case
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_custom_target(checkpointing_case
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/checkpoint.case
)

# If a prepare script exists, copy it to the test directory and run it
if (EXISTS "${BASE_EXAMPLE_DIR}/prepare.sh")
    add_custom_command(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/prepare.sh
        COMMENT "Copying and execute prepare script"
        COMMAND
            ${CMAKE_COMMAND} -E copy_if_different
            ${BASE_EXAMPLE_DIR}/prepare.sh
            ${CMAKE_CURRENT_BINARY_DIR}/prepare.sh
        COMMAND ${CMAKE_COMMAND} -E env
            LD_LIBRARY_PATH="${neko_LD_LIBRARY_PATH}" PATH="${neko_PATH}"
            ./prepare.sh > ${CMAKE_CURRENT_BINARY_DIR}/prepare.log
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    add_custom_target(checkpointing_prepare
        DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/prepare.sh
    )
else()
    message(WARNING "No prepare script found in ${BASE_EXAMPLE_DIR}.")
endif()

# ============================================================================ #
# Define a cleanup step so stale files do not leak between runs

add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/cleanup.done
    COMMENT "Cleaning checkpointing unit test artifacts"
    COMMAND ${CMAKE_COMMAND} -E rm -f
        ${CMAKE_CURRENT_BINARY_DIR}/checkpoint.log
        ${CMAKE_CURRENT_BINARY_DIR}/checkpoint_*.log
        ${CMAKE_CURRENT_BINARY_DIR}/*/checkpoint*.chkp
        ${CMAKE_CURRENT_BINARY_DIR}/*/checkpoint*.h5
        ${CMAKE_CURRENT_BINARY_DIR}/neko.log
    COMMAND ${CMAKE_COMMAND} -E touch ${CMAKE_CURRENT_BINARY_DIR}/cleanup.done
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_custom_target(checkpointing_cleanup
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/cleanup.done
)
add_dependencies(checkpointing_cleanup checkpointing_case)
if (TARGET checkpointing_prepare)
    add_dependencies(checkpointing_cleanup checkpointing_prepare)
endif()

# ============================================================================ #
# Define pFUnit checkpointing tests

# The helper module shared by the tests below, it is built once and linked into
# each of them.
set(CHECKPOINT_UTILS_MOD_DIR
    ${CMAKE_CURRENT_BINARY_DIR}/mod/checkpointing_test_utils)
add_library(checkpointing_test_utils STATIC
    ${CMAKE_CURRENT_SOURCE_DIR}/checkpointing_test_utils.f90
)
set_target_properties(checkpointing_test_utils PROPERTIES
    Fortran_MODULE_DIRECTORY ${CHECKPOINT_UTILS_MOD_DIR}
)
target_include_directories(checkpointing_test_utils
    INTERFACE ${CHECKPOINT_UTILS_MOD_DIR}
)
target_link_libraries(checkpointing_test_utils PRIVATE Neko-TOP::neko-top)

add_pfunit_ctest(checkpoint_restore_consistency
    TEST_SOURCES
        test_checkpoint_restore_consistency.pf
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    LINK_LIBRARIES checkpointing_test_utils Neko-TOP::neko-top
    MAX_PES 2
)

add_pfunit_ctest(checkpoint_end_of_run_gating
    TEST_SOURCES
        test_checkpoint_end_of_run_gating.pf
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    LINK_LIBRARIES checkpointing_test_utils Neko-TOP::neko-top
    MAX_PES 2
)

set(CHECKPOINT_RESTORE_ENV
    "NEKO_LOG_FILE=checkpoint_restore_consistency.log;OMP_NUM_THREADS=1")
set(CHECKPOINT_GATING_ENV
    "NEKO_LOG_FILE=checkpoint_end_of_run_gating.log;OMP_NUM_THREADS=1")

set_tests_properties(checkpoint_restore_consistency PROPERTIES
    LABELS "unit;checkpointing"
    ENVIRONMENT "${CHECKPOINT_RESTORE_ENV}"
    TIMEOUT 300
)


set_tests_properties(checkpoint_end_of_run_gating PROPERTIES
    LABELS "unit;checkpointing"
    ENVIRONMENT "${CHECKPOINT_GATING_ENV}"
    TIMEOUT 300
)

add_dependencies(Unit checkpoint_restore_consistency)
add_dependencies(Unit checkpoint_end_of_run_gating)

add_dependencies(checkpoint_restore_consistency
    checkpointing_case
)
add_dependencies(checkpoint_end_of_run_gating
    checkpointing_case
)

if (TARGET checkpointing_prepare)
    add_dependencies(checkpoint_restore_consistency checkpointing_prepare)
    add_dependencies(checkpoint_end_of_run_gating checkpointing_prepare)
endif()

# Disable treating unused-variable as error for pfunit tests
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(checkpoint_restore_consistency PRIVATE
        $<$<COMPILE_LANGUAGE:Fortran>:-Wno-unused-variable>
    )
    target_compile_options(checkpoint_end_of_run_gating PRIVATE
        $<$<COMPILE_LANGUAGE:Fortran>:-Wno-unused-variable>
    )
endif()

# ============================================================================ #
# Define the HDF5 checkpointing test

if (HAVE_HDF5)
    add_pfunit_ctest(checkpoint_restore_consistency_hdf5
        TEST_SOURCES
            test_checkpoint_restore_consistency_hdf5.pf
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        LINK_LIBRARIES checkpointing_test_utils Neko-TOP::neko-top
        MAX_PES 2
    )
    set(CHECKPOINT_RESTORE_HDF5_ENV
        "NEKO_LOG_FILE=checkpoint_restore_consistency_hdf5.log;OMP_NUM_THREADS=1")
    set_tests_properties(checkpoint_restore_consistency_hdf5 PROPERTIES
        LABELS "unit;checkpointing;hdf5"
        ENVIRONMENT "${CHECKPOINT_RESTORE_HDF5_ENV}"
        TIMEOUT 300
    )
    add_dependencies(Unit checkpoint_restore_consistency_hdf5)
    add_dependencies(checkpoint_restore_consistency_hdf5
        checkpointing_cleanup
        checkpointing_case
    )
    if (TARGET checkpointing_prepare)
        add_dependencies(checkpoint_restore_consistency_hdf5 checkpointing_prepare)
    endif()

    # Disable treating unused-variable as error for pfunit tests
    if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
        target_compile_options(checkpoint_restore_consistency_hdf5 PRIVATE
            $<$<COMPILE_LANGUAGE:Fortran>:-Wno-unused-variable>
        )
    endif()

    # design_checkpoint_hdf5: a design checkpoint must round-trip exactly when
    # a rank owns no design variables. It needs two ranks on every backend and
    # no case file.
    add_pfunit_ctest(design_checkpoint_hdf5
        TEST_SOURCES
            test_design_checkpoint_hdf5.pf
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        LINK_LIBRARIES checkpointing_test_utils Neko-TOP::neko-top
        MAX_PES 2
    )
    set_tests_properties(design_checkpoint_hdf5 PROPERTIES
        LABELS "unit;checkpointing;hdf5"
        ENVIRONMENT "NEKO_LOG_FILE=design_checkpoint_hdf5.log;OMP_NUM_THREADS=1"
        TIMEOUT 120
    )
    add_dependencies(Unit design_checkpoint_hdf5)

    if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
        target_compile_options(design_checkpoint_hdf5 PRIVATE
            $<$<COMPILE_LANGUAGE:Fortran>:-Wno-unused-variable>
        )
    endif()

    # checkpoint_rank_change: the MMA and the design checkpoints must stop with
    # an error on every rank when loaded on fewer or more ranks than the save.
    # checkpoint_rank_change_uneven: the same on as many ranks with another
    # layout, which needs three ranks. The error stops the run, so the driver
    # runs outside pFUnit, from a script that launches it under a time limit.
    add_executable(checkpoint_rank_change EXCLUDE_FROM_ALL
        ${CMAKE_CURRENT_SOURCE_DIR}/checkpoint_rank_change.f90
    )
    target_link_libraries(checkpoint_rank_change
        PRIVATE checkpointing_test_utils Neko-TOP::neko-top
    )
    add_dependencies(Unit checkpoint_rank_change)

    set(rank_change_script
        ${CMAKE_CURRENT_SOURCE_DIR}/check_checkpoint_rank_change.cmake)
    foreach(cases rank_count uneven)
        set(test_name checkpoint_rank_change)
        if (NOT cases STREQUAL "rank_count")
            string(APPEND test_name "_${cases}")
        endif()

        # A directory per test, as both write the same checkpoint file
        set(run_dir ${CMAKE_CURRENT_BINARY_DIR}/run/${test_name})
        file(MAKE_DIRECTORY ${run_dir})
        add_test(NAME ${test_name}
            COMMAND ${CMAKE_COMMAND}
                -DDRIVER=$<TARGET_FILE:checkpoint_rank_change>
                -DMPIEXEC=${MPIEXEC_EXECUTABLE}
                -DNUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}
                "-DPREFLAGS=${MPIEXEC_PREFLAGS}"
                -DCASES=${cases}
                -DTIME_LIMIT=60
                -P ${rank_change_script}
            WORKING_DIRECTORY ${run_dir}
        )
        set_tests_properties(${test_name} PROPERTIES
            LABELS "unit;checkpointing;hdf5"
            ENVIRONMENT "NEKO_LOG_FILE=${test_name}.log;OMP_NUM_THREADS=1"
            TIMEOUT 900
        )
    endforeach()

    if (NOT MPIEXEC_MAX_NUMPROCS OR MPIEXEC_MAX_NUMPROCS LESS 3)
        set_tests_properties(checkpoint_rank_change_uneven PROPERTIES
            DISABLED TRUE
        )
    endif()
endif()
