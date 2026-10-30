# Define the CTests for objective time windows
set(CMAKE_FOLDER "Tests/Unit/Objectives")

# Provide a list of .case files
set(test_list
    "time_window_run_length.case"
    "time_window_equivalence.case"
    "steady_unsteady_converged.case"
    "steady_unsteady_final_step.case"
    "steady_scalar_unconverged.case"
    "steady_scalar_unconverged_simtime.case"
)

# Path to the common driver source
set(OBJECTIVES_DRIVER ${CMAKE_CURRENT_SOURCE_DIR}/time_window_tester.f90)

# pFUnit source for the run-length check, which reuses
# `time_window_run_length.case`
set(STEP_COUNT_TEST ${CMAKE_CURRENT_SOURCE_DIR}/step_count.pf)
set(STEP_COUNT_CASE "time_window_run_length.case")

# ============================================================================ #
# Use a private .mod folder to avoid clashes with other test targets
set(LOCAL_MOD_DIR ${CMAKE_CURRENT_BINARY_DIR}/modules)
set(CMAKE_Fortran_MODULE_DIRECTORY ${LOCAL_MOD_DIR} CACHE PATH
    "Fortran module directory" FORCE)

# ============================================================================ #
# Preparation target
add_custom_target(objectives_test_prepare
    COMMENT "Preparing objective time-window test files"
)

# ============================================================================ #
# Neko environment
set(neko_PATH "${neko_PREFIX}/bin:$ENV{PATH}")
string(REPLACE ";" ":"
    neko_LD_LIBRARY_PATH "${neko_LIBRARY_DIRS}:$ENV{LD_LIBRARY_PATH}"
)

# Mesh generation
add_custom_target(
    objectives_prepare_script
    COMMENT "Copying and executing the objective test prepare script"
    SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/prepare.sh
    BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/prepare.sh
    COMMAND
        ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_CURRENT_SOURCE_DIR}/prepare.sh
        ${CMAKE_CURRENT_BINARY_DIR}/prepare.sh
    COMMAND ${CMAKE_COMMAND} -E env
        LD_LIBRARY_PATH="${neko_LD_LIBRARY_PATH}" PATH="${neko_PATH}"
        ./prepare.sh -N4 > ${CMAKE_CURRENT_BINARY_DIR}/prepare.log
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_dependencies(objectives_test_prepare objectives_prepare_script)

# ============================================================================ #
# The user module carrying the inflow every case here uses. It is built once
# and linked into both drivers; compiling it into each would have the two
# targets race on the same .mod file.
add_library(objectives_user_obj OBJECT EXCLUDE_FROM_ALL
    ${CMAKE_CURRENT_SOURCE_DIR}/objectives_user.f90
)
target_link_libraries(objectives_user_obj PRIVATE Neko-TOP::neko-top)

# ============================================================================ #
# Build the time-window driver executable
add_executable(time_window_tester_bin EXCLUDE_FROM_ALL ${OBJECTIVES_DRIVER})
target_link_libraries(time_window_tester_bin
    PRIVATE objectives_user_obj Neko-TOP::neko-top
)
target_include_directories(time_window_tester_bin PRIVATE ${LOCAL_MOD_DIR})
add_dependencies(Unit time_window_tester_bin)

add_pfunit_ctest(end_time_step_count
    TEST_SOURCES ${STEP_COUNT_TEST}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    LINK_LIBRARIES objectives_user_obj Neko-TOP::neko-top
    MAX_PES 1
)
target_include_directories(end_time_step_count PRIVATE ${LOCAL_MOD_DIR})
add_dependencies(Unit end_time_step_count)
add_dependencies(end_time_step_count objectives_test_prepare)

# Disable treating pFUnit's generated variables as unused.
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(end_time_step_count PRIVATE
        $<$<COMPILE_LANGUAGE:Fortran>:-Wno-unused-variable>
    )
endif()

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules CACHE PATH
    "Fortran module directory" FORCE)

# ============================================================================ #
# Loop over all case files and create tests
foreach(case_file ${test_list})
    get_filename_component(case_name ${case_file} NAME_WE)

    add_custom_target(
        ${case_name}_case_file
        COMMENT "Copying case file ${case_file} to the test directory"
        SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/${case_file}
        BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/${case_file}
        COMMAND
            ${CMAKE_COMMAND} -E copy_if_different
            ${CMAKE_CURRENT_SOURCE_DIR}/${case_file}
            ${CMAKE_CURRENT_BINARY_DIR}/${case_file}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    add_dependencies(objectives_test_prepare ${case_name}_case_file)

    add_test(NAME ${case_name}
        COMMAND $<TARGET_FILE:time_window_tester_bin> ${case_file}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    set_tests_properties(${case_name}
        PROPERTIES
            DEPENDS "time_window_tester_bin;objectives_test_prepare"
            LABELS "objectives;unit"
            ENVIRONMENT "NEKO_LOG_FILE=${case_name}.log"
            TIMEOUT 300
    )
endforeach()

# ============================================================================ #
# Check what the `steady` component reported about the scalar, from the log a
# case's test wrote. The fixture runs that test first.
function(add_steady_log_check case_name frozen at_end converged end_time)
    add_test(NAME ${case_name}_steady_log
        COMMAND ${CMAKE_COMMAND}
            -DLOG_FILE=${CMAKE_CURRENT_BINARY_DIR}/${case_name}.log
            -DFROZEN=${frozen} -DAT_END=${at_end} -DCONVERGED=${converged}
            -DEND_TIME=${end_time}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/check_steady_log.cmake
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    set_tests_properties(${case_name} PROPERTIES
        FIXTURES_SETUP ${case_name}_log
    )
    set_tests_properties(${case_name}_steady_log PROPERTIES
        FIXTURES_REQUIRED ${case_name}_log
        LABELS "objectives;unit"
    )
endfunction()

# Each case runs twice, through the steady path and the unsteady one. The
# unconverged ones stop with the scalar unconverged, so each run warns at the
# freeze and at its end. The converged one runs on until the scalar has
# converged, so each run warns at the freeze and then reports that the scalar
# converged. The last argument is the time of the component's last check in
# each run, which the end-of-run report must carry. It is the final step under
# `tsteps` control, and under `simulationtime` the last scheduled time before
# the end: 2.4, as 2.5 lies past `end_time = 2.495`.
add_steady_log_check(steady_scalar_unconverged 2 2 0 1.50000E+00)
add_steady_log_check(steady_unsteady_converged 2 0 2 2.00000E+01)
add_steady_log_check(steady_scalar_unconverged_simtime 2 2 0 2.40000E+00)

# ============================================================================ #
# Run-length check: a run whose end_time lies on a step boundary must stop
# there, on the steady and the unsteady path alike. The case file is copied
# by the loop above.
set_tests_properties(end_time_step_count
    PROPERTIES
        LABELS "objectives;unit"
        ENVIRONMENT "NEKO_LOG_FILE=end_time_step_count.log"
        TIMEOUT 300
)

# ============================================================================ #
# Final dependencies
add_dependencies(Unit objectives_test_prepare)
