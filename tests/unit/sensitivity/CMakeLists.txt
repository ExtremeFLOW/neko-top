# Define the CTest for sensitivity
set(CMAKE_FOLDER "Tests/Unit/Sensitivity")

# Provide a list of .case files
set(test_list
    "volume.case"
    "volume_filtered.case"
    "viscous_dissipation.case"
    "brinkman_dissipation.case"
    "viscous_dissipation_variable_dt.case"
    "viscous_dissipation_variable_dt_weighted.case"
    "unsteady_brinkman_dissipation.case"
    "unsteady_brinkman_dissipation_masked.case"
)

# Path to the common driver source
set(SENSITIVITY_DRIVER ${CMAKE_CURRENT_SOURCE_DIR}/problem_tester.f90)

# ============================================================================ #
# Use a private .mod folder to avoid clashes with regression targets
set(LOCAL_MOD_DIR ${CMAKE_CURRENT_BINARY_DIR}/modules)
set(CMAKE_Fortran_MODULE_DIRECTORY ${LOCAL_MOD_DIR} CACHE PATH
    "Fortran module directory" FORCE)

# ============================================================================ #
# Preparation target
add_custom_target(sensitivity_test_prepare
    COMMENT "Preparing sensitivity test files"
)

# ============================================================================ #
# Neko environment
set(neko_PATH "${neko_PREFIX}/bin:$ENV{PATH}")
string(REPLACE ";" ":"
    neko_LD_LIBRARY_PATH "${neko_LIBRARY_DIRS}:$ENV{LD_LIBRARY_PATH}"
)

# Prepare script (if exists)
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
add_dependencies(sensitivity_test_prepare sensitivity_prepare_script)

# ============================================================================ #
# Build the single driver executable
add_executable(problem_tester_sensitivity_bin EXCLUDE_FROM_ALL
    # fd_criterion is listed first: sensitivity uses it, and the Fortran
    # module it defines has to exist by then.
    ${SHARED_TEST_SOURCES}/fd_criterion.f90
    ${SHARED_TEST_SOURCES}/sensitivity.f90
    ${SENSITIVITY_DRIVER}
)
target_link_libraries(problem_tester_sensitivity_bin
    PRIVATE Neko-TOP::neko-top
)
target_include_directories(problem_tester_sensitivity_bin
    PRIVATE ${LOCAL_MOD_DIR})
add_dependencies(Unit problem_tester_sensitivity_bin)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules CACHE PATH
    "Fortran module directory" FORCE)

# ============================================================================ #
# Loop over all case files and create tests. The time limit is set outside the
# loop because the two-rank test after it uses the same one.
set(sensitivity_case_timeout 300)
foreach(case_file ${test_list})
    # Get the case name without extension
    get_filename_component(case_name ${case_file} NAME_WE)

    # Copy case file
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
    add_dependencies(sensitivity_test_prepare ${case_name}_case_file)

    # Add test for this case
    add_test(NAME ${case_name}_sensitivity
        COMMAND $<TARGET_FILE:problem_tester_sensitivity_bin> ${case_file}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    # 78 is the driver's "the assertion was skipped, not made" exit code
    # (FD_SKIP_EXIT_CODE in tests/shared/sensitivity.f90), returned by a
    # single-precision build, where the functional is not reproducible
    # enough for any finite-difference verdict to mean anything. Without
    # this the driver exits 0 and CTest reports a green gradient check that
    # checked nothing. 77 is deliberately avoided: the regression lane uses
    # it for its opt-in gate, and the two reasons must stay distinguishable.
    set_tests_properties(${case_name}_sensitivity
        PROPERTIES
            DEPENDS "problem_tester_sensitivity_bin;sensitivity_test_prepare"
            LABELS "sensitivity;unit"
            ENVIRONMENT "NEKO_LOG_FILE=neko.log"
            SKIP_RETURN_CODE 78
            TIMEOUT ${sensitivity_case_timeout}
    )
endforeach()

# ============================================================================ #
# The cases above run on a single rank, where the whole design sits on one
# rank and the probe cannot be nominated more than once. Run one of them on two
# ranks as well, so that the probe selection and the cross-rank reductions of
# the sweep are checked under MPI too. The case is copied under its own name
# because the CSV logs are appended to, not overwritten: this keeps the
# two-rank rows from being mixed into the single-rank run's CSVs.
set(case_name "viscous_dissipation_variable_dt")
set(mpi_case_name "${case_name}_np2")
add_custom_target(
    ${mpi_case_name}_case_file
    COMMENT "Copying case file ${case_name}.case to ${mpi_case_name}.case"
    SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/${case_name}.case
    BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/${mpi_case_name}.case
    COMMAND
        ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_CURRENT_SOURCE_DIR}/${case_name}.case
        ${CMAKE_CURRENT_BINARY_DIR}/${mpi_case_name}.case
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_dependencies(sensitivity_test_prepare ${mpi_case_name}_case_file)

# Like the multi-rank pFUnit tests in tests/unit/CMakeLists.txt, a device
# backend runs it on one rank only.
set(mpi_case_ranks $<IF:$<NOT:$<STREQUAL:${DEVICE_TYPE},CPU>>,1,2>)
add_test(NAME ${mpi_case_name}_sensitivity
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${mpi_case_ranks}
        ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:problem_tester_sensitivity_bin> ${MPIEXEC_POSTFLAGS}
        ${mpi_case_name}.case
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(${mpi_case_name}_sensitivity
    PROPERTIES
        LABELS "sensitivity;unit"
        ENVIRONMENT "NEKO_LOG_FILE=${mpi_case_name}.log"
        PROCESSORS ${mpi_case_ranks}
        SKIP_RETURN_CODE 78
        TIMEOUT ${sensitivity_case_timeout}
)

# ============================================================================ #
# A steady adjoint takes as many steps as the forward run, with time steps of
# its own, so its clock can stop short of end_time. It must still count as
# complete, and write no job-limit checkpoint.
set(case_name "adjoint_joblimit_variable_dt")
add_executable(adjoint_joblimit_tester_bin EXCLUDE_FROM_ALL
    ${CMAKE_CURRENT_SOURCE_DIR}/adjoint_joblimit_tester.f90
)
target_link_libraries(adjoint_joblimit_tester_bin
    PRIVATE Neko-TOP::neko-top
)
add_dependencies(Unit adjoint_joblimit_tester_bin)

add_custom_target(
    ${case_name}_case_file
    COMMENT "Copying case file ${case_name}.case to the test directory"
    SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/${case_name}.case
    BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/${case_name}.case
    COMMAND
        ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_CURRENT_SOURCE_DIR}/${case_name}.case
        ${CMAKE_CURRENT_BINARY_DIR}/${case_name}.case
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_dependencies(sensitivity_test_prepare ${case_name}_case_file)

add_test(NAME ${case_name}
    COMMAND $<TARGET_FILE:adjoint_joblimit_tester_bin> ${case_name}.case
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(${case_name}
    PROPERTIES
        DEPENDS "adjoint_joblimit_tester_bin;sensitivity_test_prepare"
        LABELS "sensitivity;unit"
        ENVIRONMENT "NEKO_LOG_FILE=${case_name}.log"
        TIMEOUT 300
)

# ============================================================================ #
# Final dependencies
add_dependencies(Unit sensitivity_test_prepare)
