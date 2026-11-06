# Define the CTest for sensitivity (regression suite)
set(CMAKE_FOLDER "Tests/Regression/Sensitivity")

# --------------------------------------------------------------------------- #
# Ensure the neko library module dir is found (from the main target)
get_target_property(NEKO_TOP_MOD_DIR neko-top Fortran_MODULE_DIRECTORY)

# --------------------------------------------------------------------------- #
# Preparation target: copies helper scripts into the binary dir
add_custom_target(reg_sensitivity_test_prepare ALL
    COMMENT "Preparing regression sensitivity test files"
)

# --------------------------------------------------------------------------- #
# Use a private .mod folder to avoid clashes with unit targets
set(LOCAL_MOD_DIR ${CMAKE_CURRENT_BINARY_DIR}/modules)
set(CMAKE_Fortran_MODULE_DIRECTORY ${LOCAL_MOD_DIR} CACHE PATH
    "Fortran module directory" FORCE)

# --------------------------------------------------------------------------- #
# Build the single driver executable for regression
add_executable(reg_problem_tester_sensitivity_bin
    ${CMAKE_CURRENT_SOURCE_DIR}/user.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/problem_tester.f90

    # Shared sensitivity module.
    ${SHARED_TEST_SOURCES}/sensitivity.f90
)
target_link_libraries(reg_problem_tester_sensitivity_bin
    PRIVATE Neko-TOP::neko-top
)

# The driver must also see the regression .mod directory
target_include_directories(reg_problem_tester_sensitivity_bin
    PRIVATE ${LOCAL_MOD_DIR})

set_target_properties(reg_problem_tester_sensitivity_bin
    PROPERTIES
        FOLDER "Tests/Regression/Sensitivity"
        RUNTIME_OUTPUT_NAME "sensitivity_regression_driver"
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

# Let the umbrella target trigger both build+prepare steps (if present)
add_dependencies(Regression reg_problem_tester_sensitivity_bin)

set(CMAKE_Fortran_MODULE_DIRECTORY ${NEKO_TOP_MOD_DIR} CACHE PATH
    "Fortran module directory" FORCE)

# --------------------------------------------------------------------------- #
# Register the regression sensitivity checks, split into two tiers.
#
# Each case compares the adjoint gradient with a one-sided finite-difference
# sweep over four perturbations, which takes six solves: the forward and
# adjoint solves at the baseline design plus one forward solve per
# perturbation. The driver itself asserts the finite-difference error
# against the case's tolerance (see run_sensitivity_regression_check.sh).
#
# - Fast tier: the unsteady cases, 5 timesteps per solve. They need no
#   opt-in.
# - Slow tier: the steady-state cases, 1000 timesteps per solve (10000 for
#   passive_scalar). They are not part of the default test budget and only
#   run when NEKO_TOP_RUN_SENSITIVITY_REGRESSION is set; otherwise
#   run_sensitivity_regression_check.sh exits 77, the SKIP_RETURN_CODE set
#   below, so CTest reports them as skipped. This mirrors the opt-in
#   pattern already used for tests/regression/valgrind (there gated on
#   `valgrind` being installed; here gated on an explicit opt-in since the
#   tool is always installed).
#
# With the default two MPI ranks on an 8-core laptop (CPU backend), a fast
# case takes a few seconds, passive_scalar about 30 minutes and each other
# slow case 4-7 minutes. The TIMEOUTs below allow about three times that or
# more: 300 s for the fast tier, 5400 s for passive_scalar and 1800 s for
# the rest.
#
# All cases share one RESOURCE_LOCK because they write into the same
# WORKING_DIRECTORY and reuse one generated mesh (see
# run_sensitivity_regression_check.sh), so ctest never runs two of them at
# once, regardless of -j.

# Runtime environment for Neko tools (genmeshbox) and libraries.
set(neko_PATH "${neko_PREFIX}/bin:$ENV{PATH}")
string(REPLACE ";" ":" neko_LD_LIBRARY_PATH
    "${neko_LIBRARY_DIRS}:$ENV{LD_LIBRARY_PATH}")

add_test(NAME sensitivity_regression_build
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR}
        --target reg_problem_tester_sensitivity_bin
)
set_tests_properties(sensitivity_regression_build PROPERTIES
    LABELS "regression;sensitivity"
    TIMEOUT 900
)

# Fast tier. The test command sets NEKO_TOP_RUN_SENSITIVITY_REGRESSION=1
# itself, so these cases run whatever the caller's environment and are
# never skipped.
set(sensitivity_regression_fast_cases
    "unsteady_viscous_dissipation"
    "unsteady_brinkman_dissipation"
)

# Slow tier. passive_scalar runs ten times as many timesteps per solve as
# the other slow cases (end_time 10.0 vs 1.0 at the same timestep), so it
# gets a larger TIMEOUT rather than a shared one sized for the cheaper
# cases.
set(sensitivity_regression_slow_cases
    "viscous_dissipation"
    "brinkman_dissipation"
    "dissipation"
    "dissipation_weights"
    "passive_scalar"
)

foreach(case_name ${sensitivity_regression_fast_cases})
    add_test(NAME ${case_name}_sensitivity_regression
        COMMAND ${CMAKE_COMMAND} -E env
            PATH=${neko_PATH}
            LD_LIBRARY_PATH=${neko_LD_LIBRARY_PATH}
            NEKO_DIR=${neko_PREFIX}
            NEKO_TOP_RUN_SENSITIVITY_REGRESSION=1
            bash ${CMAKE_CURRENT_SOURCE_DIR}/run_sensitivity_regression_check.sh
            $<TARGET_FILE:reg_problem_tester_sensitivity_bin>
            ${CMAKE_CURRENT_SOURCE_DIR}/cases/${case_name}.case
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    set_tests_properties(${case_name}_sensitivity_regression PROPERTIES
        DEPENDS "sensitivity_regression_build"
        LABELS "regression;sensitivity"
        RESOURCE_LOCK "sensitivity_regression_workspace"
        TIMEOUT 300
    )
endforeach()

foreach(case_name ${sensitivity_regression_slow_cases})
    add_test(NAME ${case_name}_sensitivity_regression
        COMMAND ${CMAKE_COMMAND} -E env
            PATH=${neko_PATH}
            LD_LIBRARY_PATH=${neko_LD_LIBRARY_PATH}
            NEKO_DIR=${neko_PREFIX}
            bash ${CMAKE_CURRENT_SOURCE_DIR}/run_sensitivity_regression_check.sh
            $<TARGET_FILE:reg_problem_tester_sensitivity_bin>
            ${CMAKE_CURRENT_SOURCE_DIR}/cases/${case_name}.case
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    if(case_name STREQUAL "passive_scalar")
        set(case_timeout 5400)
    else()
        set(case_timeout 1800)
    endif()
    set_tests_properties(${case_name}_sensitivity_regression PROPERTIES
        DEPENDS "sensitivity_regression_build"
        LABELS "regression;sensitivity;sensitivity_slow"
        RESOURCE_LOCK "sensitivity_regression_workspace"
        SKIP_RETURN_CODE 77
        TIMEOUT ${case_timeout}
    )
endforeach()

# --------------------------------------------------------------------------- #
# Report a slow tier that was not opted in as a failure.
#
# CTest's summary line ("X% tests passed, Y tests failed out of Z") does not
# count a SKIP_RETURN_CODE test as failed, so a run with the slow tier
# skipped reads as a full pass unless the "did not run" list below it is
# read too. This test puts the opt-in state into the pass/fail count: it
# passes when NEKO_TOP_RUN_SENSITIVITY_REGRESSION is set to a non-empty
# value and otherwise fails, under a name that says why, printing how to
# opt in. It is labelled "sensitivity" and "sensitivity_slow" but not
# "regression", so `ctest -L regression` never selects it, while
# `ctest -L sensitivity` and an unfiltered `ctest` do. Those two therefore
# exit non-zero unless the slow tier is opted in, and opting in runs the
# whole slow tier, about 50 minutes by the timings above.
# `ctest -LE sensitivity_slow` runs everything else, the fast tier included.
list(LENGTH sensitivity_regression_slow_cases n_slow_cases)
set(slow_tier_sentinel
    "sensitivity_slow_tier_opted_in__${n_slow_cases}_cases_else_all_skipped")
add_test(NAME ${slow_tier_sentinel}
    COMMAND bash -c "test -n \"$NEKO_TOP_RUN_SENSITIVITY_REGRESSION\" || \
        (echo 'Set NEKO_TOP_RUN_SENSITIVITY_REGRESSION=1 to opt in.' >&2 && \
        exit 1)"
)
set_tests_properties(${slow_tier_sentinel} PROPERTIES
    LABELS "sensitivity;sensitivity_slow"
    TIMEOUT 10
)

# --------------------------------------------------------------------------- #
