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

    # Shared sensitivity modules. fd_criterion is listed first: sensitivity
    # uses it, and the Fortran module it defines has to exist by then.
    ${SHARED_TEST_SOURCES}/fd_criterion.f90
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
# Register the (opt-in, non-blocking) regression sensitivity checks.
#
# These run realistic ~1000-timestep forward+adjoint solves per perturbation
# and are not part of the default/PR-blocking test budget — they only run
# when NEKO_TOP_RUN_SENSITIVITY_REGRESSION is set (see
# run_sensitivity_regression_check.sh, which exits 77 — CTest's "skipped"
# code — otherwise). This mirrors the opt-in pattern already used for
# tests/regression/valgrind (there gated on `valgrind` being installed;
# here gated on an explicit opt-in since the tool is always installed).

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

# passive_scalar's case runs 10x longer than the others (end_time 10.0 vs
# 1.0 at the same timestep) — give it a much larger TIMEOUT rather than a
# shared one sized for the cheaper cases.
#
# two_volume_constraints and objective_and_constraint are the multi-target
# cases: the first checks two constraints in one run, which is what the
# (n_constraints, n_design) sensitivity layout has to be right for, and the
# second checks an objective and a constraint together, which the old
# single-target dispatch refused outright. Both are copies of the cheap
# unsteady case (five timesteps), so they cost seconds rather than the
# minutes the steady cases take.
#
# scalar_mixing needs its own mesh and its own resolution (see
# run_sensitivity_regression_check.sh): it is the only non-periodic-y case
# here, since it splits a real y in [0,1] outlet into quarters. Its
# end_time (5.0) and Pe (1000) were re-measured from scratch after the Neko
# inlet fix c528651e9 -- before it, `expression_dirichlet` never set
# bc_type and this case had no scalar inflow at all, so every parameter
# derived against it, including the 185.0 end_time this comment used to
# justify, was void. Its dt is the same 1.0e-3 the rest of the suite uses;
# the smaller dt an older comment described was suppressing an instability
# whose real cause was the extrapolation order. See the case file's own
# header comment for all of those derivations. It remains the most
# expensive case in this suite, but by a factor of a few rather than of
# hundreds.
set(sensitivity_regression_cases
    "viscous_dissipation"
    "brinkman_dissipation"
    "dissipation"
    "dissipation_weights"
    "passive_scalar"
    "unsteady_viscous_dissipation"
    "unsteady_brinkman_dissipation"
    "two_volume_constraints"
    "objective_and_constraint"
    "scalar_mixing"
)

foreach(case_name ${sensitivity_regression_cases})
    # scalar_mixing's mesh (26x13 elements, see
    # run_sensitivity_regression_check.sh) is ~2.7x the element count of the
    # 16x8 shared box the other cases use, and its end_time is 5x the
    # steady siblings' -- the most expensive case in this suite. A higher
    # MPI rank count
    # was tried here to keep its wall time tractable (8 measured ~2.75x
    # faster than 2 on an isolated benchmark) and REJECTED for a measured
    # reason: this container's /dev/shm is capped at 64 MiB (see the
    # neko-container-environment skill), and OpenMPI's vader BTL requests a
    # 16 MiB sm_segment per rank -- 8 ranks demands 128 MiB, which cannot
    # fit regardless of stale-segment cleanup. Measured directly: an 8-rank
    # run of this exact case failed after 123.6s with repeated "not enough
    # space for .../sm_segment..." warnings and OSC UCX errors, not a clean
    # FD-tolerance verdict. Stays at the suite default (2 ranks, 32 MiB,
    # comfortable headroom) rather than any rank count above 4.
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
    elseif(case_name STREQUAL "scalar_mixing")
        # MEASURED, not guessed: CTest timed the complete four-point sweep
        # -- baseline, adjoint, 4 perturbed solves, 1 restore -- at 2336 s
        # and 2372 s on two consecutive runs at 2 ranks on an otherwise
        # idle box, i.e. ~40 min. 7200 s is ~3x that, leaving room for
        # machine load without letting a hung run sit for hours. The 48 h
        # that used to sit here was sized for the pre-inlet-fix end_time
        # of 185, which is gone.
        set(case_timeout 7200)
    else()
        set(case_timeout 1800)
    endif()
    # Two distinct reasons this lane can decline to run, and CTest's
    # SKIP_RETURN_CODE holds exactly one code (a "77;78" list matches only
    # the first — measured with cmake 4.2.3), so the second reason is
    # matched on the driver's banner instead:
    #   77 — the opt-in gate above (run_sensitivity_regression_check.sh);
    #   the banner — a single-precision build, where the driver exits with
    #   FD_SKIP_EXIT_CODE = 78 having made no assertion. Both must report
    #   SKIP rather than a pass that checked nothing.
    set_tests_properties(${case_name}_sensitivity_regression PROPERTIES
        DEPENDS "sensitivity_regression_build"
        LABELS "regression;sensitivity"
        RESOURCE_LOCK "sensitivity_regression_workspace"
        SKIP_RETURN_CODE 77
        SKIP_REGULAR_EXPRESSION "FD sweep: SKIPPED"
        TIMEOUT ${case_timeout}
    )
endforeach()

# --------------------------------------------------------------------------- #
