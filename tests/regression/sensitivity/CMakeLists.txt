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
set(sensitivity_regression_cases
    "viscous_dissipation"
    "brinkman_dissipation"
    "dissipation"
    "dissipation_weights"
    "passive_scalar"
    "unsteady_viscous_dissipation"
    "unsteady_brinkman_dissipation"
)

foreach(case_name ${sensitivity_regression_cases})
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
        LABELS "regression;sensitivity"
        RESOURCE_LOCK "sensitivity_regression_workspace"
        SKIP_RETURN_CODE 77
        TIMEOUT ${case_timeout}
    )
endforeach()

# --------------------------------------------------------------------------- #
