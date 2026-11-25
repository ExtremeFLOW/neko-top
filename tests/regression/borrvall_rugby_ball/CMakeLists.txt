# Borrvall & Petersson (2003) rugby-ball tests.
#
# Two tiers share one driver and one run script:
#
# - borrvall_rugby_ball_regression (labels "regression;rugby_ball"):
#   the gamma = 0.8 physics on a 4x4 mesh at polynomial order 3 for 23 MMA
#   iterations, with the RAMP continuation shortened so that the
#   q = 100 -> 10 switch is exercised (regression_0p8.case). The whole
#   optimization history is compared against reference_data/ by
#   compare_history.py. About a minute on two ranks.
#
# - borrvall_rugby_ball_validation_0p8 (labels "validation;rugby_ball"):
#   the committed borrvall_rugby_0p8.case on the 10x10 mesh, run to
#   completion and checked against the paper (Table III). It takes about
#   80 minutes on eight ranks, so it is only registered when the
#   NEKO_TOP_VALIDATION_TESTS option is ON, and it does not carry the
#   "regression" label, so `ctest -L regression` never selects it.
set(CMAKE_FOLDER "Tests/Regression/BorrvallRugbyBall")

# Keep module outputs isolated for this target.
get_target_property(NEKO_TOP_MOD_DIR neko-top Fortran_MODULE_DIRECTORY)
set(LOCAL_MOD_DIR ${CMAKE_CURRENT_BINARY_DIR}/modules)
set(CMAKE_Fortran_MODULE_DIRECTORY ${LOCAL_MOD_DIR} CACHE PATH
    "Fortran module directory" FORCE)

# Runtime environment for Neko tools (genmeshbox) and libraries.
set(neko_PATH "${neko_PREFIX}/bin:$ENV{PATH}")
string(REPLACE ";" ":"
    neko_LD_LIBRARY_PATH "${neko_LIBRARY_DIRS}:$ENV{LD_LIBRARY_PATH}")

# The plain topology-optimisation driver, as used by examples/rugby_ball.
add_executable(reg_rugby_ball_bin
    ${Neko-TOP_SOURCE_DIR}/sources/drivers/topopt.f90
)
target_link_libraries(reg_rugby_ball_bin PRIVATE Neko-TOP::neko-top)
target_include_directories(reg_rugby_ball_bin PRIVATE ${LOCAL_MOD_DIR})
set_target_properties(reg_rugby_ball_bin PROPERTIES
    FOLDER "Tests/Regression/BorrvallRugbyBall"
    RUNTIME_OUTPUT_NAME "reg_rugby_ball_bin"
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_dependencies(Regression reg_rugby_ball_bin)

# Restore module output setting after configuring this target.
set(CMAKE_Fortran_MODULE_DIRECTORY ${NEKO_TOP_MOD_DIR} CACHE PATH
    "Fortran module directory" FORCE)

add_test(NAME borrvall_rugby_ball_build
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR}
        --target reg_rugby_ball_bin
)
set_tests_properties(borrvall_rugby_ball_build PROPERTIES
    LABELS "regression;rugby_ball"
    TIMEOUT 900
)

# --------------------------------------------------------------------------- #
# Fast regression tier.

set(rugby_regression_dir ${CMAKE_CURRENT_BINARY_DIR}/regression_0p8)
file(MAKE_DIRECTORY ${rugby_regression_dir})

add_test(NAME borrvall_rugby_ball_regression
    COMMAND ${CMAKE_COMMAND} -E env
        PATH=${neko_PATH}
        LD_LIBRARY_PATH=${neko_LD_LIBRARY_PATH}
        NEKO_DIR=${neko_PREFIX}
        bash ${CMAKE_CURRENT_SOURCE_DIR}/run_rugby_ball_check.sh
        $<TARGET_FILE:reg_rugby_ball_bin>
        ${CMAKE_CURRENT_SOURCE_DIR}/regression_0p8.case
        4 2
        python3 ${CMAKE_CURRENT_SOURCE_DIR}/compare_history.py
        --reference
        ${CMAKE_CURRENT_SOURCE_DIR}/reference_data/regression_0p8.csv
    WORKING_DIRECTORY ${rugby_regression_dir}
)
set_tests_properties(borrvall_rugby_ball_regression PROPERTIES
    DEPENDS "borrvall_rugby_ball_build"
    LABELS "regression;rugby_ball"
    TIMEOUT 1800
)

# --------------------------------------------------------------------------- #
# Validation tier (opt-in).

if (NEKO_TOP_VALIDATION_TESTS)
    set(rugby_validation_dir ${CMAKE_CURRENT_BINARY_DIR}/validation_0p8)
    file(MAKE_DIRECTORY ${rugby_validation_dir})

    # Paper, Table III: gamma = 0.8 gives J = 31.75. The case runs a fixed
    # 100 MMA iterations, with q = 10 from iteration 51.
    add_test(NAME borrvall_rugby_ball_validation_0p8
        COMMAND ${CMAKE_COMMAND} -E env
            PATH=${neko_PATH}
            LD_LIBRARY_PATH=${neko_LD_LIBRARY_PATH}
            NEKO_DIR=${neko_PREFIX}
            bash ${CMAKE_CURRENT_SOURCE_DIR}/run_rugby_ball_check.sh
            $<TARGET_FILE:reg_rugby_ball_bin>
            ${CMAKE_CURRENT_SOURCE_DIR}/borrvall_rugby_0p8.case
            10 8
            python3 ${CMAKE_CURRENT_SOURCE_DIR}/check_paper.py
            --objective 31.75 --objective-rtol 0.02
            --iterations 100 --final-q 10.0
        WORKING_DIRECTORY ${rugby_validation_dir}
    )
    set_tests_properties(borrvall_rugby_ball_validation_0p8 PROPERTIES
        DEPENDS "borrvall_rugby_ball_build"
        LABELS "validation;rugby_ball"
        TIMEOUT 14400
    )
endif()

# --------------------------------------------------------------------------- #
