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
    ${CMAKE_CURRENT_SOURCE_DIR}/sensitivity.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/user.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/problem_tester.f90
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
