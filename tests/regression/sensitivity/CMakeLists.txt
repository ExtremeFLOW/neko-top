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
# Regression sensitivity utilities (Fortran modules)
add_library(reg_sensitivity_utils
    ${CMAKE_CURRENT_SOURCE_DIR}/sensitivity.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/user.f90
)

# Use a private .mod folder to avoid clashes with unit targets
set(LOCAL_MOD_DIR ${CMAKE_CURRENT_BINARY_DIR}/modules)
set_target_properties(reg_sensitivity_utils
    PROPERTIES Fortran_MODULE_DIRECTORY ${LOCAL_MOD_DIR})

# Make both NEKO and this target's .mod directory visible to dependents
target_include_directories(reg_sensitivity_utils
    PUBLIC
        ${NEKO_TOP_MOD_DIR}
        ${LOCAL_MOD_DIR}
)

target_link_libraries(reg_sensitivity_utils
    PUBLIC
        Neko-TOP::neko-top
        PkgConfig::neko
        PkgConfig::json-fortran
        MPI::MPI_Fortran
        $<$<BOOL:${BLAS_FOUND}>:BLAS::BLAS>
        $<$<BOOL:${LAPACK_FOUND}>:LAPACK::LAPACK>
)

# --------------------------------------------------------------------------- #
# Build the single driver executable for regression
add_executable(reg_problem_tester_sensitivity_bin
    ${CMAKE_CURRENT_SOURCE_DIR}/problem_tester.f90
)
target_link_libraries(reg_problem_tester_sensitivity_bin
    reg_sensitivity_utils
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

# --------------------------------------------------------------------------- #
