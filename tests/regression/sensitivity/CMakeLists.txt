# Define the CTest for sensitivity (regression suite)
set(CMAKE_FOLDER "Tests/Regression/Sensitivity")

# --------------------------------------------------------------------------- #
# .case files to run in the regression suite (long-running)
# NOTE: 
# These have all been commented out such that they don't run along side
# the unit tests.
# In order to test these one can uncomment the tests, run:
# ./setup.sh -t
# And then within the build directory run:
# python plotter_FD_check.py
# To compare against the reference values
set(test_list
    #"minimum_dissipation.case"
    #"lube_term.case"
    #"min_dis_and_lube.case"
    #"min_dis_and_lube_weights.case"
    # "passive_scalar.case"
)

# Common driver source (same as unit)
set(SENSITIVITY_DRIVER ${CMAKE_CURRENT_SOURCE_DIR}/problem_tester.f90)

# Ensure the neko library module dir is found (from the main target)
get_target_property(NEKO_TOP_MODDIR neko-top Fortran_MODULE_DIRECTORY)

# --------------------------------------------------------------------------- #
# Preparation target: copies helper scripts into the binary dir
add_custom_target(reg_sensitivity_test_prepare ALL
    COMMENT "Preparing regression sensitivity test files"
)

add_custom_command(
    TARGET reg_sensitivity_test_prepare
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            ${CMAKE_CURRENT_SOURCE_DIR}/plotter_FD_check.py
            ${CMAKE_CURRENT_BINARY_DIR}/plotter_FD_check.py
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            ${CMAKE_CURRENT_SOURCE_DIR}/steady_state_plotter.py
            ${CMAKE_CURRENT_BINARY_DIR}/steady_state_plotter.py
)

# --------------------------------------------------------------------------- #
# Regression sensitivity utilities (Fortran modules)
add_library(reg_sensitivity_utils
    ${CMAKE_CURRENT_SOURCE_DIR}/sensitivity.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/user.f90
)

# Use a private .mod folder to avoid clashes with unit targets
set(reg_sens_moddir ${CMAKE_CURRENT_BINARY_DIR}/modules)
set_target_properties(reg_sensitivity_utils
    PROPERTIES Fortran_MODULE_DIRECTORY ${reg_sens_moddir})

# Make both NEKO and this target's .mod directory visible to dependents
target_include_directories(reg_sensitivity_utils
    PUBLIC
        ${NEKO_TOP_MODDIR}
        ${reg_sens_moddir}
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

add_dependencies(reg_sensitivity_test_prepare reg_sensitivity_utils)

# --------------------------------------------------------------------------- #
# Neko runtime environment for prepare.sh (mirrors unit)
set(neko_PATH "${neko_PREFIX}/bin:$ENV{PATH}")
string(REPLACE ";" ":"
    neko_LD_LIBRARY_PATH "${neko_LIBRARY_DIRS}:$ENV{LD_LIBRARY_PATH}"
)

# Copy and execute the prepare script in the binary dir
add_custom_target(
    reg_sensitivity_prepare_script ALL
    COMMENT "Copying and executing regression sensitivity prepare.sh"
    SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/prepare.sh
    BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/prepare.sh
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            ${CMAKE_CURRENT_SOURCE_DIR}/prepare.sh
            ${CMAKE_CURRENT_BINARY_DIR}/prepare.sh
    COMMAND ${CMAKE_COMMAND} -E env
            LD_LIBRARY_PATH="${neko_LD_LIBRARY_PATH}" PATH="${neko_PATH}"
            ./prepare.sh -N8 > ${CMAKE_CURRENT_BINARY_DIR}/prepare.log
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_dependencies(reg_sensitivity_test_prepare reg_sensitivity_prepare_script)

# --------------------------------------------------------------------------- #
# Build the single driver executable for regression
add_executable(reg_problem_tester_sensitivity_bin
    ${SENSITIVITY_DRIVER}
)
target_link_libraries(reg_problem_tester_sensitivity_bin
    reg_sensitivity_utils
)

# Ensure driver builds after the utils (target exists now)
add_dependencies(reg_problem_tester_sensitivity_bin reg_sensitivity_utils)

# The driver must also see the regression .mod directory
target_include_directories(reg_problem_tester_sensitivity_bin
    PRIVATE ${reg_sens_moddir})

# Let the umbrella target trigger both build+prepare steps (if present)
if (TARGET Unit)
  add_dependencies(Unit reg_problem_tester_sensitivity_bin)
  add_dependencies(Unit reg_sensitivity_test_prepare)
endif()

# --------------------------------------------------------------------------- #
# Per-case setup and CTest entries
foreach(case_file ${test_list})
    # Copy matching reference CSV if it exists:
    set(ref_src "${CMAKE_CURRENT_SOURCE_DIR}/ref_FD_check_${case_file}.csv")
    set(ref_dst "${CMAKE_CURRENT_BINARY_DIR}/ref_FD_check_${case_file}.csv")

    if(EXISTS "${ref_src}")
        add_custom_target(
            ${case_name}_ref_fd_reg ALL
            COMMENT "Copying reference CSV for ${case_file}"
            BYPRODUCTS "${ref_dst}"
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                    "${ref_src}" "${ref_dst}"
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        )
        # Ensure the prepare step brings these in before plotting/tests
        add_dependencies(reg_sensitivity_test_prepare ${case_name}_ref_fd_reg)
    else()
        message(STATUS "No reference CSV found for ${case_file} (${ref_src})")
    endif()

    # name without extension
    get_filename_component(case_name ${case_file} NAME_WE)

    # Copy each .case file into the binary dir (namespaced target to avoid clashes)
    add_custom_target(
        ${case_name}_case_file_reg ALL
        COMMENT "Copying case file ${case_file} to the regression sensitivity directory"
        SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/${case_file}
        BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/${case_file}
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                ${CMAKE_CURRENT_SOURCE_DIR}/${case_file}
                ${CMAKE_CURRENT_BINARY_DIR}/${case_file}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    add_dependencies(reg_sensitivity_test_prepare ${case_name}_case_file_reg)

    # Test (namespaced to avoid colliding with unit tests)
    add_test(NAME ${case_name}_sensitivity_reg
        COMMAND $<TARGET_FILE:reg_problem_tester_sensitivity_bin> ${case_file}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )

    # Make sure the test waits for driver+prepare, label as regression
    set_tests_properties(${case_name}_sensitivity_reg
        PROPERTIES
            DEPENDS "reg_problem_tester_sensitivity_bin;reg_sensitivity_test_prepare"
            LABELS "sensitivity;regression"
            ENVIRONMENT "NEKO_LOG_FILE=neko.log"
            TIMEOUT 300000
    )
endforeach()

# --------------------------------------------------------------------------- #
