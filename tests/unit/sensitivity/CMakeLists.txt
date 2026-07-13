# Define the CTest for sensitivity
set(CMAKE_FOLDER "Tests/Unit/Sensitivity")

# Provide a list of .case files
set(test_list
    "volume.case"
    "volume_filtered.case"
    "viscous_dissipation.case"
    "brinkman_dissipation.case"
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
    ${CMAKE_CURRENT_SOURCE_DIR}/sensitivity.f90
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
# Loop over all case files and create tests
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
    set_tests_properties(${case_name}_sensitivity
        PROPERTIES
            DEPENDS "problem_tester_sensitivity_bin;sensitivity_test_prepare"
            LABELS "sensitivity;unit"
            ENVIRONMENT "NEKO_LOG_FILE=neko.log"
            TIMEOUT 300
    )
endforeach()

# ============================================================================ #
# Final dependencies
add_dependencies(Unit sensitivity_test_prepare)
