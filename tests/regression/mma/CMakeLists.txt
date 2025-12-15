# Define the CTest for sensitivity (regression suite)
set(CMAKE_FOLDER "Tests/Regression/MMA")

# --------------------------------------------------------------------------- #
# Ensure the neko library module dir is found (from the main target)
get_target_property(NEKO_TOP_MOD_DIR neko-top Fortran_MODULE_DIRECTORY)

# --------------------------------------------------------------------------- #
# Regression sensitivity utilities (Fortran modules)

# Use a private .mod folder to avoid clashes with unit targets
set(LOCAL_MOD_DIR ${CMAKE_CURRENT_BINARY_DIR}/modules)

# --------------------------------------------------------------------------- #
# Build the single driver executable for regression

add_executable(reg_mma_bin
    ${CMAKE_CURRENT_SOURCE_DIR}/driver.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/example_problem.f90
)
target_link_libraries(reg_mma_bin PRIVATE Neko-TOP::neko-top)

# Additional regression: MMA HDF5 read/write test
add_executable(reg_mma_hdf5_bin
    ${CMAKE_CURRENT_SOURCE_DIR}/test_mma_hdf5.f90
)
target_link_libraries(reg_mma_hdf5_bin PRIVATE Neko-TOP::neko-top)

set_target_properties(reg_mma_hdf5_bin
    PROPERTIES
        FOLDER "Tests/Regression/MMA"
        RUNTIME_OUTPUT_NAME "reg_mma_hdf5_bin"
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
# Expose local .mod dir to the test target
target_include_directories(reg_mma_hdf5_bin PRIVATE ${LOCAL_MOD_DIR})

# Register as a CTest
add_test(NAME reg_mma_hdf5_test
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/reg_mma_hdf5_bin
)

# Ensure CTest builds the test target before running and run in source dir
set_tests_properties(reg_mma_hdf5_test PROPERTIES
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS reg_mma_hdf5_bin
    TIMEOUT 30
)

add_dependencies(Regression reg_mma_hdf5_bin)

# The driver must also see the regression .mod directory
target_include_directories(reg_mma_bin
    PRIVATE ${LOCAL_MOD_DIR})

set_target_properties(reg_mma_bin
    PROPERTIES
        FOLDER "Tests/Regression/MMA"
        RUNTIME_OUTPUT_NAME "reg_mma_bin"
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

# Let the umbrella target trigger both build+prepare steps (if present)
add_dependencies(Regression reg_mma_bin)

# --------------------------------------------------------------------------- #
