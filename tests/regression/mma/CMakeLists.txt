# Define the CTest for sensitivity (regression suite)
set(CMAKE_FOLDER "Tests/Regression/MMA")

# --------------------------------------------------------------------------- #
# Ensure the neko library module dir is found (from the main target)
get_target_property(NEKO_TOP_MOD_DIR neko-top Fortran_MODULE_DIRECTORY)

# --------------------------------------------------------------------------- #
# Regression sensitivity utilities (Fortran modules)

# Use a private .mod folder to avoid clashes with unit targets
set(LOCAL_MOD_DIR ${CMAKE_CURRENT_BINARY_DIR}/modules)
set(CMAKE_Fortran_MODULE_DIRECTORY ${LOCAL_MOD_DIR} CACHE PATH
    "Fortran module directory" FORCE)

# --------------------------------------------------------------------------- #
# Build the single driver executable for regression

add_executable(reg_mma_bin
    ${CMAKE_CURRENT_SOURCE_DIR}/driver.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/example_problem.f90
)
target_link_libraries(reg_mma_bin PRIVATE Neko-TOP::neko-top)

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

set(CMAKE_Fortran_MODULE_DIRECTORY ${NEKO_TOP_MOD_DIR} CACHE PATH
    "Fortran module directory" FORCE)

# --------------------------------------------------------------------------- #
