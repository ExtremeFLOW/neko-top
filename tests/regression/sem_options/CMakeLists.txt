# Define the CTest for SEM option regression runs.
set(CMAKE_FOLDER "Tests/Regression/SEMOptions")

# --------------------------------------------------------------------------- #
# Build a standard topopt driver for the regression case files.
add_executable(reg_sem_options_bin
    ${Neko-TOP_SOURCE_DIR}/sources/drivers/topopt.f90
)
target_link_libraries(reg_sem_options_bin PRIVATE Neko-TOP::neko-top)

set_target_properties(reg_sem_options_bin
    PROPERTIES
        FOLDER "Tests/Regression/SEMOptions"
        RUNTIME_OUTPUT_NAME "reg_sem_options_bin"
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

# Let the umbrella target trigger this build.
add_dependencies(Regression reg_sem_options_bin)

# --------------------------------------------------------------------------- #
