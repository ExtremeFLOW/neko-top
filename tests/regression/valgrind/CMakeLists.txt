# Define the CTest for valgrind regression checks.
set(CMAKE_FOLDER "Tests/Regression/Valgrind")

# Keep module outputs isolated for this target.
get_target_property(NEKO_TOP_MOD_DIR neko-top Fortran_MODULE_DIRECTORY)
set(LOCAL_MOD_DIR ${CMAKE_CURRENT_BINARY_DIR}/modules)
set(CMAKE_Fortran_MODULE_DIRECTORY ${LOCAL_MOD_DIR} CACHE PATH
    "Fortran module directory" FORCE)

# Runtime environment for Neko tools and libraries.
set(neko_PATH "${neko_PREFIX}/bin:$ENV{PATH}")
string(REPLACE ";" ":"
    neko_LD_LIBRARY_PATH "${neko_LIBRARY_DIRS}:$ENV{LD_LIBRARY_PATH}")

# Build a dedicated executable for this regression test.
add_executable(reg_valgrind_bin
    ${Neko-TOP_SOURCE_DIR}/sources/drivers/topopt-user.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/user.f90
)
target_link_libraries(reg_valgrind_bin PRIVATE Neko-TOP::neko-top)
target_include_directories(reg_valgrind_bin PRIVATE ${LOCAL_MOD_DIR})
set_target_properties(reg_valgrind_bin PROPERTIES
    FOLDER "Tests/Regression/Valgrind"
    RUNTIME_OUTPUT_NAME "reg_valgrind_bin"
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

# Keep umbrella regression target aware of this executable.
add_dependencies(Regression reg_valgrind_bin)

# Restore module output setting after configuring this target.
set(CMAKE_Fortran_MODULE_DIRECTORY ${NEKO_TOP_MOD_DIR} CACHE PATH
    "Fortran module directory" FORCE)

# Register valgrind regression test.
add_test(NAME unsteady_mixer_valgrind_build
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target reg_valgrind_bin
)
set_tests_properties(unsteady_mixer_valgrind_build PROPERTIES
    LABELS "regression;valgrind"
    TIMEOUT 900
)

add_test(NAME unsteady_mixer_valgrind_regression
    COMMAND ${CMAKE_COMMAND} -E env
        PATH=${neko_PATH}
        LD_LIBRARY_PATH=${neko_LD_LIBRARY_PATH}
        NEKO_DIR=${neko_PREFIX}
        NEKO_LOG_FILE=neko_valgrind.log
        NEKO_TOP_LOG_FILE=neko_valgrind_nekotop.log
        OMP_NUM_THREADS=8
        OPENBLAS_NUM_THREADS=1
        VALGRIND_MESH_NX=8
        VALGRIND_MESH_NY=2
        VALGRIND_MESH_NZ=2
        bash ${CMAKE_CURRENT_SOURCE_DIR}/run_valgrind_check.sh
        $<TARGET_FILE:reg_valgrind_bin>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/steady_mixer_valgrind_cpu.case
        ${CMAKE_CURRENT_SOURCE_DIR}/valgrind_runtime.supp
        ${CMAKE_CURRENT_SOURCE_DIR}/prepare.sh
        ${CMAKE_CURRENT_SOURCE_DIR}/check_valgrind.py
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(unsteady_mixer_valgrind_regression PROPERTIES
    DEPENDS "unsteady_mixer_valgrind_build"
    LABELS "regression;valgrind"
    RESOURCE_LOCK "valgrind_regression_workspace"
    SKIP_RETURN_CODE 77
    TIMEOUT 3600
)

add_test(NAME unsteady_mixer_valgrind_unsteady_regression
    COMMAND ${CMAKE_COMMAND} -E env
        PATH=${neko_PATH}
        LD_LIBRARY_PATH=${neko_LD_LIBRARY_PATH}
        NEKO_DIR=${neko_PREFIX}
        NEKO_LOG_FILE=neko_valgrind.log
        NEKO_TOP_LOG_FILE=neko_valgrind_nekotop.log
        OMP_NUM_THREADS=8
        OPENBLAS_NUM_THREADS=1
        VALGRIND_MESH_NX=8
        VALGRIND_MESH_NY=2
        VALGRIND_MESH_NZ=2
        bash ${CMAKE_CURRENT_SOURCE_DIR}/run_valgrind_check.sh
        $<TARGET_FILE:reg_valgrind_bin>
        ${CMAKE_CURRENT_SOURCE_DIR}/cases/unsteady_mixer_valgrind_cpu.case
        ${CMAKE_CURRENT_SOURCE_DIR}/valgrind_runtime.supp
        ${CMAKE_CURRENT_SOURCE_DIR}/prepare.sh
        ${CMAKE_CURRENT_SOURCE_DIR}/check_valgrind.py
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(unsteady_mixer_valgrind_unsteady_regression PROPERTIES
    DEPENDS "unsteady_mixer_valgrind_build"
    LABELS "regression;valgrind"
    RESOURCE_LOCK "valgrind_regression_workspace"
    SKIP_RETURN_CODE 77
    TIMEOUT 3600
)
