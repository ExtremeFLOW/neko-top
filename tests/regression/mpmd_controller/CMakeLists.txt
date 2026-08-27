# Define the CTest for the Neko/Python MPMD controller regression.
set(CMAKE_FOLDER "Tests/Regression/MPMDController")

get_target_property(NEKO_TOP_MOD_DIR neko-top Fortran_MODULE_DIRECTORY)
set(LOCAL_MOD_DIR ${CMAKE_CURRENT_BINARY_DIR}/modules)
set(CMAKE_Fortran_MODULE_DIRECTORY ${LOCAL_MOD_DIR} CACHE PATH
    "Fortran module directory" FORCE)

set(neko_PATH "${neko_PREFIX}/bin:$ENV{PATH}")
string(REPLACE ";" ":"
    neko_LD_LIBRARY_PATH "${neko_LIBRARY_DIRS}:$ENV{LD_LIBRARY_PATH}")

add_executable(reg_mpmd_controller_bin
    ${CMAKE_CURRENT_SOURCE_DIR}/driver.f90
)
target_link_libraries(reg_mpmd_controller_bin PRIVATE Neko-TOP::neko-top)
target_include_directories(reg_mpmd_controller_bin PRIVATE ${LOCAL_MOD_DIR})
set_target_properties(reg_mpmd_controller_bin PROPERTIES
    FOLDER "Tests/Regression/MPMDController"
    RUNTIME_OUTPUT_NAME "reg_mpmd_controller_bin"
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_dependencies(Regression reg_mpmd_controller_bin)

set(CMAKE_Fortran_MODULE_DIRECTORY ${NEKO_TOP_MOD_DIR} CACHE PATH
    "Fortran module directory" FORCE)

find_program(MPMD_PYTHON_EXECUTABLE NAMES python3 python)

if(MPMD_PYTHON_EXECUTABLE)
    # This local-only regression requires mpi4py and is not run by CI.
    add_test(NAME mpmd_controller_build
        COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR}
            --target reg_mpmd_controller_bin
    )
    set_tests_properties(mpmd_controller_build PROPERTIES
        LABELS "regression;mpmd"
        TIMEOUT 900
    )

    add_test(NAME mpmd_controller_regression
        COMMAND ${CMAKE_COMMAND} -E env
            PATH=${neko_PATH}
            LD_LIBRARY_PATH=${neko_LD_LIBRARY_PATH}
            PYTHON_BIN=${MPMD_PYTHON_EXECUTABLE}
            NEKO_LAUNCHER=mpirun
            NEKO_STARTUP_DELAY=0
            bash ${CMAKE_CURRENT_SOURCE_DIR}/run.sh
            ${CMAKE_CURRENT_SOURCE_DIR}/unused.case
            ${CMAKE_CURRENT_SOURCE_DIR}/controller_peer.py
            $<TARGET_FILE:reg_mpmd_controller_bin>
            ${CMAKE_CURRENT_BINARY_DIR}/mpmd_controller.log
    )
    set_tests_properties(mpmd_controller_regression PROPERTIES
        DEPENDS "mpmd_controller_build"
        LABELS "regression;mpmd"
        RESOURCE_LOCK "mpmd_controller_workspace"
        TIMEOUT 120
    )
endif()
