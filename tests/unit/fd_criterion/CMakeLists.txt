# CMakeLists.txt for the strict finite-difference criterion's own unit
# tests (tests/shared/fd_criterion.f90), kept separate from
# tests/unit/sensitivity (which exercises the harness end-to-end through
# real case files) because this module is deliberately MPI-free and I/O-free
# and is tested directly, with no mesh, no MPI, and no simulation.
#
# fd_criterion_test.pf: pFUnit tests against fd_evaluate's public verdict --
# the five real recorded sweeps, and every guard that returns a verdict
# rather than aborting. Named fd_criterion_test, not fd_criterion, because
# the pFUnit module name must equal the .pf basename, and this file itself
# `use`s the production module `fd_criterion` -- naming both the same would
# make the test module use itself.
#
# fd_criterion_reject: a small standalone probe (NOT a pFUnit test) for the
# four guards that instead reject via `neko_error` (an `error stop`, so
# uncatchable in-process by pFUnit -- see the .pf file's header). Each
# scenario is wired up as its own ctest entry with WILL_FAIL TRUE: the
# guard being present makes the probe abort (non-zero exit), which is the
# PASSING outcome for a WILL_FAIL test.

set(CMAKE_FOLDER "Tests/Unit/FDCriterion")

# ============================================================================ #
# The pFUnit tests.

add_pfunit_ctest(fd_criterion_test
    TEST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/fd_criterion_test.pf
    OTHER_SOURCES ${SHARED_TEST_SOURCES}/fd_criterion.f90
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    LINK_LIBRARIES Neko-TOP::neko-top
    MAX_PES 1
)

add_dependencies(Unit fd_criterion_test)

set_tests_properties(fd_criterion_test PROPERTIES
    LABELS "unit;fd_criterion"
    ENVIRONMENT "NEKO_LOG_FILE=fd_criterion_test.log"
    TIMEOUT 120
)

# Disable treating unused-variable as error for pfunit tests
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(fd_criterion_test PRIVATE
        $<$<COMPILE_LANGUAGE:Fortran>:-Wno-unused-variable>
    )
endif()

# ============================================================================ #
# The four "must abort" guards, probed outside pFUnit.

add_executable(fd_criterion_reject EXCLUDE_FROM_ALL
    ${SHARED_TEST_SOURCES}/fd_criterion.f90
    ${SHARED_TEST_SOURCES}/sensitivity.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/fd_criterion_reject.f90
)
target_link_libraries(fd_criterion_reject PRIVATE Neko-TOP::neko-top)
set_target_properties(fd_criterion_reject PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules
    FOLDER "Tests/Unit/FDCriterion"
)
add_dependencies(Unit fd_criterion_reject)

# Two parallel lists, not a list of pairs: a ";" embedded in a quoted list
# element does not survive an unquoted foreach() expansion -- it flattens
# the whole thing back into single elements, so list(GET ... 1) fails.
set(fd_reject_ids "1" "2" "3" "4")
set(fd_reject_names
    fd_criterion_reject_duplicate_magnitude
    fd_criterion_reject_nan
    fd_criterion_reject_mixed_sign
    fd_criterion_reject_order_tolerance
)

list(LENGTH fd_reject_ids n_fd_reject_scenarios)
math(EXPR fd_reject_last "${n_fd_reject_scenarios} - 1")

foreach(i RANGE ${fd_reject_last})
    list(GET fd_reject_ids ${i} scenario_id)
    list(GET fd_reject_names ${i} scenario_name)

    add_test(NAME ${scenario_name}
        COMMAND $<TARGET_FILE:fd_criterion_reject> ${scenario_id}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    set_tests_properties(${scenario_name} PROPERTIES
        LABELS "unit;fd_criterion"
        WILL_FAIL TRUE
        TIMEOUT 60
    )
endforeach()
