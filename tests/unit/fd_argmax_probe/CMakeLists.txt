# CMakeLists.txt for the fd_argmax_probe_index unit tests
# (tests/shared/sensitivity.f90).
#
# Kept in its own directory rather than alongside tests/unit/fd_criterion:
# that directory's whole premise is "no MPI, no mesh, no simulation" (see
# its own CMakeLists.txt), while this routine's entire contract is about
# cross-rank coordination, so a dedicated MAX_PES 2 target belongs here
# instead.
#
# fd_argmax_probe_test.pf: six tests, mixing npes=[2] (the cross-rank
# tie-break the routine exists for) with one npes=[1] test (the
# single-rank reduction to a plain maxloc) -- see the .pf file's own doc
# comment for how the two coexist in one MAX_PES 2 executable.
# OTHER_SOURCES repeats the fd_criterion.f90 + sensitivity.f90 pair from
# tests/unit/fd_criterion/CMakeLists.txt: sensitivity `use`s fd_criterion,
# so both have to be compiled into this executable too, fd_criterion first
# so its module exists by then.

set(CMAKE_FOLDER "Tests/Unit/FDArgmaxProbe")

add_pfunit_ctest(fd_argmax_probe_test
    TEST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/fd_argmax_probe_test.pf
    OTHER_SOURCES ${SHARED_TEST_SOURCES}/fd_criterion.f90
        ${SHARED_TEST_SOURCES}/sensitivity.f90
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    LINK_LIBRARIES Neko-TOP::neko-top
    MAX_PES 2
)

add_dependencies(Unit fd_argmax_probe_test)

set_tests_properties(fd_argmax_probe_test PROPERTIES
    LABELS "unit;sensitivity"
    ENVIRONMENT "NEKO_LOG_FILE=fd_argmax_probe_test.log"
    TIMEOUT 120
)

# Disable treating unused-variable as error for pfunit tests
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(fd_argmax_probe_test PRIVATE
        $<$<COMPILE_LANGUAGE:Fortran>:-Wno-unused-variable>
    )
endif()
