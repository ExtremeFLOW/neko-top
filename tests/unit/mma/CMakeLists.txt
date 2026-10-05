# ============================================================================ #
# mma-backend_consistency
# ============================================================================ #
# Registered here, like tests/unit/optimizer/ and tests/unit/checkpointing/,
# rather than through the generic PFUNIT_TEST_SOURCES loop in
# tests/unit/CMakeLists.txt, because it needs settings that loop does not
# give every test:
#
#   - MAX_PES 2 unconditionally. The generic loop's MAX_PES expression
#     ($<IF:$<STREQUAL:${DEVICE_TYPE},CPU>,2,1>) gives only 1 rank on
#     device builds, which would fail every npes=2 pFUnit variant in this
#     file ("Insufficient processes to run this test."). Neko's
#     device_init() only requires device_count() == 1 per rank (see
#     src/device/device.F90), not a dedicated device per rank, so 2
#     ranks sharing one GPU is supported.
#   - A longer TIMEOUT: the test drives a cpu-backend and a device-backend
#     mma_t through several outer MMA iterations each (see the test
#     file's own "Runtime" docstring for measured figures), so it runs
#     well past a typical unit test's TIMEOUT even when the GPU is
#     otherwise idle, and a contended or thermally throttled GPU can
#     slow it further.
#   - RUN_SERIAL TRUE: `ctest -j` is otherwise free to run this test at
#     the same time as another GPU-bound test on the same device. That
#     changes wall time, not numerics (both backends' computation here
#     is deterministic regardless of what else shares the GPU, so a
#     contended run is slower, not wrong) -- but the slowdown is exactly
#     what the long TIMEOUT above has to budget for. A RESOURCE_LOCK
#     naming a resource no other test locks would not achieve even that:
#     CTest only serialises tests that share a lock name, so it would
#     not stop `ctest -j` from running this test alongside any other GPU
#     test. RUN_SERIAL instead tells CTest never to run this test
#     concurrently with any other test, which does.

add_pfunit_ctest(mma-backend_consistency
    TEST_SOURCES backend_consistency.pf
    LINK_LIBRARIES Neko-TOP::neko-top
    MAX_PES 2
)

set(mma_backend_consistency_log
    "${CMAKE_BINARY_DIR}/tests/unit/mma-backend_consistency.log")

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(mma-backend_consistency PRIVATE
        $<$<COMPILE_LANGUAGE:Fortran>:-Wno-unused-variable>
    )
endif()

set_target_properties(mma-backend_consistency PROPERTIES
    FOLDER "Tests/Unit"
)
set_tests_properties(mma-backend_consistency PROPERTIES
    LABELS "unit;mma"
    ENVIRONMENT "NEKO_LOG_FILE=${mma_backend_consistency_log}"
    TIMEOUT 600
    RUN_SERIAL TRUE
)
add_dependencies(Unit mma-backend_consistency)
