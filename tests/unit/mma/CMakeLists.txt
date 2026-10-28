# CMakeLists.txt for the MMA unit tests.
#
# mma_kkt_nan: a NaN sensitivity on one rank must make the max-norm of the CPU
# KKT residual NaN on every rank, for both subsolvers, so that the optimiser
# cannot read a NaN design as converged. It needs two ranks on every backend,
# so it is registered here rather than in the shared list in the parent
# CMakeLists.txt, which runs device builds with a single rank.
#
# mma_checkpoint_cpu_backend: saving an HDF5 checkpoint of a cpu-backend MMA
# must not change its state, which on a device build lives only on the host.
# It needs HDF5 and two ranks on every backend, for the same reason.

set(CMAKE_FOLDER "Tests/Unit/MMA")

# ============================================================================ #
# Define the tests themselves

add_pfunit_ctest(mma_kkt_nan
    TEST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/mma_kkt_nan.pf
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    LINK_LIBRARIES Neko-TOP::neko-top
    MAX_PES 2
)

add_dependencies(Unit
    mma_kkt_nan
)

# Set test properties
set_tests_properties(mma_kkt_nan PROPERTIES
    LABELS "unit;mma"
    ENVIRONMENT "NEKO_LOG_FILE=mma_kkt_nan.log"
    TIMEOUT 120
)

# Disable treating unused-variable as error for pfunit tests
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(mma_kkt_nan PRIVATE
        $<$<COMPILE_LANGUAGE:Fortran>:-Wno-unused-variable>
    )
endif()

if (HAVE_HDF5)
    add_pfunit_ctest(mma_checkpoint_cpu_backend
        TEST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/mma_checkpoint_cpu_backend.pf
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        LINK_LIBRARIES Neko-TOP::neko-top
        MAX_PES 2
    )

    add_dependencies(Unit
        mma_checkpoint_cpu_backend
    )

    set_tests_properties(mma_checkpoint_cpu_backend PROPERTIES
        LABELS "unit;mma"
        ENVIRONMENT "NEKO_LOG_FILE=mma_checkpoint_cpu_backend.log"
        TIMEOUT 120
    )

    if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
        target_compile_options(mma_checkpoint_cpu_backend PRIVATE
            $<$<COMPILE_LANGUAGE:Fortran>:-Wno-unused-variable>
        )
    endif()
endif()
