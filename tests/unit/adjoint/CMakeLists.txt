# CMakeLists.txt for the adjoint unit tests.
#
# adv_adjoint_no_dealias_transpose: the non-dealiased adjoint advection terms
# (adv_lin_no_dealias_t) are verified to be the exact discrete transposes of
# the forward (linearised) advection terms, for the fluid and for a passive
# scalar, via a Euclidean dot-product identity on a two-element hex mesh.
#
# adjoint_curl_curl_transpose: verifies that adjoint_curl_curl_load is the
# exact transpose of the primal's Karniadakis-Israeli-Orszag curl-curl
# pressure-residual contribution, built independently from Neko's own
# pnpn_prs_res_cpu, for random continuous fields at P = 2 and P = 5.

set(CMAKE_FOLDER "Tests/Unit/Adjoint")

# ============================================================================ #
# Define the tests themselves

add_pfunit_ctest(adv_adjoint_no_dealias_transpose
    TEST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/adv_adjoint_no_dealias_transpose.pf
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    LINK_LIBRARIES Neko-TOP::neko-top
    MAX_PES 1
)

add_pfunit_ctest(adjoint_curl_curl_transpose
    TEST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/adjoint_curl_curl_transpose.pf
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    LINK_LIBRARIES Neko-TOP::neko-top
    MAX_PES 2
)

add_dependencies(Unit
    adv_adjoint_no_dealias_transpose
    adjoint_curl_curl_transpose
)

# Set test properties
set_tests_properties(adv_adjoint_no_dealias_transpose PROPERTIES
    LABELS "unit;adjoint"
    ENVIRONMENT "NEKO_LOG_FILE=adv_adjoint_no_dealias_transpose.log"
    TIMEOUT 120
)

set_tests_properties(adjoint_curl_curl_transpose PROPERTIES
    LABELS "unit;adjoint"
    ENVIRONMENT "NEKO_LOG_FILE=adjoint_curl_curl_transpose.log"
    TIMEOUT 120
)

# Disable treating unused-variable as error for pfunit tests
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(adv_adjoint_no_dealias_transpose PRIVATE
        $<$<COMPILE_LANGUAGE:Fortran>:-Wno-unused-variable>
    )
    target_compile_options(adjoint_curl_curl_transpose PRIVATE
        $<$<COMPILE_LANGUAGE:Fortran>:-Wno-unused-variable>
    )
endif()
