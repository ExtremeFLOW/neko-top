# CMakeLists.txt for adjoint-operator unit tests.
#
# adjoint_curl_curl_transpose: verifies that adjoint_curl_curl_load is the
# exact transpose of the primal's Karniadakis-Israeli-Orszag curl-curl
# pressure-residual contribution, built independently from Neko's own
# pnpn_prs_res_cpu, for random continuous fields at P = 2 and P = 5.

set(CMAKE_FOLDER "Tests/Unit/Adjoint")

add_pfunit_ctest(adjoint_curl_curl_transpose
    TEST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/adjoint_curl_curl_transpose.pf
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    LINK_LIBRARIES Neko-TOP::neko-top
    MAX_PES 2
)

add_dependencies(Unit
    adjoint_curl_curl_transpose
)

set_tests_properties(adjoint_curl_curl_transpose PROPERTIES
    LABELS "unit;adjoint"
    ENVIRONMENT "NEKO_LOG_FILE=adjoint_curl_curl_transpose.log"
    TIMEOUT 120
)

# Disable treating unused-variable as error for pfunit tests
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(adjoint_curl_curl_transpose PRIVATE
        $<$<COMPILE_LANGUAGE:Fortran>:-Wno-unused-variable>
    )
endif()
