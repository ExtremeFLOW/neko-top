# CMakeLists.txt for the mapping-function unit tests.
#
# Finite-difference verification of the RAMP and Borrvall & Peterson mappings:
# each forward mapping is central-differenced and compared against the
# corresponding backward (chain-rule) mapping.

set(CMAKE_FOLDER "Tests/Unit/Mapping")

# ============================================================================ #
# Define the test itself

add_pfunit_ctest(mapping_finite_difference
    TEST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/mapping_finite_difference.pf
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    LINK_LIBRARIES Neko-TOP::neko-top
    MAX_PES 1
)

# Set test properties
set_tests_properties(mapping_finite_difference PROPERTIES
    LABELS "unit;mapping"
    ENVIRONMENT "NEKO_LOG_FILE=mapping_finite_difference.log"
    TIMEOUT 120
)

# Disable treating unused-variable as error for pfunit tests
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(mapping_finite_difference PRIVATE
        $<$<COMPILE_LANGUAGE:Fortran>:-Wno-unused-variable>
    )
endif()
