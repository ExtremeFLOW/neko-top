# ============================================================================ #
# Neko-TOP tests
# ============================================================================ #
# This CMakeLists.txt file is used to define the tests for the Neko-TOP
# library. The tests are defined using the pFUnit testing framework, which
# is a Fortran testing framework that is also used by Neko.
#
# ============================================================================ #
# Setup the Neko-TOP tests
# ============================================================================ #
# The unit tests should be self contained and not depend on any external
# libraries or files outside of the Neko-TOP library. Each test is defined in a
# ".pf" file, which contains the test code and any necessary setup. The tests
# are compiled and linked against the Neko-TOP library and the MPI library.
#
# New tests should be added to the `PFUNIT_TEST_SOURCES` list below. Each
# test is then added to the CMake project using the `add_pfunit_ctest`
# command, which creates a CTest target for the test. The tests are organized
# into a "Tests/Unit" folder in the CMake project for better organization.
#
# In addition to the simple unit tests using PFUnit, this file also includes
# the subdirectory for sensitivity tests and other more complex tests.

if (BUILD_TESTING)
    # Ensure that the BUILD_TESTING option is enabled.
    message(STATUS "Building Neko-TOP unit tests.")

else()
    return()
endif()

add_custom_target(Unit COMMENT "Build all unit tests.")
set_target_properties(Unit PROPERTIES
    FOLDER "Tests/Unit"
    EXCLUDE_FROM_ALL $<IF:$<BOOL:BUILD_TESTING>,FALSE,TRUE>
)
add_dependencies(Tests Unit)

# Create a file list of the pFUnit test sources. These should be relative paths.
set(PFUNIT_TEST_SOURCES
    # Add the .pf unit test files here.
)

# Add the more complex tests that require additional setup or dependencies.
add_subdirectory(sensitivity)

# ============================================================================ #
# Add the listed tests to the CMake project.

if (NOT PFUNIT_TEST_SOURCES)
    return()
endif()


# Find the pFUnit package. This is required to build the unit tests.
find_package(PFUNIT REQUIRED HINTS ${PFUNIT_DIR})

foreach(pf_file ${PFUNIT_TEST_SOURCES})
    # Ensure the path is relative to the tests/unit directory.
    get_filename_component(pf_file_abs "${pf_file}"
        ABSOLUTE BASE_DIR ${CMAKE_SOURCE_DIR}/tests/unit)

    # Check if the file exists before proceeding.
    if(NOT EXISTS ${pf_file_abs})
        message(FATAL_ERROR "Test file ${pf_file} does not exist.")
    endif()

    # Get the name and relative path of the test file.
    get_filename_component(pf_file_path ${pf_file} DIRECTORY)
    get_filename_component(pf_file_name_we ${pf_file} NAME_WE)

    # Create the test name by combining the path and name. if the path is not empty,
    # it will be included in the test name.
    if (pf_file_path)
        string(REPLACE "/" "-" pf_file_path "${pf_file_path}")
        set(test_name "${pf_file_path}-${pf_file_name_we}")
    else()
        set(test_name "${pf_file_name_we}")
    endif()

    # Add the test to the CMake project.
    add_pfunit_ctest(${test_name}
        TEST_SOURCES ${pf_file}
        LINK_LIBRARIES PkgConfig::neko MPI::MPI_Fortran neko-top-a
    )
    set_target_properties(${test_name} PROPERTIES
        FOLDER "Tests/Unit"
    )
    set_tests_properties(${test_name} PROPERTIES
        LABELS "unit"
    )
    add_dependencies(Unit ${test_name})
endforeach()
