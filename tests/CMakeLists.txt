# ============================================================================ #
# Neko-TOP tests
# ============================================================================ #
# This CMakeLists.txt file is used to define the tests for the Neko-TOP
# library. The tests are defined using the pFUnit testing framework, which
# is a Fortran testing framework that is used to test the Neko-TOP library.

# ============================================================================ #
# External dependencies
# ============================================================================ #

# Setup of testing framework using pFUnit as the library
# and CMake as the build system.

find_package(MPI REQUIRED)
find_package(PkgConfig REQUIRED)

pkg_check_modules(neko REQUIRED IMPORTED_TARGET neko)

# ============================================================================ #
# Setup the Neko-TOP tests
# ============================================================================ #
# A couple of different types of tests can be setup and will be defined in this
# file. The tests are defined in the following ways:
#
# 1. Simple tests: These test are just a single pfunit test that is completely
#    self-contained and does not require any additional files.
# 2. Complex tests: These tests are more complex and require additional files
#    to be included in the test or require additional setup.
#    Examples of complex tests are MPI tests and tests that require additional
#    files to be included.
#
# Simple tests are located and added automatically, complex tests must have an
# associated CMakeLists.txt file in the test directory. The CMakeLists.txt file
# should define the test and any additional files that are required for the
# test.

# Recursively find all folders in the test directory.
file(GLOB_RECURSE folder_list
    LIST_DIRECTORIES YES
    RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/*
)

foreach(folder ${folder_list})
    if (IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${folder})

        # Check if a CMakeLists.txt file exists in the directory.
        if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${folder}/CMakeLists.txt)
            add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/${folder})

        else()

            # Assign test name as the directory name, with "/" replaced by "-".
            string(REPLACE "/" "-" test_name ${folder})

            # Find all .pf files in the directory and add them as tests.
            file(GLOB pf_files ${CMAKE_CURRENT_SOURCE_DIR}/${folder}/*.pf)
            foreach(pf_file ${pf_files})
                get_filename_component(pf_file_name ${pf_file} NAME_WE)
                file(RELATIVE_PATH pf_file ${CMAKE_CURRENT_SOURCE_DIR} ${pf_file})
                message(STATUS "Adding test: ${pf_file_name}")
                message(STATUS "Adding test: ${pf_file}")
                add_pfunit_ctest(${test_name}-${pf_file_name}
                    TEST_SOURCES ${pf_file}
                    LINK_LIBRARIES PkgConfig::neko MPI::MPI_Fortran Neko-TOP
                )
            endforeach()
        endif()
    endif()
endforeach()
