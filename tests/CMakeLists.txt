# ============================================================================ #
# Neko-TOP tests
# ============================================================================ #
# This CMakeLists.txt file is used to define the tests for the Neko-TOP
# library. The tests are defined using the pFUnit testing framework, which
# is a Fortran testing framework that is used to test the Neko-TOP library.
#
# ============================================================================ #
# External dependencies
# ============================================================================ #

# Setup of testing framework using pFUnit as the library
# and CMake as the build system.
#
find_package(MPI REQUIRED)
find_package(PkgConfig REQUIRED)

pkg_check_modules(neko REQUIRED IMPORTED_TARGET neko)

# ============================================================================ #
# Setup the Neko-TOP tests
# ============================================================================ #
# A couple of different types of tests can be setup and will be defined in this
# file. The tests are defined in the following way:
#
# 1. Sequential tests: These tests are run on a single processor and are used
#    to test the basic functionality of the Neko-TOP library.
# 2. Parallel tests: These tests are run on multiple processors and are used
#    to test the parallel functionality of the Neko-TOP library.
# 3. Performance tests: These tests are used to test the performance of the
#    Neko-TOP library and are run on multiple processors.

# ---------------------------------------------------------------------------- #
# Sequential tests

# Recursively find all .pf files in the sequential directory.
file(GLOB_RECURSE sequential_tests ${CMAKE_CURRENT_SOURCE_DIR}/sequential/*.pf)

foreach(test ${sequential_tests})
  get_filename_component(test_name ${test} NAME_WE)

  # Add the neko tests
  add_pfunit_ctest(neko-top_seq_${test_name}
    TEST_SOURCES ${test}
    LINK_LIBRARIES  PkgConfig::neko MPI::MPI_Fortran Neko-TOP
  )
endforeach()

# ---------------------------------------------------------------------------- #
# Parallel tests

# Recursively find all .pf files in the parallel directory.
file(GLOB parallel_tests ${CMAKE_CURRENT_SOURCE_DIR}/parallel/*.pf)

foreach(test ${parallel_tests})
  get_filename_component(test_name ${test} NAME_WE)

  # Add the neko tests
  add_pfunit_ctest(neko-top_par_${test_name}
    TEST_SOURCES ${test}
    LINK_LIBRARIES  PkgConfig::neko MPI::MPI_Fortran Neko-TOP
    MAX_PES 4
  )
endforeach()