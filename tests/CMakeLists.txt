# Setup of testing framework using pFUnit as the library
# and CMake as the build system.
#
include(CTest)

if (NOT ${BUILD_TESTING})
  message(STATUS "Tests are disabled")
  return()
endif()

find_package(MPI REQUIRED)
find_package(PkgConfig REQUIRED)

pkg_check_modules(neko REQUIRED IMPORTED_TARGET neko)

# ============================================================================ #
# Setup the Neko tests
# ============================================================================ #

# Find all .pf files in the test directory
set(NEKO_DIR ${CMAKE_SOURCE_DIR}/external/neko CACHE PATH
"Path to the Neko installation")

# Recursively find all .pf files in the test directory of Neko
file(GLOB_RECURSE neko_tests ${NEKO_DIR}/tests/*.pf)

foreach(test ${neko_tests})
  get_filename_component(test_name ${test} NAME_WE)

  add_pfunit_ctest(neko_${test_name}
    TEST_SOURCES ${test}
    LINK_LIBRARIES  PkgConfig::neko MPI::MPI_Fortran
    MAX_PES 2
  )
endforeach()

# ============================================================================ #
# Setup the Local Neko-top tests
# ============================================================================ #

# # Recursively find all .pf files in the test directory of Neko
# file(GLOB_RECURSE local_neko_tests ${CMAKE_CURRENT_SOURCE_DIR}/*.pf)

# foreach(test ${local_neko_tests})
#   get_filename_component(test_name ${test} NAME_WE)

#   add_pfunit_ctest(neko-top_${test_name}
#     TEST_SOURCES ${test}
#     LINK_LIBRARIES  PkgConfig::neko MPI::MPI_Fortran
#     MAX_PES 2
#   )
# endforeach()
