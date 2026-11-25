# ============================================================================ #
# Lint: every add_pfunit_ctest() TEST_SOURCES .pf path must be relative to its
# own CMakeLists.txt, not an absolute path.
# ============================================================================ #
# pFUnit's add_pfunit_sources() (see PFUNIT_DIR/include/add_pfunit_sources.cmake)
# treats an absolute TEST_SOURCES path as already being inside the build tree,
# and writes the generated .F90 driver as its sibling. An absolute path that
# actually points into the SOURCE tree -- which is exactly what
# ${CMAKE_CURRENT_SOURCE_DIR}/foo.pf is -- therefore writes the generated
# driver into the source tree instead of the build tree. Two build
# directories on one source tree then share and race to rewrite it.
#
# There is nothing to introspect after configure for this: add_pfunit_ctest()
# consumes its TEST_SOURCES argument immediately, and the literal string
# isn't preserved anywhere a built tree can be queried for. So this scans the
# CMakeLists.txt SOURCE TEXT directly, rather than the build graph the way
# check_test_depends.cmake does.
#
# Usage:
#   cmake -DTESTS_DIR=<path to tests/> -P check_pfunit_test_sources.cmake

cmake_minimum_required(VERSION 3.21)

if (NOT DEFINED TESTS_DIR)
    message(FATAL_ERROR "Usage: cmake -DTESTS_DIR=<path to tests/> "
        "-P ${CMAKE_SCRIPT_MODE_FILE}")
endif()

file(GLOB_RECURSE cmake_files "${TESTS_DIR}/CMakeLists.txt")
list(LENGTH cmake_files n_files)
if (n_files EQUAL 0)
    message(FATAL_ERROR "No CMakeLists.txt files found under ${TESTS_DIR}.")
endif()

set(n_checked 0)
set(invalid "")
foreach(cmake_file ${cmake_files})
    file(READ "${cmake_file}" contents)

    # Collapse line breaks so a TEST_SOURCES keyword and its value, which are
    # often on separate lines, land in one regex match.
    string(REGEX REPLACE "[\r\n]+" " " flattened "${contents}")

    # Every standalone TEST_SOURCES keyword (not a suffix of a longer name
    # such as SHARED_TEST_SOURCES or PFUNIT_TEST_SOURCES -- CMake's regex
    # engine has no \b, hence the explicit non-identifier-or-start guard)
    # immediately followed by its first value token, up to the next
    # whitespace or a closing paren. Every real instance of this bug -- past
    # and present -- has the absolute path as that first token.
    string(REGEX MATCHALL "(^|[^A-Za-z_])TEST_SOURCES[ \t]+[^ \t)]+"
        matches "${flattened}")

    foreach(match ${matches})
        string(REGEX REPLACE "^.*TEST_SOURCES[ \t]+" "" token "${match}")
        if (NOT token MATCHES "\\.pf$")
            continue()
        endif()
        math(EXPR n_checked "${n_checked} + 1")
        if (token MATCHES "^/" OR token MATCHES "CMAKE_CURRENT_SOURCE_DIR")
            string(APPEND invalid "\n  ${cmake_file}: TEST_SOURCES \"${token}\"")
        endif()
    endforeach()
endforeach()

if (n_checked EQUAL 0)
    message(FATAL_ERROR
        "No add_pfunit_ctest TEST_SOURCES .pf entries found under "
        "${TESTS_DIR}.")
endif()

if (NOT invalid STREQUAL "")
    message(FATAL_ERROR
        "add_pfunit_ctest TEST_SOURCES entries that write their generated "
        "driver into the SOURCE tree (absolute, or "
        "\${CMAKE_CURRENT_SOURCE_DIR}-prefixed):${invalid}\n"
        "Pass a path relative to the CMakeLists.txt's own directory "
        "instead, so pFUnit's add_pfunit_sources() generates the driver "
        "into the build tree.")
endif()

message(STATUS "Checked ${n_checked} add_pfunit_ctest TEST_SOURCES .pf "
    "entries across ${n_files} files: none are absolute.")
