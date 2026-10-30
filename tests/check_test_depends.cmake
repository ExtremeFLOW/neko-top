# ============================================================================ #
# Lint: every entry of a test's DEPENDS property must name an existing test.
# ============================================================================ #
# The CTest DEPENDS test property takes test names. CTest silently ignores any
# entry that is not a test, such as a build-target name, so that entry orders
# nothing. This script reads the registered tests with
# `ctest --show-only=json-v1` and fails, listing every such entry, if there is
# one.
#
# Usage:
#   cmake -DCTEST_COMMAND=<ctest> -DTEST_DIR=<build dir> \
#         -P check_test_depends.cmake

cmake_minimum_required(VERSION 3.21)

if (NOT DEFINED CTEST_COMMAND OR NOT DEFINED TEST_DIR)
    message(FATAL_ERROR "Usage: cmake -DCTEST_COMMAND=<ctest> "
        "-DTEST_DIR=<build dir> -P ${CMAKE_SCRIPT_MODE_FILE}")
endif()

execute_process(
    COMMAND ${CTEST_COMMAND} --test-dir ${TEST_DIR} --show-only=json-v1
    OUTPUT_VARIABLE ctest_json
    ERROR_VARIABLE ctest_error
    RESULT_VARIABLE ctest_result
)
if (NOT ctest_result EQUAL 0)
    message(FATAL_ERROR
        "ctest --show-only=json-v1 failed (${ctest_result}):\n${ctest_error}")
endif()

string(JSON tests GET "${ctest_json}" tests)
string(JSON n_tests LENGTH "${tests}")
if (n_tests EQUAL 0)
    message(FATAL_ERROR "No tests found in ${TEST_DIR}.")
endif()
math(EXPR last_test "${n_tests} - 1")

# Collect the names of all registered tests.
set(test_names "")
foreach(i RANGE ${last_test})
    string(JSON name GET "${tests}" ${i} name)
    list(APPEND test_names "${name}")
endforeach()

# Check every DEPENDS entry against them.
set(n_checked 0)
set(invalid "")
foreach(i RANGE ${last_test})
    string(JSON name GET "${tests}" ${i} name)
    string(JSON properties ERROR_VARIABLE no_properties
        GET "${tests}" ${i} properties)
    if (no_properties)
        continue()
    endif()

    string(JSON n_properties LENGTH "${properties}")
    if (n_properties EQUAL 0)
        continue()
    endif()
    math(EXPR last_property "${n_properties} - 1")

    foreach(j RANGE ${last_property})
        string(JSON property GET "${properties}" ${j} name)
        if (NOT property STREQUAL "DEPENDS")
            continue()
        endif()

        string(JSON depends GET "${properties}" ${j} value)
        string(JSON n_depends LENGTH "${depends}")
        if (n_depends EQUAL 0)
            continue()
        endif()
        math(EXPR last_depend "${n_depends} - 1")

        foreach(k RANGE ${last_depend})
            string(JSON depend GET "${depends}" ${k})
            math(EXPR n_checked "${n_checked} + 1")
            if (NOT depend IN_LIST test_names)
                string(APPEND invalid "\n  ${name}: DEPENDS \"${depend}\"")
            endif()
        endforeach()
    endforeach()
endforeach()

if (NOT invalid STREQUAL "")
    message(FATAL_ERROR
        "DEPENDS entries that are not test names (CTest ignores them):"
        "${invalid}\n"
        "Order build steps with add_dependencies() and test setup steps with "
        "FIXTURES_SETUP/FIXTURES_REQUIRED.")
endif()

message(STATUS "Checked ${n_checked} DEPENDS entries across ${n_tests} tests: "
    "all name existing tests.")
