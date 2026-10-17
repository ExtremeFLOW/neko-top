# Check what the `steady` simulation component reported about the scalar.
#
# The component freezes only the fluid, so it reports on a scalar that has not
# converged. This counts its three diagnostics in the log a test run wrote,
# and fails unless each count is exactly the expected one. An exact count also
# catches a diagnostic repeated every step instead of issued once per run.
#
#   cmake -DLOG_FILE=<log> -DFROZEN=<n> -DAT_END=<n> -DCONVERGED=<n> \
#         [-DEND_TIME=<t>] -P check_steady_log.cmake
#
# FROZEN     Warnings that the fluid was frozen with the scalar not converged.
# AT_END     Warnings that a run is ending with the scalar still not converged.
# CONVERGED  Messages that the scalar converged after the fluid was frozen.
# END_TIME   The time the runs end at, as the log prints it (e.g.
#            1.50000E+00). If given, the last two only count when issued at
#            that time, so that one issued before the end is caught.

foreach(var LOG_FILE FROZEN AT_END CONVERGED)
    if (NOT DEFINED ${var})
        message(FATAL_ERROR "${var} is not set")
    endif()
endforeach()

if (NOT EXISTS "${LOG_FILE}")
    message(FATAL_ERROR "${LOG_FILE} does not exist")
endif()
file(READ "${LOG_FILE}" log)

# The times are printed as `ES12.5`, so `.` and `+` are the only characters to
# escape.
if (DEFINED END_TIME)
    string(REGEX REPLACE "([.+])" "\\\\\\1" t "${END_TIME}")
else()
    set(t "[^ )]+")
endif()

set(FROZEN_PATTERN "WARNING: steady: fluid frozen at t = ")
set(AT_END_PATTERN "WARNING: steady: run ending at t = ${t} ")
set(CONVERGED_PATTERN
    "steady: the scalar has converged since the fluid was frozen")
string(APPEND CONVERGED_PATTERN " [^\n]* at t = ${t}\\)")

set(failed FALSE)
foreach(kind FROZEN AT_END CONVERGED)
    string(REGEX MATCHALL "${${kind}_PATTERN}" matches "${log}")
    list(LENGTH matches count)
    message(STATUS "${kind}: found ${count}, expected ${${kind}}")
    if (NOT count EQUAL "${${kind}}")
        set(failed TRUE)
    endif()
endforeach()

if (failed)
    message(FATAL_ERROR "Unexpected steady diagnostics in ${LOG_FILE}")
endif()
