# ============================================================================ #
# Lint: the Neko-TOP library must call MPI only through the mpi_f08 interfaces.
# ============================================================================ #
# `implicit none` does not forbid calling a procedure without an explicit
# interface. An MPI routine missing from the `use mpi_f08, only:` list of the
# file calling it is therefore still accepted, but the call links to the
# legacy F77 binding, e.g. `mpi_barrier_` instead of `mpi_barrier_f08_`. That
# binding expects integer handles and a mandatory `ierror` argument; a call
# without `ierror` makes it store through whatever address is left in the
# argument register, corrupting memory. A submodule does not take `implicit
# none` from its parent module, so there MPI_Wtime is even typed integer.
#
# This script lists the undefined symbols of the library with `nm` and fails,
# naming each object file, if any of them is an F77 MPI binding: a lower-case
# `mpi_..._` symbol other than the mpi_f08 bindings (`..._f08_`, `..._f08ts_`),
# the `mpi_fortran_..._` constants such as MPI_IN_PLACE, and the
# `mpi_sizeof_..._` procedures of MPI_Sizeof, which has no F77 binding.
#
# Usage:
#   cmake -DNM=<nm> -DLIBRARY=<library> -P check_mpi_bindings.cmake

cmake_minimum_required(VERSION 3.21)

if (NOT DEFINED NM OR NOT DEFINED LIBRARY)
    message(FATAL_ERROR "Usage: cmake -DNM=<nm> -DLIBRARY=<library> "
        "-P ${CMAKE_SCRIPT_MODE_FILE}")
endif()
if (NOT NM)
    message(FATAL_ERROR "No nm tool was found (CMAKE_NM is empty or NOTFOUND).")
endif()
if (NOT EXISTS "${LIBRARY}")
    message(FATAL_ERROR "Library ${LIBRARY} does not exist. Build it first.")
endif()

execute_process(
    COMMAND ${NM} -A -u ${LIBRARY}
    OUTPUT_VARIABLE nm_output
    ERROR_VARIABLE nm_error
    RESULT_VARIABLE nm_result
)
if (NOT nm_result EQUAL 0)
    message(FATAL_ERROR
        "${NM} -A -u ${LIBRARY} failed (${nm_result}):\n${nm_error}")
endif()

# Each line reads "<library>:<object>: U <symbol>". Mach-O prefixes the symbol
# with an underscore.
string(REGEX MATCHALL "[^\n]+" nm_lines "${nm_output}")
set(n_checked 0)
set(invalid "")
foreach(line IN LISTS nm_lines)
    if (NOT line MATCHES "^(.*):[ \t]+U[ \t]+([^ \t]+)$")
        continue()
    endif()
    string(REPLACE "${LIBRARY}:" "" object "${CMAKE_MATCH_1}")
    set(symbol "${CMAKE_MATCH_2}")
    math(EXPR n_checked "${n_checked} + 1")

    if (symbol MATCHES "^_?mpi_[a-z0-9_]*_$"
            AND NOT symbol MATCHES "_f08(ts)?_$"
            AND NOT symbol MATCHES "^_?mpi_fortran_"
            AND NOT symbol MATCHES "^_?mpi_sizeof_")
        string(APPEND invalid "\n  ${object}: ${symbol}")
    endif()
endforeach()

if (n_checked EQUAL 0)
    message(FATAL_ERROR "${NM} reported no undefined symbols in ${LIBRARY}.")
endif()

if (NOT invalid STREQUAL "")
    message(FATAL_ERROR
        "MPI calls without an explicit interface, linking to the F77 MPI "
        "bindings, in ${LIBRARY}:${invalid}\n"
        "Add each routine to the `use mpi_f08, only:` list of the file that "
        "calls it.")
endif()

message(STATUS "Checked ${n_checked} undefined symbols in ${LIBRARY}: "
    "no F77 MPI bindings.")
