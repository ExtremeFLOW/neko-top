# Check that the MMA and the design HDF5 checkpoints refuse, on every rank, to
# load into a run they do not match, rather than overrunning a buffer, failing
# on some ranks only or hanging.
#
#   cmake -DDRIVER=<exe> -DMPIEXEC=<launcher> -DNUMPROC_FLAG=<flag> \
#         [-DPREFLAGS=<flags>] -DCASES=<rank_count|uneven> \
#         -DTIME_LIMIT=<seconds> -P check_checkpoint_rank_change.cmake
#
# DRIVER      The checkpoint_rank_change executable.
# MPIEXEC     The MPI launcher, with NUMPROC_FLAG and PREFLAGS as for pFUnit.
# CASES       `rank_count`: loads on fewer and on more ranks than the save.
#             `uneven`: a load on as many ranks with another layout, which
#             needs three ranks.
# TIME_LIMIT  The time limit of each run, in seconds. A run that does not stop
#             within it hangs, which fails the test.
#
# Each case saves both checkpoints on some ranks and loads each of them on as
# many, which must complete, so that the checkpoint is known to be sound. It
# then loads each on others, which must stop with the expected error on every
# rank: the driver makes each rank wait for the others before it stops, so
# the output holds one message per rank, and an error raised on some ranks
# only hangs.

foreach(var DRIVER MPIEXEC NUMPROC_FLAG CASES TIME_LIMIT)
    if (NOT DEFINED ${var})
        message(FATAL_ERROR "${var} is not set")
    endif()
endforeach()

set(checkpoint "checkpoint_rank_change.h5")
set(rank_count_error "checkpoint was written with a different number of ranks")
set(failed FALSE)

# Run the driver on `ranks` ranks and return its exit code and output.
function(run_driver ranks mode layout result_var output_var)
    execute_process(
        COMMAND ${MPIEXEC} ${PREFLAGS} ${NUMPROC_FLAG} ${ranks}
            ${DRIVER} ${mode} ${checkpoint} ${layout}
        RESULT_VARIABLE result
        OUTPUT_VARIABLE output
        ERROR_VARIABLE output
        TIMEOUT ${TIME_LIMIT}
    )
    set(${result_var} "${result}" PARENT_SCOPE)
    set(${output_var} "${output}" PARENT_SCOPE)
endfunction()

# Count the occurrences of `pattern` in `text`.
function(count_matches pattern text count_var)
    string(REGEX MATCHALL "${pattern}" matches "${text}")
    list(LENGTH matches count)
    set(${count_var} ${count} PARENT_SCOPE)
endfunction()

# Save on `save_ranks` ranks. Then load each checkpoint on as many ranks,
# which must complete, and on `load_ranks` ranks with the given layout, which
# must stop with `expected` from every loading rank.
function(check_case save_ranks load_ranks layout expected)
    set(case "saved on ${save_ranks}, loaded on ${load_ranks} ${layout}")
    message(STATUS "Case: ${case}")

    file(REMOVE ${checkpoint})
    run_driver(${save_ranks} save even result output)
    count_matches("save completed on rank" "${output}" completed)
    if (NOT result EQUAL 0 OR NOT completed EQUAL save_ranks)
        message(STATUS "${output}")
        message(STATUS "FAIL: save (exit: ${result}, completed on "
            "${completed} of ${save_ranks} ranks)")
        set(failed TRUE PARENT_SCOPE)
        return()
    endif()

    foreach(kind mma design)
        run_driver(${save_ranks} load_${kind} even result output)
        count_matches("load_${kind} completed on rank" "${output}" completed)
        if (NOT result EQUAL 0 OR NOT completed EQUAL save_ranks)
            message(STATUS "${output}")
            message(STATUS "FAIL: ${kind} load as saved (exit: ${result}, "
                "completed on ${completed} of ${save_ranks} ranks)")
            set(failed TRUE PARENT_SCOPE)
            continue()
        endif()
        message(STATUS "${kind}: loads as saved, on ${save_ranks} ranks")

        run_driver(${load_ranks} load_${kind} ${layout} result output)
        count_matches("${kind}: ${expected}" "${output}" errors)
        count_matches("load_${kind} completed" "${output}" completed)
        if (NOT result MATCHES "^[0-9]+$")
            # A hang past the time limit, or a crash
            message(STATUS "${output}")
            message(STATUS "FAIL: ${kind} load: ${result}")
            set(failed TRUE PARENT_SCOPE)
        elseif (result EQUAL 0 OR NOT errors EQUAL load_ranks
                OR NOT completed EQUAL 0)
            message(STATUS "${output}")
            message(STATUS "FAIL: ${kind} load (exit: ${result}, expected "
                "error on ${errors} of ${load_ranks} ranks, completed on "
                "${completed})")
            set(failed TRUE PARENT_SCOPE)
        else()
            message(STATUS "${kind}: error on ${errors} of ${load_ranks} "
                "ranks, as expected")
        endif()
    endforeach()
endfunction()

if (CASES STREQUAL "rank_count")
    # Fewer ranks than the save: the stored per-rank sizes no longer fit.
    check_case(2 1 even "${rank_count_error}")

    # More ranks than the save.
    check_case(1 2 even "${rank_count_error}")
elseif (CASES STREQUAL "uneven")
    # The same number of ranks, but one rank holds a different share. Only
    # some ranks see their own size differ, so all must agree before reading
    # on.
    check_case(3 3 uneven "mismatch in n during HDF5 read")
else()
    message(FATAL_ERROR "Unknown CASES: '${CASES}'")
endif()

file(REMOVE ${checkpoint})

if (failed)
    message(FATAL_ERROR "A checkpoint was loaded into a run it does not match")
endif()
