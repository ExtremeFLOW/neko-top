#!/bin/bash

# Generic helpers for coupled Neko/Python MPMD launches. This includes
# additional flags to be passed;
# PYTHON_BIN     The directory where python is installed
# NEKO_LAUNCHER  The mpmd launcher, options are mpirun or srun

# Find which python should be used
function mpmd_python_exe() {
    if [ -n "${PYTHON_BIN:-}" ]; then
        if [ -x "${PYTHON_BIN}" ]; then
            printf '%s\n' "${PYTHON_BIN}"
            return 0
        fi

        if command -v "${PYTHON_BIN}" >/dev/null 2>&1; then
            command -v "${PYTHON_BIN}"
            return 0
        fi

        echo "Error: PYTHON_BIN is set but not executable: ${PYTHON_BIN}" >&2
        return 1
    fi

    command -v python3 2>/dev/null || command -v python 2>/dev/null
}

# Make sure the python being used contains mpi4py
function mpmd_validate_mpi4py() {
    local pyexe

    pyexe=$(mpmd_python_exe) || {
        echo "Error: could not find python3 or python in PATH." >&2
        return 1
    }

    if ! "${pyexe}" -c 'from mpi4py import MPI; print(MPI.Get_version())'
    then
        echo "Error: mpi4py is not usable with ${pyexe}." >&2
        return 1
    fi

    export PYTHON_BIN="${pyexe}"
}

# Select between mpirun or srun
function mpmd_selected_launcher() {
    local requested=${NEKO_LAUNCHER:-auto}

    case "${requested}" in
        auto)
            if command -v mpirun >/dev/null 2>&1; then
                printf 'mpirun\n'
                return 0
            fi

            if command -v srun >/dev/null 2>&1; then
                printf 'srun\n'
                return 0
            fi

            echo "Error: neither mpirun nor srun was found in PATH." >&2
            return 1
            ;;
        mpirun|srun)
            if ! command -v "${requested}" >/dev/null 2>&1; then
                echo "Error: ${requested} was requested via NEKO_LAUNCHER" >&2
                echo "but is not available in PATH." >&2
                return 1
            fi

            printf '%s\n' "${requested}"
            ;;
        *)
            echo "Error: NEKO_LAUNCHER must be one of: auto, mpirun, srun." >&2
            return 1
            ;;
    esac
}

function mpmd_launch_shared() {
    local case_path=$1
    local py_script=$2
    local neko_exe=$3
    local py_ranks=$4
    local neko_ranks=$5
    local log_file=$6
    local total_ranks=$((py_ranks + neko_ranks))
    local abs_case_path
    local abs_case_dir
    local abs_py_script
    local abs_neko_exe
    local pyexe
    local launcher
    local startup_delay
    local py_cmd
    local neko_cmd
    local conf_file
    local rc
    local i
    local -a mpirun_cmd

    if [ "${py_ranks}" -lt 1 ] || [ "${neko_ranks}" -lt 1 ]; then
        echo "Error: both Python and Neko need at least one MPI rank." >&2
        return 1
    fi

    abs_case_path=$(realpath "${case_path}") || return 1
    abs_case_dir=$(dirname "${abs_case_path}")
    abs_py_script=$(realpath "${py_script}") || return 1
    abs_neko_exe=$(realpath "${neko_exe}") || return 1
    pyexe=$(mpmd_python_exe) || return 1
    launcher=$(mpmd_selected_launcher) || return 1
    startup_delay=${NEKO_STARTUP_DELAY:-0}

    # The launcher assigns MPMD roles; it does not split communicators. Neko
    # and the Python peer must perform matching MPI collectives after startup.
    # Python ranks are placed first, so the first Neko world rank is py_ranks.
    printf -v py_cmd \
        'cd %q && exec /usr/bin/env NEKO_COMM_ID=1 NEKO_CTRL_PEER_ROOT=%q ' \
        "${abs_case_dir}" "${py_ranks}"
    printf -v py_cmd '%s%q %q %q' \
        "${py_cmd}" "${pyexe}" "${abs_py_script}" "${abs_case_path}"

    printf -v neko_cmd \
        'cd %q && sleep %q && exec /usr/bin/env NEKO_COMM_ID=0 ' \
        "${abs_case_dir}" "${startup_delay}"
    printf -v neko_cmd '%sNEKO_CTRL_PEER_ROOT=0 %q %q' \
        "${neko_cmd}" "${abs_neko_exe}" "${abs_case_path}"

    if [ "${launcher}" = "mpirun" ]; then
        mpirun_cmd=(
            mpirun
            --tag-output
            -n "${py_ranks}"
            /bin/bash
            -c
            "${py_cmd}"
            :
            -n "${neko_ranks}"
            /bin/bash
            -c
            "${neko_cmd}"
        )

        printf 'Launching shared MPI job:'
        printf ' %q' "${mpirun_cmd[@]}"
        printf '\n'
        "${mpirun_cmd[@]}" > "${log_file}" 2>&1
        return $?
    fi

    conf_file=$(mktemp "${TMPDIR:-/tmp}/mpmd.XXXXXX.conf") || return 1
    {
        for ((i = 0; i < py_ranks; i++)); do
            printf '%s /bin/bash -c %q\n' "${i}" "${py_cmd}"
        done
        for ((i = py_ranks; i < total_ranks; i++)); do
            printf '%s /bin/bash -c %q\n' "${i}" "${neko_cmd}"
        done
    } > "${conf_file}"

    srun --unbuffered --multi-prog "${conf_file}" > "${log_file}" 2>&1
    rc=$?
    rm -f "${conf_file}"
    return ${rc}
}
