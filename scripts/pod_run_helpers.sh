#!/bin/bash

function pod_find_repo_root() {
    local start_dir=$1
    local root_dir

    root_dir="${start_dir}"
    while [ ! -d "${root_dir}/sources" ] && [ "${root_dir}" != "/" ]; do
        root_dir="$(dirname "${root_dir}")"
    done

    if [ ! -d "${root_dir}/sources" ]; then
        echo "Error: could not locate repo root from ${start_dir}" >&2
        return 1
    fi

    printf '%s\n' "${root_dir}"
}

function pod_resolve_case_file() {
    local script_dir=$1
    local requested=${2:-}
    local case_files

    if [ -n "${requested}" ]; then
        printf '%s\n' "${requested}"
        return 0
    fi

    shopt -s nullglob
    case_files=( "${script_dir}"/*.case )
    shopt -u nullglob

    if [ ${#case_files[@]} -eq 1 ]; then
        printf '%s\n' "$(basename "${case_files[0]}")"
        return 0
    fi

    if [ ${#case_files[@]} -eq 0 ]; then
        echo "Error: no .case file found in ${script_dir}" >&2
    else
        echo "Error: multiple .case files found in ${script_dir}" >&2
    fi
    return 1
}

function pod_resolve_adios2_root() {
    local root_dir=$1
    local candidate

    for candidate in \
        "${ADIOS2_PATH:-}" \
        "${ADIOS2_DIR:-}" \
        "${root_dir}/external/adios2" \
        "${MAIN_DIR:-}/external/adios2"; do
        if [ -n "${candidate}" ] && \
           [ -x "${candidate}/bin/adios2-config" ]; then
            printf '%s\n' "${candidate}"
            return 0
        fi
    done

    if command -v adios2-config >/dev/null 2>&1; then
        dirname "$(dirname "$(realpath "$(command -v adios2-config)")")"
        return 0
    fi

    return 1
}

function pod_ensure_adios2_python() {
    local root_dir=$1
    local adios2_root
    local pyver

    if python3 -c 'import adios2.bindings' >/dev/null 2>&1; then
        return 0
    fi

    adios2_root=$(pod_resolve_adios2_root "${root_dir}") || {
        echo "Error: could not locate ADIOS2 for the POD Python driver." >&2
        return 1
    }

    pyver=$(python3 -c \
        'import sys; print(f"{sys.version_info.major}.'\
'{sys.version_info.minor}")')

    export ADIOS2_PATH="${adios2_root}"
    export PYTHONPATH="${ADIOS2_PATH}/lib/python${pyver}/site-packages"\
"${PYTHONPATH:+:${PYTHONPATH}}"
    export LD_LIBRARY_PATH="${ADIOS2_PATH}/lib:${ADIOS2_PATH}/lib64"\
"${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

    python3 -c 'import adios2.bindings' >/dev/null 2>&1 || {
        echo "Error: adios2.bindings is still not importable." >&2
        return 1
    }

    echo "Using ADIOS2_PATH=${ADIOS2_PATH}"
}

function pod_print_runtime_env() {
    echo "Using Python:      $(which python3)"
    echo "Using mpirun:      $(which mpirun)"
    echo "CONDA_PREFIX:      ${CONDA_PREFIX:-<unset>}"
    echo "PYTHONPATH (start):${PYTHONPATH:-<unset>}"
    echo "LD_LIBRARY_PATH:   ${LD_LIBRARY_PATH:-<unset>}"
    echo "---------------------------------------------------------"
}

function pod_case_is_pod() {
    local case_path=$1
    grep -q '"type"[[:space:]]*:[[:space:]]*"pod"' "${case_path}"
}

function pod_launch_shared_mpi() {
    local case_path=$1
    local py_script=$2
    local neko_exe=$3
    local py_ranks=$4
    local neko_ranks=$5
    local log_file=$6
    local mpirun_cmd

    mpirun_cmd=(
        mpirun
        --tag-output
        -n "${py_ranks}"
        /usr/bin/env
        NEKO_COMM_ID=1
        NEKO_CTRL_PEER_ROOT="${py_ranks}"
        python3 "${py_script}" "${case_path}"
        :
        -n "${neko_ranks}"
        /usr/bin/env
        NEKO_COMM_ID=0
        NEKO_CTRL_PEER_ROOT=0
        "${neko_exe}" "${case_path}"
    )

    echo "Launching shared MPI job:"
    printf '  %q' "${mpirun_cmd[@]}"
    printf '\n'
    echo "Output:            ${log_file}"
    "${mpirun_cmd[@]}" > "${log_file}" 2>&1
}

function pod_launch_neko_only() {
    local case_path=$1
    local neko_exe=$2
    local neko_ranks=$3
    local log_file=$4
    local mpirun_cmd

    mpirun_cmd=(
        mpirun
        --tag-output
        -n "${neko_ranks}"
        "${neko_exe}" "${case_path}"
    )

    echo "Launching MPI job:"
    printf '  %q' "${mpirun_cmd[@]}"
    printf '\n'
    echo "Output:            ${log_file}"
    "${mpirun_cmd[@]}" > "${log_file}" 2>&1
}
