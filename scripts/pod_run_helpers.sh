#!/bin/bash

function pod_python_exe() {
    command -v python
}

function pod_prepend_pathlike() {
    local var_name=$1
    local new_entry=$2
    local current_value

    [ -n "${new_entry}" ] || return 0
    [ -d "${new_entry}" ] || return 0

    current_value="${!var_name-}"

    case ":${current_value}:" in
        *":${new_entry}:"*) ;;
        *)
            if [ -n "${current_value}" ]; then
                printf -v "${var_name}" '%s:%s' "${new_entry}" "${current_value}"
            else
                printf -v "${var_name}" '%s' "${new_entry}"
            fi
            export "${var_name}"
            ;;
    esac
}

function pod_append_pathlike() {
    local var_name=$1
    local new_entry=$2
    local current_value

    [ -n "${new_entry}" ] || return 0
    [ -d "${new_entry}" ] || return 0

    current_value="${!var_name-}"

    case ":${current_value}:" in
        *":${new_entry}:"*) ;;
        *)
            if [ -n "${current_value}" ]; then
                printf -v "${var_name}" '%s:%s' "${current_value}" "${new_entry}"
            else
                printf -v "${var_name}" '%s' "${new_entry}"
            fi
            export "${var_name}"
            ;;
    esac
}

function pod_check_python_import() {
    local module_name=$1
    local pyexe

    pyexe=$(pod_python_exe) || {
        echo "Error: could not find python in PATH." >&2
        return 1
    }

    "${pyexe}" - <<PY
import traceback
try:
    import ${module_name}
except Exception:
    traceback.print_exc()
    raise SystemExit(1)
PY
}

function pod_find_repo_root() {
    local start_dir=$1
    local root_dir

    root_dir=$(realpath "${start_dir}")

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
        if [ -n "${candidate}" ] && [ -x "${candidate}/bin/adios2-config" ]; then
            printf '%s\n' "$(realpath "${candidate}")"
            return 0
        fi
    done

    if command -v adios2-config >/dev/null 2>&1; then
        dirname "$(dirname "$(realpath "$(command -v adios2-config)")")"
        return 0
    fi

    return 1
}

function pod_collect_mpi_libdirs() {
    local compiler
    local launcher
    local show_line
    local token
    local dir

    for compiler in "${MPICC:-}" mpicc mpicc.openmpi mpicc.mpich; do
        [ -n "${compiler}" ] || continue
        command -v "${compiler}" >/dev/null 2>&1 || continue

        show_line=$("${compiler}" -show 2>/dev/null || true)
        [ -n "${show_line}" ] || continue

        for token in ${show_line}; do
            case "${token}" in
                -L*)
                    dir="${token#-L}"
                    [ -d "${dir}" ] && printf '%s\n' "${dir}"
                    ;;
            esac
        done
    done

    for launcher in "${MPI_HOME:-}" "${I_MPI_ROOT:-}" "${OPENMPI_HOME:-}" "${MPICH_DIR:-}"; do
        [ -n "${launcher}" ] || continue
        [ -d "${launcher}/lib" ] && printf '%s\n' "${launcher}/lib"
        [ -d "${launcher}/lib64" ] && printf '%s\n' "${launcher}/lib64"
    done

    for launcher in mpirun mpiexec srun; do
        command -v "${launcher}" >/dev/null 2>&1 || continue
        dir="$(dirname "$(dirname "$(realpath "$(command -v "${launcher}")")")")"
        [ -d "${dir}/lib" ] && printf '%s\n' "${dir}/lib"
        [ -d "${dir}/lib64" ] && printf '%s\n' "${dir}/lib64"
    done
}

function pod_ensure_mpi_python() {
    local dir
    local mpi4py_site
    local pyver
    local pyexe

    pyexe=$(pod_python_exe) || {
        echo "Error: could not find python in PATH." >&2
        return 1
    }

    if "${pyexe}" -c 'from mpi4py import MPI' >/dev/null 2>&1; then
        return 0
    fi

    pyver=$("${pyexe}" -c 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}")')

    if [ -n "${VIRTUAL_ENV:-}" ]; then
        pod_prepend_pathlike LD_LIBRARY_PATH "${VIRTUAL_ENV}/lib"
        pod_prepend_pathlike PYTHONPATH "${VIRTUAL_ENV}/lib/python${pyver}/site-packages"
    fi

    if [ -n "${CONDA_PREFIX:-}" ]; then
        pod_prepend_pathlike LD_LIBRARY_PATH "${CONDA_PREFIX}/lib"
        pod_prepend_pathlike PYTHONPATH "${CONDA_PREFIX}/lib/python${pyver}/site-packages"
    fi

    while IFS= read -r dir; do
        [ -n "${dir}" ] || continue
        pod_prepend_pathlike LD_LIBRARY_PATH "${dir}"
    done < <(pod_collect_mpi_libdirs | awk '!seen[$0]++')

    mpi4py_site=$("${pyexe}" -c 'import sysconfig; print(sysconfig.get_paths().get("purelib",""))' 2>/dev/null || true)
    [ -n "${mpi4py_site}" ] && pod_prepend_pathlike PYTHONPATH "${mpi4py_site}"

    "${pyexe}" -c 'from mpi4py import MPI' >/dev/null 2>&1 || {
        echo "Error: mpi4py is installed but cannot load the MPI runtime library." >&2
        echo "python:          $(command -v python 2>/dev/null || echo '<not found>')" >&2
        echo "mpicc:           $(command -v mpicc 2>/dev/null || echo '<not found>')" >&2
        echo "mpirun:          $(command -v mpirun 2>/dev/null || echo '<not found>')" >&2
        echo "LD_LIBRARY_PATH: ${LD_LIBRARY_PATH:-<unset>}" >&2
        pod_check_python_import mpi4py >&2 || true
        return 1
    }
}

function pod_ensure_adios2_python() {
    local root_dir=$1
    local adios2_root
    local pyver
    local dir
    local site_dir
    local pyexe

    pyexe=$(pod_python_exe) || {
        echo "Error: could not find python in PATH." >&2
        return 1
    }

    adios2_root=$(pod_resolve_adios2_root "${root_dir}") || {
        echo "Error: could not locate ADIOS2 for the POD Python driver." >&2
        return 1
    }

    pyver=$("${pyexe}" -c 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}")')

    export ADIOS2_PATH="${adios2_root}"
    export ADIOS2_DIR="${adios2_root}"

    pod_prepend_pathlike LD_LIBRARY_PATH "${ADIOS2_PATH}/lib"
    pod_prepend_pathlike LD_LIBRARY_PATH "${ADIOS2_PATH}/lib64"
    pod_prepend_pathlike PATH "${ADIOS2_PATH}/bin"

    pod_prepend_pathlike PYTHONPATH "${ADIOS2_PATH}/lib/python${pyver}/site-packages"
    pod_prepend_pathlike PYTHONPATH "${ADIOS2_PATH}/lib64/python${pyver}/site-packages"

    while IFS= read -r dir; do
        [ -n "${dir}" ] || continue
        site_dir="$(dirname "$(dirname "${dir}")")"
        pod_prepend_pathlike PYTHONPATH "${site_dir}"
    done < <(find "${ADIOS2_PATH}" -type d -path "*/python${pyver}/site-packages/adios2/bindings" 2>/dev/null | awk '!seen[$0]++')

    "${pyexe}" -c 'import adios2.bindings' >/dev/null 2>&1 && return 0

    echo "Error: adios2.bindings is still not importable." >&2
    echo "ADIOS2_PATH:     ${ADIOS2_PATH}" >&2
    echo "PYTHONPATH:      ${PYTHONPATH:-<unset>}" >&2
    echo "LD_LIBRARY_PATH: ${LD_LIBRARY_PATH:-<unset>}" >&2
    echo "Detailed Python error:" >&2
    pod_check_python_import adios2.bindings >&2 || true
    return 1
}

function pod_prepare_python_runtime() {
    local root_dir=$1
    local pyexe

    pyexe=$(pod_python_exe) || {
        echo "Error: could not find python in PATH." >&2
        return 1
    }

    pod_ensure_mpi_python || return 1
    pod_ensure_adios2_python "${root_dir}" || return 1

    # Optional high-level dependency checks for POD Python side
    "${pyexe}" -c 'import numpy' >/dev/null 2>&1 || {
        echo "Error: numpy is not installed in the active Python environment." >&2
        echo "Active python: ${pyexe}" >&2
        echo "Install it with: python -m pip install numpy" >&2
        return 1
    }
}

function pod_print_runtime_env() {
    echo "Using Python:       $(command -v python 2>/dev/null || echo '<not found>')"
    echo "Using mpicc:        $(command -v mpicc 2>/dev/null || echo '<not found>')"
    echo "Using mpirun:       $(command -v mpirun 2>/dev/null || echo '<not found>')"
    echo "Using srun:         $(command -v srun 2>/dev/null || echo '<not found>')"
    echo "CONDA_PREFIX:       ${CONDA_PREFIX:-<unset>}"
    echo "VIRTUAL_ENV:        ${VIRTUAL_ENV:-<unset>}"
    echo "ADIOS2_PATH:        ${ADIOS2_PATH:-<unset>}"
    echo "PYTHONPATH:         ${PYTHONPATH:-<unset>}"
    echo "LD_LIBRARY_PATH:    ${LD_LIBRARY_PATH:-<unset>}"
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

    local total_ranks=$((py_ranks + neko_ranks))
    local conf_file
    local abs_case_path
    local abs_py_script
    local abs_neko_exe

    abs_case_path=$(realpath "${case_path}")
    abs_py_script=$(realpath "${py_script}")
    abs_neko_exe=$(realpath "${neko_exe}")
    conf_file="$(mktemp "${TMPDIR:-/tmp}/pod_mpmd.XXXXXX.conf")"

    {
        for ((i=0; i<py_ranks; i++)); do
            echo "${i} /usr/bin/env NEKO_COMM_ID=1 NEKO_CTRL_PEER_ROOT=${py_ranks} python ${abs_py_script} ${abs_case_path}"
        done

        for ((i=py_ranks; i<total_ranks; i++)); do
        	   echo "${i} /bin/bash -lc 'sleep 20; exec /usr/bin/env NEKO_COMM_ID=0 NEKO_CTRL_PEER_ROOT=0 \"${abs_neko_exe}\" \"${abs_case_path}\"'"
        done
    } > "${conf_file}"

    echo "Running MPMD job with ${total_ranks} ranks"
    echo "Config:  ${conf_file}"
    echo "Output:  ${log_file}"

    srun --multi-prog "${conf_file}" > "${log_file}" 2>&1
    local rc=$?

    rm -f "${conf_file}"
    return ${rc}
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
