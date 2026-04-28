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

function pod_parse_first_int() {
    local value=${1:-}

    if [[ "${value}" =~ ^([0-9]+) ]]; then
        printf '%s\n' "${BASH_REMATCH[1]}"
        return 0
    fi

    return 1
}

function pod_csv_count() {
    local csv=${1:-}
    local -a values=()

    if [ -z "${csv}" ]; then
        printf '0\n'
        return 0
    fi

    IFS=',' read -r -a values <<< "${csv}"
    printf '%s\n' "${#values[@]}"
}

function pod_current_allowed_cpu_list() {
    local cpu_list

    cpu_list=$(awk -F: '/^Cpus_allowed_list:/ {gsub(/[[:space:]]/, "", $2); print $2; exit}' \
        /proc/self/status 2>/dev/null)

    if [ -n "${cpu_list}" ]; then
        printf '%s\n' "${cpu_list}"
        return 0
    fi

    return 1
}

function pod_cpu_list_contains() {
    local cpu_list=$1
    local cpu=$2
    local -a spans=()
    local span
    local first
    local last

    if [ -z "${cpu_list}" ] || [ -z "${cpu}" ]; then
        return 1
    fi

    IFS=',' read -r -a spans <<< "${cpu_list}"
    for span in "${spans[@]}"; do
        span=${span//[[:space:]]/}

        if [[ "${span}" =~ ^([0-9]+)-([0-9]+)$ ]]; then
            first=${BASH_REMATCH[1]}
            last=${BASH_REMATCH[2]}

            if [ "${cpu}" -ge "${first}" ] && [ "${cpu}" -le "${last}" ]; then
                return 0
            fi
        elif [[ "${span}" =~ ^([0-9]+)$ ]]; then
            if [ "${cpu}" -eq "${BASH_REMATCH[1]}" ]; then
                return 0
            fi
        fi
    done

    return 1
}

function pod_requested_cores_available() {
    local allowed_cpu_list=$1
    local requested_csv=$2
    local -a requested=()
    local cpu

    if [ -z "${requested_csv}" ]; then
        return 0
    fi

    IFS=',' read -r -a requested <<< "${requested_csv}"
    for cpu in "${requested[@]}"; do
        cpu=${cpu//[[:space:]]/}
        [ -z "${cpu}" ] && continue

        if ! pod_cpu_list_contains "${allowed_cpu_list}" "${cpu}"; then
            return 1
        fi
    done

    return 0
}

function pod_launch_shared_mpi_bound() {
    local case_path=$1
    local py_script=$2
    local neko_exe=$3
    local py_ranks=$4
    local neko_ranks=$5
    local log_file=$6

    local node_count
    local pyexe
    local startup_delay
    local neko_ranks_per_node
    local py_ranks_per_node
    local gpu_core_count
    local cpu_core_count
    local total_ranks
    local expected_total_ranks
    local ranks_per_node
    local abs_case_path
    local abs_py_script
    local abs_neko_exe
    local launcher_file
    local conf_file
    local allowed_cpu_list
    local taskset_mode
    local use_taskset
    local srun_cpu_bind
    local errexit_was_set
    local rc
    local node_idx
    local node_rank_base
    local world_rank
    local i

    node_count=$(pod_parse_first_int "${SLURM_JOB_NUM_NODES:-${SLURM_NNODES:-}}") || {
        echo "Error: could not determine the number of Slurm nodes." >&2
        return 1
    }

    pyexe=$(pod_python_exe) || {
        echo "Error: could not find python in PATH." >&2
        return 1
    }

    startup_delay=${POD_NEKO_STARTUP_DELAY:-20}
    total_ranks=$((py_ranks + neko_ranks))

    neko_ranks_per_node=${NEKO_RANKS_PER_NODE:-}
    if [ -z "${neko_ranks_per_node}" ]; then
        if [ $((neko_ranks % node_count)) -ne 0 ]; then
            echo "Error: NEKO_RANKS=${neko_ranks} is not divisible by ${node_count} nodes." >&2
            return 1
        fi
        neko_ranks_per_node=$((neko_ranks / node_count))
    fi

    py_ranks_per_node=${PY_RANKS_PER_NODE:-}
    if [ -z "${py_ranks_per_node}" ]; then
        if [ $((py_ranks % node_count)) -ne 0 ]; then
            echo "Error: PY_RANKS=${py_ranks} is not divisible by ${node_count} nodes." >&2
            return 1
        fi
        py_ranks_per_node=$((py_ranks / node_count))
    fi

    if [ $((neko_ranks_per_node * node_count)) -ne "${neko_ranks}" ]; then
        echo "Error: NEKO_RANKS_PER_NODE=${neko_ranks_per_node} does not match NEKO_RANKS=${neko_ranks} across ${node_count} nodes." >&2
        return 1
    fi

    if [ $((py_ranks_per_node * node_count)) -ne "${py_ranks}" ]; then
        echo "Error: PY_RANKS_PER_NODE=${py_ranks_per_node} does not match PY_RANKS=${py_ranks} across ${node_count} nodes." >&2
        return 1
    fi

    gpu_core_count=$(pod_csv_count "${GPU_CORES:-}")
    if [ "${gpu_core_count}" -lt "${neko_ranks_per_node}" ]; then
        echo "Error: GPU_CORES defines ${gpu_core_count} cores but ${neko_ranks_per_node} GPU ranks per node were requested." >&2
        return 1
    fi

    cpu_core_count=$(pod_csv_count "${CPU_CORES:-}")
    if [ "${cpu_core_count}" -lt "${py_ranks_per_node}" ]; then
        echo "Error: CPU_CORES defines ${cpu_core_count} cores but ${py_ranks_per_node} Python ranks per node were requested." >&2
        return 1
    fi

    ranks_per_node=$((neko_ranks_per_node + py_ranks_per_node))
    expected_total_ranks=$((ranks_per_node * node_count))
    if [ "${expected_total_ranks}" -ne "${total_ranks}" ]; then
        echo "Error: expected ${expected_total_ranks} total ranks from the per-node layout, but got ${total_ranks}." >&2
        return 1
    fi

    taskset_mode=${POD_TASKSET_MODE:-auto}
    use_taskset=0
    srun_cpu_bind=cores

    case "${taskset_mode}" in
        auto)
            if allowed_cpu_list=$(pod_current_allowed_cpu_list); then
                if pod_requested_cores_available "${allowed_cpu_list}" "${GPU_CORES:-}" &&
                   pod_requested_cores_available "${allowed_cpu_list}" "${CPU_CORES:-}"; then
                    use_taskset=1
                else
                    echo "Warning: requested taskset CPU bindings are outside the current CPU allocation (${allowed_cpu_list})."
                    echo "         Falling back to Slurm CPU binding for this run."
                fi
            else
                echo "Warning: could not determine the current CPU allocation."
                echo "         Falling back to Slurm CPU binding for this run."
            fi
            ;;
        always)
            allowed_cpu_list=$(pod_current_allowed_cpu_list || true)
            if [ -n "${allowed_cpu_list}" ] && \
               { ! pod_requested_cores_available "${allowed_cpu_list}" "${GPU_CORES:-}" || \
                 ! pod_requested_cores_available "${allowed_cpu_list}" "${CPU_CORES:-}"; }; then
                echo "Error: requested taskset CPU bindings are outside the current CPU allocation (${allowed_cpu_list})." >&2
                echo "       Use POD_TASKSET_MODE=auto to allow Slurm binding, or request a larger CPU allocation." >&2
                return 1
            fi
            use_taskset=1
            ;;
        never)
            use_taskset=0
            ;;
        *)
            echo "Error: POD_TASKSET_MODE must be one of: auto, always, never." >&2
            return 1
            ;;
    esac

    if [ "${use_taskset}" -eq 1 ]; then
        srun_cpu_bind=none
    fi

    abs_case_path=$(realpath "${case_path}")
    abs_py_script=$(realpath "${py_script}")
    abs_neko_exe=$(realpath "${neko_exe}")
    launcher_file="$(mktemp "${TMPDIR:-/tmp}/pod_mpmd_launcher.XXXXXX.sh")"
    conf_file="$(mktemp "${TMPDIR:-/tmp}/pod_mpmd.XXXXXX.conf")"

    export POD_CASE_FILE="${abs_case_path}"
    export POD_PY_SCRIPT="${abs_py_script}"
    export POD_NEKO_EXE="${abs_neko_exe}"
    export POD_PYTHON_EXE="${pyexe}"
    export POD_NEKO_STARTUP_DELAY="${startup_delay}"
    export NEKO_RANKS_PER_NODE="${neko_ranks_per_node}"
    export PY_RANKS_PER_NODE="${py_ranks_per_node}"
    export POD_USE_TASKSET="${use_taskset}"

    cat > "${launcher_file}" <<'EOF'
#!/bin/bash
set -euo pipefail

local_rank=${SLURM_LOCALID:?}
role=${1:?}

IFS=',' read -r -a gpu_cores_arr <<< "${GPU_CORES:-}"
IFS=',' read -r -a cpu_cores_arr <<< "${CPU_CORES:-}"

if [ "${role}" = "neko" ]; then
    if [ "${local_rank}" -ge "${NEKO_RANKS_PER_NODE}" ]; then
        echo "Error: local Neko rank ${local_rank} exceeds NEKO_RANKS_PER_NODE=${NEKO_RANKS_PER_NODE}." >&2
        exit 1
    fi

    export ROCR_VISIBLE_DEVICES="${local_rank}"
    sleep "${POD_NEKO_STARTUP_DELAY}"

    if [ "${POD_USE_TASKSET:-0}" = "1" ]; then
        core="${gpu_cores_arr[local_rank]}"
        exec taskset -c "${core}" "${POD_NEKO_EXE}" "${POD_CASE_FILE}"
    fi

    exec "${POD_NEKO_EXE}" "${POD_CASE_FILE}"
fi

python_rank=$((local_rank - NEKO_RANKS_PER_NODE))
if [ "${python_rank}" -lt 0 ] || [ "${python_rank}" -ge "${PY_RANKS_PER_NODE}" ]; then
    echo "Error: local Python rank ${local_rank} is outside the configured CPU-rank range." >&2
    exit 1
fi

if [ "${POD_USE_TASKSET:-0}" = "1" ]; then
    core="${cpu_cores_arr[python_rank]}"
    exec taskset -c "${core}" "${POD_PYTHON_EXE}" "${POD_PY_SCRIPT}" "${POD_CASE_FILE}"
fi

exec "${POD_PYTHON_EXE}" "${POD_PY_SCRIPT}" "${POD_CASE_FILE}"
EOF

    chmod +x "${launcher_file}"

    {
        for ((node_idx=0; node_idx<node_count; node_idx++)); do
            node_rank_base=$((node_idx * ranks_per_node))

            for ((i=0; i<neko_ranks_per_node; i++)); do
                world_rank=$((node_rank_base + i))
                echo "${world_rank} /usr/bin/env NEKO_COMM_ID=0 NEKO_CTRL_PEER_ROOT=${neko_ranks_per_node} ${launcher_file} neko"
            done

            for ((i=0; i<py_ranks_per_node; i++)); do
                world_rank=$((node_rank_base + neko_ranks_per_node + i))
                echo "${world_rank} /usr/bin/env NEKO_COMM_ID=1 NEKO_CTRL_PEER_ROOT=0 ${launcher_file} python"
            done
        done
    } > "${conf_file}"

    echo "Running bound MPMD job with ${total_ranks} ranks"
    echo "Nodes:   ${node_count}"
    echo "Layout:  ${neko_ranks_per_node} Neko GPU ranks/node + ${py_ranks_per_node} Python ranks/node"
    if [ "${use_taskset}" -eq 1 ]; then
        echo "CPU bind: taskset (${GPU_CORES} | ${CPU_CORES})"
    else
        echo "CPU bind: Slurm (--cpu-bind=${srun_cpu_bind})"
    fi
    echo "Config:  ${conf_file}"
    echo "Output:  ${log_file}"

    errexit_was_set=0
    if [[ $- == *e* ]]; then
        errexit_was_set=1
        set +e
    fi

    srun --cpu-bind="${srun_cpu_bind}" --distribution=block:block --multi-prog "${conf_file}" > "${log_file}" 2>&1
    rc=$?

    if [ "${errexit_was_set}" -eq 1 ]; then
        set -e
    fi

    if [ "${rc}" -ne 0 ]; then
        echo "Error: shared MPMD launch failed. See ${log_file}." >&2
    fi

    rm -f "${launcher_file}" "${conf_file}"
    return ${rc}
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
    local pyexe
    local startup_delay
    local errexit_was_set

    abs_case_path=$(realpath "${case_path}")
    abs_py_script=$(realpath "${py_script}")
    abs_neko_exe=$(realpath "${neko_exe}")

    if [ -n "${SLURM_JOB_ID:-}" ] && [ -n "${GPU_CORES:-}" ] && [ -n "${CPU_CORES:-}" ]; then
        pod_launch_shared_mpi_bound "${case_path}" "${py_script}" "${neko_exe}" \
            "${py_ranks}" "${neko_ranks}" "${log_file}"
        return $?
    fi

    pyexe=$(pod_python_exe) || {
        echo "Error: could not find python in PATH." >&2
        return 1
    }
    startup_delay=${POD_NEKO_STARTUP_DELAY:-20}
    conf_file="$(mktemp "${TMPDIR:-/tmp}/pod_mpmd.XXXXXX.conf")"

    {
        for ((i=0; i<py_ranks; i++)); do
            echo "${i} /usr/bin/env NEKO_COMM_ID=1 NEKO_CTRL_PEER_ROOT=${py_ranks} ${pyexe} ${abs_py_script} ${abs_case_path}"
        done

        for ((i=py_ranks; i<total_ranks; i++)); do
            echo "${i} /bin/bash -lc 'sleep ${startup_delay}; exec /usr/bin/env NEKO_COMM_ID=0 NEKO_CTRL_PEER_ROOT=0 \"${abs_neko_exe}\" \"${abs_case_path}\"'"
        done
    } > "${conf_file}"

    echo "Running MPMD job with ${total_ranks} ranks"
    echo "Config:  ${conf_file}"
    echo "Output:  ${log_file}"

    errexit_was_set=0
    if [[ $- == *e* ]]; then
        errexit_was_set=1
        set +e
    fi

    srun --multi-prog "${conf_file}" > "${log_file}" 2>&1
    local rc=$?

    if [ "${errexit_was_set}" -eq 1 ]; then
        set -e
    fi

    if [ "${rc}" -ne 0 ]; then
        echo "Error: shared MPMD launch failed. See ${log_file}." >&2
    fi

    rm -f "${conf_file}"
    return ${rc}
}

function pod_launch_neko_only() {
    local case_path=$1
    local neko_exe=$2
    local neko_ranks=$3
    local log_file=$4
    local mpirun_cmd
    local errexit_was_set

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

    errexit_was_set=0
    if [[ $- == *e* ]]; then
        errexit_was_set=1
        set +e
    fi

    "${mpirun_cmd[@]}" > "${log_file}" 2>&1
    local rc=$?

    if [ "${errexit_was_set}" -eq 1 ]; then
        set -e
    fi

    if [ "${rc}" -ne 0 ]; then
        echo "Error: MPI launch failed. See ${log_file}." >&2
    fi

    return ${rc}
}
