#!/bin/bash

function mpmd_python_exe() {
    if [ -n "${PYTHON_BIN:-}" ]; then
        if [ -x "${PYTHON_BIN}" ]; then
            printf '%s\n' "${PYTHON_BIN}"
            return 0
        elif command -v "${PYTHON_BIN}" >/dev/null 2>&1; then
            command -v "${PYTHON_BIN}"
            return 0
        fi

        echo "Error: PYTHON_BIN is set but not executable: ${PYTHON_BIN}" >&2
        return 1
    fi

    command -v python3 2>/dev/null || command -v python 2>/dev/null
}

function mpmd_find_repo_root() {
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

function mpmd_runtime_env_file() {
    local root_dir=$1
    local repo_root

    repo_root=$(mpmd_find_repo_root "${root_dir}") || return 1
    printf '%s\n' "${repo_root}/build/pod_runtime.env"
}

function mpmd_source_runtime_env() {
    local root_dir=$1
    local runtime_env

    runtime_env=$(mpmd_runtime_env_file "${root_dir}") || return 1
    if [ ! -f "${runtime_env}" ]; then
        echo "Error: Python runtime env not found: ${runtime_env}" >&2
        echo "Run ./setup.sh -e from the repo root after activating the" >&2
        echo "target Python environment." >&2
        return 1
    fi

    # shellcheck disable=SC1090
    source "${runtime_env}"
    export MPMD_RUNTIME_ENV_FILE="${runtime_env}"
    export POD_RUNTIME_ENV_FILE="${runtime_env}"
}

function mpmd_validate_python_runtime() {
    local root_dir=$1
    local repo_root
    local validator
    local pyexe

    repo_root=$(mpmd_find_repo_root "${root_dir}") || return 1
    validator="${repo_root}/scripts/python/validate_pod_runtime.py"
    if [ ! -f "${validator}" ]; then
        echo "Error: Python runtime validator not found: ${validator}" >&2
        return 1
    fi

    pyexe=$(mpmd_python_exe) || {
        echo "Error: could not find python3 or python in PATH." >&2
        return 1
    }

    if ! "${pyexe}" "${validator}"; then
        echo "Error: the active Python runtime is not valid." >&2
        echo "Python: ${pyexe}" >&2
        echo "Runtime env: ${MPMD_RUNTIME_ENV_FILE:-<unset>}" >&2
        return 1
    fi

    export PYTHON_BIN="${pyexe}"
}

function mpmd_prepare_python_runtime() {
    local root_dir=$1

    mpmd_source_runtime_env "${root_dir}" || return 1
    mpmd_validate_python_runtime "${root_dir}" || return 1
}

function mpmd_ensure_mpi_python() {
    local root_dir=${1:-${MAIN_DIR:-.}}
    mpmd_prepare_python_runtime "${root_dir}"
}

function mpmd_ensure_adios2_python() {
    local root_dir=${1:-${MAIN_DIR:-.}}
    mpmd_prepare_python_runtime "${root_dir}"
}

function mpmd_startup_delay() {
    if [ -n "${NEKO_STARTUP_DELAY:-}" ]; then
        printf '%s\n' "${NEKO_STARTUP_DELAY}"
        return 0
    fi

    printf '%s\n' "${POD_NEKO_STARTUP_DELAY:-20}"
}

function mpmd_taskset_mode() {
    if [ -n "${NEKO_TASKSET_MODE:-}" ]; then
        printf '%s\n' "${NEKO_TASKSET_MODE}"
        return 0
    fi

    printf '%s\n' "${POD_TASKSET_MODE:-auto}"
}

function mpmd_selected_launcher() {
    local requested=${NEKO_LAUNCHER:-auto}

    case "${requested}" in
        auto)
            if [ -n "${SLURM_JOB_ID:-}" ] && command -v srun >/dev/null 2>&1
            then
                printf 'srun\n'
                return 0
            fi

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
        srun|mpirun)
            if ! command -v "${requested}" >/dev/null 2>&1; then
                echo "Error: ${requested} was requested via" >&2
                echo "NEKO_LAUNCHER but is not available in PATH." >&2
                return 1
            fi

            printf '%s\n' "${requested}"
            return 0
            ;;
        *)
            echo "Error: NEKO_LAUNCHER must be one of:" >&2
            echo "auto, srun, mpirun." >&2
            return 1
            ;;
    esac
}

function mpmd_print_runtime_env() {
    echo "Using Python:       $(mpmd_python_exe 2>/dev/null || echo \
'<not found>')"
    echo "Using mpicc:        $(command -v mpicc 2>/dev/null || echo \
'<not found>')"
    echo "Using mpirun:       $(command -v mpirun 2>/dev/null || echo \
'<not found>')"
    echo "Using srun:         $(command -v srun 2>/dev/null || echo \
'<not found>')"
    echo "NEKO_LAUNCHER:      ${NEKO_LAUNCHER:-auto}"
    echo "Selected launcher:  $(mpmd_selected_launcher 2>/dev/null || echo \
'<unavailable>')"
    echo "Python runtime env: ${MPMD_RUNTIME_ENV_FILE:-<unset>}"
    echo "CONDA_PREFIX:       ${CONDA_PREFIX:-<unset>}"
    echo "VIRTUAL_ENV:        ${VIRTUAL_ENV:-<unset>}"
    echo "ADIOS2_PATH:        ${ADIOS2_PATH:-<unset>}"
    echo "PYTHONPATH:         ${PYTHONPATH:-<unset>}"
    echo "LD_LIBRARY_PATH:    ${LD_LIBRARY_PATH:-<unset>}"
    echo "---------------------------------------------------------"
}

function mpmd_check_python_import() {
    local module_name=$1
    local pyexe

    pyexe=$(mpmd_python_exe) || {
        echo "Error: could not find python3 or python in PATH." >&2
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

function mpmd_resolve_case_file() {
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

function mpmd_case_needs_python_driver() {
    local case_path=$1
    grep -q '"type"[[:space:]]*:[[:space:]]*"pod"' "${case_path}"
}

function mpmd_parse_first_int() {
    local value=${1:-}

    if [[ "${value}" =~ ^([0-9]+) ]]; then
        printf '%s\n' "${BASH_REMATCH[1]}"
        return 0
    fi

    return 1
}

function mpmd_csv_count() {
    local csv=${1:-}
    local -a values=()

    if [ -z "${csv}" ]; then
        printf '0\n'
        return 0
    fi

    IFS=',' read -r -a values <<< "${csv}"
    printf '%s\n' "${#values[@]}"
}

function mpmd_current_allowed_cpu_list() {
    local cpu_list

    cpu_list=$(
        awk -F: \
            '/^Cpus_allowed_list:/ {gsub(/[[:space:]]/, "", $2);'\
            ' print $2; exit}' \
            /proc/self/status 2>/dev/null
    )

    if [ -n "${cpu_list}" ]; then
        printf '%s\n' "${cpu_list}"
        return 0
    fi

    return 1
}

function mpmd_cpu_list_contains() {
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

            if [ "${cpu}" -ge "${first}" ] && [ "${cpu}" -le "${last}" ]
            then
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

function mpmd_requested_cores_available() {
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

        if ! mpmd_cpu_list_contains "${allowed_cpu_list}" "${cpu}"; then
            return 1
        fi
    done

    return 0
}

function mpmd_launch_shared_bound() {
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
    local abs_case_dir
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

    node_count=$(
        mpmd_parse_first_int \
            "${SLURM_JOB_NUM_NODES:-${SLURM_NNODES:-}}"
    ) || {
        echo "Error: could not determine the number of Slurm nodes." >&2
        return 1
    }

    pyexe=$(mpmd_python_exe) || {
        echo "Error: could not find python in PATH." >&2
        return 1
    }

    startup_delay=$(mpmd_startup_delay)
    total_ranks=$((py_ranks + neko_ranks))

    neko_ranks_per_node=${NEKO_RANKS_PER_NODE:-}
    if [ -z "${neko_ranks_per_node}" ]; then
        if [ $((neko_ranks % node_count)) -ne 0 ]; then
            echo "Error: NEKO_RANKS=${neko_ranks} is not divisible by" >&2
            echo "${node_count} nodes." >&2
            return 1
        fi
        neko_ranks_per_node=$((neko_ranks / node_count))
    fi

    py_ranks_per_node=${PY_RANKS_PER_NODE:-}
    if [ -z "${py_ranks_per_node}" ]; then
        if [ $((py_ranks % node_count)) -ne 0 ]; then
            echo "Error: PY_RANKS=${py_ranks} is not divisible by" >&2
            echo "${node_count} nodes." >&2
            return 1
        fi
        py_ranks_per_node=$((py_ranks / node_count))
    fi

    if [ $((neko_ranks_per_node * node_count)) -ne "${neko_ranks}" ]; then
        echo "Error: NEKO_RANKS_PER_NODE=${neko_ranks_per_node} does" >&2
        echo "not match NEKO_RANKS=${neko_ranks} across ${node_count}" >&2
        echo "nodes." >&2
        return 1
    fi

    if [ $((py_ranks_per_node * node_count)) -ne "${py_ranks}" ]; then
        echo "Error: PY_RANKS_PER_NODE=${py_ranks_per_node} does not" >&2
        echo "match PY_RANKS=${py_ranks} across ${node_count} nodes." >&2
        return 1
    fi

    gpu_core_count=$(mpmd_csv_count "${GPU_CORES:-}")
    if [ "${gpu_core_count}" -lt "${neko_ranks_per_node}" ]; then
        echo "Error: GPU_CORES defines ${gpu_core_count} cores but" >&2
        echo "${neko_ranks_per_node} GPU ranks per node were" >&2
        echo "requested." >&2
        return 1
    fi

    cpu_core_count=$(mpmd_csv_count "${CPU_CORES:-}")
    if [ "${cpu_core_count}" -lt "${py_ranks_per_node}" ]; then
        echo "Error: CPU_CORES defines ${cpu_core_count} cores but" >&2
        echo "${py_ranks_per_node} Python ranks per node were" >&2
        echo "requested." >&2
        return 1
    fi

    ranks_per_node=$((neko_ranks_per_node + py_ranks_per_node))
    expected_total_ranks=$((ranks_per_node * node_count))
    if [ "${expected_total_ranks}" -ne "${total_ranks}" ]; then
        echo "Error: expected ${expected_total_ranks} total ranks from" >&2
        echo "the per-node layout, but got ${total_ranks}." >&2
        return 1
    fi

    taskset_mode=$(mpmd_taskset_mode)
    use_taskset=0
    srun_cpu_bind=cores

    case "${taskset_mode}" in
        auto)
            if allowed_cpu_list=$(mpmd_current_allowed_cpu_list); then
                if mpmd_requested_cores_available \
                    "${allowed_cpu_list}" "${GPU_CORES:-}" &&
                    mpmd_requested_cores_available \
                    "${allowed_cpu_list}" "${CPU_CORES:-}"
                then
                    use_taskset=1
                else
                    echo "Warning: requested taskset CPU bindings are" >&2
                    echo "outside the current CPU allocation" >&2
                    echo "(${allowed_cpu_list}). Falling back to Slurm" >&2
                    echo "CPU binding for this run." >&2
                fi
            else
                echo "Warning: could not determine the current CPU" >&2
                echo "allocation. Falling back to Slurm CPU binding" >&2
                echo "for this run." >&2
            fi
            ;;
        always)
            allowed_cpu_list=$(mpmd_current_allowed_cpu_list || true)
            if [ -n "${allowed_cpu_list}" ] && {
                ! mpmd_requested_cores_available \
                    "${allowed_cpu_list}" "${GPU_CORES:-}" ||
                ! mpmd_requested_cores_available \
                    "${allowed_cpu_list}" "${CPU_CORES:-}"
            }; then
                echo "Error: requested taskset CPU bindings are" >&2
                echo "outside the current CPU allocation" >&2
                echo "(${allowed_cpu_list})." >&2
                echo "Use NEKO_TASKSET_MODE=auto to allow Slurm" >&2
                echo "binding, or request a larger CPU allocation." >&2
                return 1
            fi
            use_taskset=1
            ;;
        never)
            use_taskset=0
            ;;
        *)
            echo "Error: NEKO_TASKSET_MODE must be one of:" >&2
            echo "auto, always, never." >&2
            return 1
            ;;
    esac

    if [ "${use_taskset}" -eq 1 ]; then
        srun_cpu_bind=none
    fi

    abs_case_path=$(realpath "${case_path}")
    abs_case_dir=$(dirname "${abs_case_path}")
    abs_py_script=$(realpath "${py_script}")
    abs_neko_exe=$(realpath "${neko_exe}")
    launcher_file=$(mktemp "${TMPDIR:-/tmp}/mpmd_launcher.XXXXXX.sh")
    conf_file=$(mktemp "${TMPDIR:-/tmp}/mpmd.XXXXXX.conf")

    export MPMD_CASE_FILE="${abs_case_path}"
    export MPMD_CASE_DIR="${abs_case_dir}"
    export MPMD_PY_SCRIPT="${abs_py_script}"
    export MPMD_NEKO_EXE="${abs_neko_exe}"
    export MPMD_PYTHON_EXE="${pyexe}"
    export NEKO_STARTUP_DELAY="${startup_delay}"
    export NEKO_RANKS_PER_NODE="${neko_ranks_per_node}"
    export PY_RANKS_PER_NODE="${py_ranks_per_node}"
    export MPMD_USE_TASKSET="${use_taskset}"

    cat > "${launcher_file}" <<'EOF'
#!/bin/bash
set -euo pipefail

local_rank=${SLURM_LOCALID:?}
role=${1:?}

IFS=',' read -r -a gpu_cores_arr <<< "${GPU_CORES:-}"
IFS=',' read -r -a cpu_cores_arr <<< "${CPU_CORES:-}"

if [ "${role}" = "neko" ]; then
    if [ "${local_rank}" -ge "${NEKO_RANKS_PER_NODE}" ]; then
        echo "Error: local Neko rank ${local_rank} exceeds" >&2
        echo "NEKO_RANKS_PER_NODE=${NEKO_RANKS_PER_NODE}." >&2
        exit 1
    fi

    export ROCR_VISIBLE_DEVICES="${local_rank}"
    sleep "${NEKO_STARTUP_DELAY}"
    cd "${MPMD_CASE_DIR}"

    if [ "${MPMD_USE_TASKSET:-0}" = "1" ]; then
        core="${gpu_cores_arr[local_rank]}"
        exec taskset -c "${core}" "${MPMD_NEKO_EXE}" \
            "${MPMD_CASE_FILE}"
    fi

    exec "${MPMD_NEKO_EXE}" "${MPMD_CASE_FILE}"
fi

python_rank=$((local_rank - NEKO_RANKS_PER_NODE))
if [ "${python_rank}" -lt 0 ] || [ "${python_rank}" -ge \
"${PY_RANKS_PER_NODE}" ]; then
    echo "Error: local Python rank ${local_rank} is outside the" >&2
    echo "configured CPU-rank range." >&2
    exit 1
fi

cd "${MPMD_CASE_DIR}"

if [ "${MPMD_USE_TASKSET:-0}" = "1" ]; then
    core="${cpu_cores_arr[python_rank]}"
    exec taskset -c "${core}" "${MPMD_PYTHON_EXE}" \
        "${MPMD_PY_SCRIPT}" "${MPMD_CASE_FILE}"
fi

exec "${MPMD_PYTHON_EXE}" "${MPMD_PY_SCRIPT}" "${MPMD_CASE_FILE}"
EOF

    chmod +x "${launcher_file}"

    {
        for ((node_idx=0; node_idx<node_count; node_idx++)); do
            node_rank_base=$((node_idx * ranks_per_node))

            for ((i=0; i<neko_ranks_per_node; i++)); do
                world_rank=$((node_rank_base + i))
                echo "${world_rank} /usr/bin/env NEKO_COMM_ID=0" \
                    "NEKO_CTRL_PEER_ROOT=${neko_ranks_per_node}" \
                    "${launcher_file} neko"
            done

            for ((i=0; i<py_ranks_per_node; i++)); do
                world_rank=$((node_rank_base + neko_ranks_per_node + i))
                echo "${world_rank} /usr/bin/env NEKO_COMM_ID=1" \
                    "NEKO_CTRL_PEER_ROOT=0 ${launcher_file} python"
            done
        done
    } > "${conf_file}"

    echo "Running bound MPMD job with ${total_ranks} ranks"
    echo "Nodes:   ${node_count}"
    echo "Layout:  ${neko_ranks_per_node} Neko GPU ranks/node +"
    echo "         ${py_ranks_per_node} Python ranks/node"
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

    srun --unbuffered --cpu-bind="${srun_cpu_bind}" \
        --distribution=block:block --multi-prog "${conf_file}" \
        > "${log_file}" 2>&1
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

function mpmd_launch_shared() {
    local case_path=$1
    local py_script=$2
    local neko_exe=$3
    local py_ranks=$4
    local neko_ranks=$5
    local log_file=$6

    local total_ranks=$((py_ranks + neko_ranks))
    local conf_file
    local abs_case_path
    local abs_case_dir
    local abs_py_script
    local abs_neko_exe
    local pyexe
    local startup_delay
    local launcher
    local errexit_was_set
    local rc
    local i
    local neko_cmd
    local py_cmd
    local -a mpirun_cmd

    abs_case_path=$(realpath "${case_path}")
    abs_case_dir=$(dirname "${abs_case_path}")
    abs_py_script=$(realpath "${py_script}")
    abs_neko_exe=$(realpath "${neko_exe}")
    launcher=$(mpmd_selected_launcher) || return 1

    if [ "${launcher}" = "srun" ] &&
        [ -n "${SLURM_JOB_ID:-}" ] &&
        [ -n "${GPU_CORES:-}" ] &&
        [ -n "${CPU_CORES:-}" ]
    then
        mpmd_launch_shared_bound \
            "${case_path}" "${py_script}" "${neko_exe}" \
            "${py_ranks}" "${neko_ranks}" "${log_file}"
        return $?
    fi

    pyexe=$(mpmd_python_exe) || {
        echo "Error: could not find python in PATH." >&2
        return 1
    }
    startup_delay=$(mpmd_startup_delay)

    if [ "${launcher}" = "srun" ]; then
        conf_file=$(mktemp "${TMPDIR:-/tmp}/mpmd.XXXXXX.conf")

        {
            for ((i=0; i<py_ranks; i++)); do
                echo "${i} /bin/bash -lc 'cd \"${abs_case_dir}\" &&" \
                    "exec /usr/bin/env NEKO_COMM_ID=1" \
                    "NEKO_CTRL_PEER_ROOT=${py_ranks} \"${pyexe}\"" \
                    "\"${abs_py_script}\" \"${abs_case_path}\"'"
            done

            for ((i=py_ranks; i<total_ranks; i++)); do
                echo "${i} /bin/bash -lc 'cd \"${abs_case_dir}\" &&" \
                    "sleep ${startup_delay};" \
                    "exec /usr/bin/env NEKO_COMM_ID=0" \
                    "NEKO_CTRL_PEER_ROOT=0 \"${abs_neko_exe}\"" \
                    "\"${abs_case_path}\"'"
            done
        } > "${conf_file}"

        echo "Running MPMD job with ${total_ranks} ranks"
        echo "Launcher: srun --multi-prog"
        echo "Config:   ${conf_file}"
        echo "Output:   ${log_file}"

        errexit_was_set=0
        if [[ $- == *e* ]]; then
            errexit_was_set=1
            set +e
        fi

        srun --unbuffered --multi-prog "${conf_file}" \
            > "${log_file}" 2>&1
        rc=$?

        if [ "${errexit_was_set}" -eq 1 ]; then
            set -e
        fi

        if [ "${rc}" -ne 0 ]; then
            echo "Error: shared MPMD launch failed. See ${log_file}." >&2
        fi

        rm -f "${conf_file}"
        return ${rc}
    fi

    printf -v py_cmd \
        'cd %q && exec /usr/bin/env NEKO_COMM_ID=1 ' \
        "${abs_case_dir}"
    printf -v py_cmd \
        '%sNEKO_CTRL_PEER_ROOT=%q %q %q %q' \
        "${py_cmd}" "${py_ranks}" "${pyexe}" "${abs_py_script}" \
        "${abs_case_path}"

    printf -v neko_cmd \
        'cd %q && sleep %q; exec /usr/bin/env NEKO_COMM_ID=0 ' \
        "${abs_case_dir}" "${startup_delay}"
    printf -v neko_cmd \
        '%sNEKO_CTRL_PEER_ROOT=0 %q %q' \
        "${neko_cmd}" "${abs_neko_exe}" "${abs_case_path}"

    mpirun_cmd=(
        mpirun
        --tag-output
        -n "${py_ranks}"
        /bin/bash
        -lc
        "${py_cmd}"
        :
        -n "${neko_ranks}"
        /bin/bash
        -lc
        "${neko_cmd}"
    )

    echo "Launching shared MPI job:"
    printf '  %q' "${mpirun_cmd[@]}"
    printf '\n'
    echo "Output:            ${log_file}"

    errexit_was_set=0
    if [[ $- == *e* ]]; then
        errexit_was_set=1
        set +e
    fi

    "${mpirun_cmd[@]}" > "${log_file}" 2>&1
    rc=$?

    if [ "${errexit_was_set}" -eq 1 ]; then
        set -e
    fi

    if [ "${rc}" -ne 0 ]; then
        echo "Error: shared MPMD launch failed. See ${log_file}." >&2
    fi

    return ${rc}
}

function mpmd_launch_neko_only() {
    local case_path=$1
    local neko_exe=$2
    local neko_ranks=$3
    local log_file=$4
    local launcher
    local errexit_was_set
    local rc
    local -a run_cmd

    launcher=$(mpmd_selected_launcher) || return 1

    if [ "${launcher}" = "srun" ]; then
        run_cmd=(
            srun
            --unbuffered
            -n "${neko_ranks}"
            "${neko_exe}" "${case_path}"
        )
    else
        run_cmd=(
            mpirun
            --tag-output
            -n "${neko_ranks}"
            "${neko_exe}" "${case_path}"
        )
    fi

    echo "Launching ${launcher} job:"
    printf '  %q' "${run_cmd[@]}"
    printf '\n'
    echo "Output:            ${log_file}"

    errexit_was_set=0
    if [[ $- == *e* ]]; then
        errexit_was_set=1
        set +e
    fi

    "${run_cmd[@]}" > "${log_file}" 2>&1
    rc=$?

    if [ "${errexit_was_set}" -eq 1 ]; then
        set -e
    fi

    if [ "${rc}" -ne 0 ]; then
        echo "Error: ${launcher} launch failed. See ${log_file}." >&2
    fi

    return ${rc}
}
