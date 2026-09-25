#!/bin/bash
set -euo pipefail

function help() {
    cat <<'EOF'
Generate the LUMI-G mpmd.conf/select_gpu pair used by low_Re.

Usage:
  generate_mpmd_conf.sh --case-file PATH [options]

Options:
  --nodes N                  Number of nodes. Defaults to $SLURM_NNODES or 1.
  --case-file PATH           Case file passed to Neko and the Python driver.
  --neko-exe PATH            Neko executable to place in mpmd.conf.
  --python-bin PATH          Python executable to place in mpmd.conf.
  --python-script PATH       Python POD driver script.
  --neko-ranks-per-node N    Default: 8.
  --python-ranks-per-node N  Default: 48.
  --output PATH              Output mpmd.conf path. Default: ./mpmd.conf.
  --select-gpu PATH          Output select_gpu path. Default: ./select_gpu.
  -h, --help                 Show this help message and exit.
EOF
}

repo_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../../../.." && pwd)

nodes=${SLURM_NNODES:-1}
case_file=""
neko_exe=${NEKO_EXE:-./neko}
python_bin=${PYTHON_BIN:-$(command -v python3 || command -v python || echo python3)}
python_script=${PYTHON_SCRIPT:-"${repo_root}/scripts/python/pod_state_recover.py"}
neko_ranks_per_node=${NEKO_RANKS_PER_NODE:-8}
python_ranks_per_node=${PY_RANKS_PER_NODE:-48}
output_path=./mpmd.conf
select_gpu_path=./select_gpu

while (($# > 0)); do
    case "$1" in
    --nodes)
        nodes=$2
        shift 2
        ;;
    --case-file)
        case_file=$2
        shift 2
        ;;
    --neko-exe)
        neko_exe=$2
        shift 2
        ;;
    --python-bin)
        python_bin=$2
        shift 2
        ;;
    --python-script)
        python_script=$2
        shift 2
        ;;
    --neko-ranks-per-node)
        neko_ranks_per_node=$2
        shift 2
        ;;
    --python-ranks-per-node)
        python_ranks_per_node=$2
        shift 2
        ;;
    --output)
        output_path=$2
        shift 2
        ;;
    --select-gpu)
        select_gpu_path=$2
        shift 2
        ;;
    -h|--help)
        help
        exit 0
        ;;
    *)
        echo "Error: unknown option: $1" >&2
        help >&2
        exit 1
        ;;
    esac
done

if [ -z "${case_file}" ]; then
    echo "Error: --case-file is required." >&2
    exit 1
fi

if [ "${nodes}" -lt 1 ]; then
    echo "Error: --nodes must be at least 1." >&2
    exit 1
fi

if [ "${neko_ranks_per_node}" -lt 1 ]; then
    echo "Error: --neko-ranks-per-node must be at least 1." >&2
    exit 1
fi

if [ "${python_ranks_per_node}" -lt 1 ]; then
    echo "Error: --python-ranks-per-node must be at least 1." >&2
    exit 1
fi

tasks_per_node=$((neko_ranks_per_node + python_ranks_per_node))
if [ "${tasks_per_node}" -ne 56 ]; then
    echo "Error: this LUMI-G helper expects 56 tasks per node." >&2
    echo "Got ${neko_ranks_per_node} Neko ranks + ${python_ranks_per_node} Python ranks." >&2
    exit 1
fi

cat <<'EOF' > "${select_gpu_path}"
#!/bin/bash

export ROCR_VISIBLE_DEVICES=${SLURM_LOCALID}
export NEKO_GS_COMM=${NEKO_GS_COMM:-MPI}
sleep "${NEKO_STARTUP_DELAY:-20}"
exec "$@"
EOF

chmod +x "${select_gpu_path}"

: > "${output_path}"

for ((node=0; node<nodes; node++)); do
    node_base=$((node * tasks_per_node))

    for ((local_rank=0; local_rank<neko_ranks_per_node; local_rank++)); do
        rank=$((node_base + local_rank))
        printf '%d /usr/bin/env NEKO_COMM_ID=0 NEKO_CTRL_PEER_ROOT=%d %s %s %s\n' \
            "${rank}" "${neko_ranks_per_node}" \
            "${select_gpu_path}" "${neko_exe}" "${case_file}" >> "${output_path}"
    done

    for ((local_rank=0; local_rank<python_ranks_per_node; local_rank++)); do
        rank=$((node_base + neko_ranks_per_node + local_rank))
        printf '%d /usr/bin/env NEKO_COMM_ID=1 NEKO_CTRL_PEER_ROOT=0 %s %s %s\n' \
            "${rank}" "${python_bin}" "${python_script}" "${case_file}" \
            >> "${output_path}"
    done
done

echo "Generated ${output_path} for ${nodes} LUMI-G node(s)"
echo "Neko ranks/node:   ${neko_ranks_per_node}"
echo "Python ranks/node: ${python_ranks_per_node}"
echo "Case file:         ${case_file}"
echo "Neko executable:   ${neko_exe}"
echo "Python command:    ${python_bin} ${python_script}"
