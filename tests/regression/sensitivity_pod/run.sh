#!/usr/bin/env bash
set -euo pipefail

function help() {
    echo "run.sh [-N#] [-np=#] [-np_py=#]"
    echo "  Builds a mesh and runs sensitivity with checkpointing and POD sweeps."
    echo ""
    echo "Options:"
    echo "  -h, --help     Show this help message and exit."
    echo "  -N#            Number of cells in y (x uses 2*N). Default: 8"
    echo "  -np=#          MPI ranks for the Neko driver. Default: 4"
    echo "  -np_py=#       MPI ranks for the Python POD helper. Default: 1"
    echo ""
    echo "Environment overrides:"
    echo "  POD_MODES       Space-separated list, e.g. '1 2 3 4'."
    echo "  POD_BATCH_SIZE  POD batch size (default: 2)."
    echo "  POD_I_STREAM    Stream interval (default: 1)."
    echo "  POD_DTYPE       single|double (default: double)."
    echo "  POD_WRITE_MODES true|false (default: false)."
    echo "  PYTHON          Python executable (default: python3)."
    echo "  PY_STARTUP_DELAY Seconds to wait before launching Neko (default: 2)."
    exit 0
}

CURRENT_DIR=$(pwd)
WORKING_DIR=$(dirname "$0")
cd "$WORKING_DIR" || exit 1

# Defaults
N=8
NP=4
NP_PY=1

for arg in "$@"; do
    if [ "${arg:0:2}" == "--" ]; then
        case ${arg:2} in
            help) help ;;
            *) echo "Invalid option: $arg" >&2; help ;;
        esac
    elif [ "${arg:0:1}" == "-" ]; then
        case ${arg:1:1} in
            h) help ;;
            N) N=${arg:2} ;;
            n)
               if [ "${arg:0:4}" == "-np=" ]; then
                   NP=${arg:4}
               elif [ "${arg:0:7}" == "-np_py=" ]; then
                   NP_PY=${arg:7}
               else
                   echo "Invalid option: $arg" >&2
                   help
               fi
               ;;
            *) echo "Invalid option: ${arg:1}" >&2; help ;;
        esac
    fi
done

NEKO_RANKS=${NEKO_RANKS:-$NP}
PY_RANKS=${PY_RANKS:-$NP_PY}
POD_MODES=${POD_MODES:-"1 2 3 4 5"}
POD_BATCH_SIZE=${POD_BATCH_SIZE:-2}
POD_I_STREAM=${POD_I_STREAM:-1}
POD_DTYPE=${POD_DTYPE:-double}
POD_WRITE_MODES=${POD_WRITE_MODES:-false}
PYTHON=${PYTHON:-python3}
PY_STARTUP_DELAY=${PY_STARTUP_DELAY:-2}
RUN_TAG=${RUN_TAG:-$(date +%Y%m%d_%H%M%S)}

REPO_ROOT=$(cd "$WORKING_DIR/../../.." && pwd)
RESULTS_DIR="$WORKING_DIR/results/${RUN_TAG}"
mkdir -p "$RESULTS_DIR"

DRIVER="$WORKING_DIR/sensitivity_pod_regression_driver"
if [ ! -x "$DRIVER" ]; then
    echo "Driver not found: $DRIVER" >&2
    echo "Build with: cmake --build . --target Regression" >&2
    exit 1
fi

PY_SCRIPT="$WORKING_DIR/pod_state_recover.py"
if [ ! -f "$PY_SCRIPT" ]; then
    echo "Python script not found: $PY_SCRIPT" >&2
    exit 1
fi

# Ensure Neko tools are in PATH
if [ "${NEKO_DIR:-}" ]; then
    PATH="$NEKO_DIR/bin:$PATH"
fi

if [[ -z $(which genmeshbox) ]]; then
    echo "Neko tool 'genmeshbox' not found." >&2
    echo "Please ensure Neko is installed and in your PATH." >&2
    echo "Alternatively, set the NEKO_DIR environment variable." >&2
    exit 1
fi

# ADIOS2 runtime for Python (optional)
if [ -z "${ADIOS2_PATH:-}" ]; then
    if [ -n "${ADIOS2_DIR:-}" ]; then
        ADIOS2_PATH="$ADIOS2_DIR"
    elif [ -d "$REPO_ROOT/external/adios2" ]; then
        ADIOS2_PATH="$REPO_ROOT/external/adios2"
    fi
fi

if [ -n "${ADIOS2_PATH:-}" ]; then
    pyver=$($PYTHON -c 'import sys; print("{}.{}".format(sys.version_info.major, sys.version_info.minor))')
    export PYTHONPATH="${ADIOS2_PATH}/lib/python${pyver}/site-packages:${PYTHONPATH:-}"
    export LD_LIBRARY_PATH="${ADIOS2_PATH}/lib:${ADIOS2_PATH}/lib64:${LD_LIBRARY_PATH:-}"
    echo "Using ADIOS2_PATH=${ADIOS2_PATH}"
fi

echo "Results:        $RESULTS_DIR"
echo "Neko ranks:     $NEKO_RANKS"
echo "Python ranks:   $PY_RANKS"
echo "POD modes:      $POD_MODES"
echo "POD batch size: $POD_BATCH_SIZE"
echo "POD i_stream:   $POD_I_STREAM"
echo "POD dtype:      $POD_DTYPE"
echo "POD write:      $POD_WRITE_MODES"

Nx=$((N * 2))
Ny=$N
Nz=1
z=$($PYTHON -c "print(1.0 / $N)")

MESH_DIR="$RESULTS_DIR/mesh"
mkdir -p "$MESH_DIR"
(
    cd "$MESH_DIR"
    echo "Generating mesh with dimensions: $Nx $Ny $Nz"
    genmeshbox 0 2 0 1 0 "$z" "$Nx" "$Ny" "$Nz" .false. .true. .true.
)

function generate_case() {
    local base_case=$1
    local out_case=$2
    local mode=$3
    local n_modes=${4:-0}

    $PYTHON - <<PY
import json, re
base_case = "${base_case}"
out_case = "${out_case}"
mode = "${mode}"

with open(base_case, "r", encoding="utf-8") as fh:
    raw = fh.read()
cleaned = re.sub(r",\s*([}\]])", r"\1", raw)
case = json.loads(cleaned)

sr = case.get("state_recovery", {})
if mode == "checkpoint":
    sr["type"] = "checkpoint"
else:
    sr["type"] = "pod"
    sr["batch_size"] = int("${POD_BATCH_SIZE}")
    sr["n_modes"] = int(n_modes)
    sr["i_stream"] = int("${POD_I_STREAM}")
    sr["dtype"] = "${POD_DTYPE}"
    write_modes_str = "${POD_WRITE_MODES}".strip().lower()
    sr["write_modes"] = write_modes_str in ("1", "true", "yes", "y", "t")
    sr["enabled"] = True
case["state_recovery"] = sr

with open(out_case, "w", encoding="utf-8") as fh:
    json.dump(case, fh, indent=2)
PY
}

function run_checkpoint_case() {
    local case_file=$1
    local case_name=$2
    local run_dir=$3

    mkdir -p "$run_dir"
    cp "$MESH_DIR/box.nmsh" "$run_dir/box.nmsh"
    generate_case "$case_file" "$run_dir/${case_name}.case" "checkpoint"

    echo "Running checkpoint: ${case_name}"
    (
        cd "$run_dir"
        export NEKO_LOG_FILE="neko_${case_name}_checkpoint.log"
        mpirun -n "$NEKO_RANKS" "$DRIVER" "${case_name}.case" > neko.log 2>&1
    )
}

function run_pod_case() {
    local case_file=$1
    local case_name=$2
    local run_dir=$3
    local n_modes=$4

    mkdir -p "$run_dir"
    cp "$MESH_DIR/box.nmsh" "$run_dir/box.nmsh"
    generate_case "$case_file" "$run_dir/${case_name}.case" "pod" "$n_modes"

    echo "Running POD: ${case_name} (n_modes=${n_modes})"
    (
        cd "$run_dir"
        export NEKO_LOG_FILE="neko_${case_name}_pod_${n_modes}.log"
        mpirun --tag-output -n "$PY_RANKS" "$PYTHON" "$PY_SCRIPT" "${case_name}.case" > python.log 2>&1 &
        py_pid=$!
        sleep "$PY_STARTUP_DELAY"
        mpirun --tag-output -n "$NEKO_RANKS" "$DRIVER" "${case_name}.case" > neko.log 2>&1
        kill "$py_pid" 2>/dev/null || true
        wait "$py_pid" 2>/dev/null || true
    )
}

for case_file in cases/*.case; do
    case_base=$(basename "$case_file")
    case_name=${case_base%.case}
    case_dir="$RESULTS_DIR/${case_name}"

    run_checkpoint_case "$case_file" "$case_name" "$case_dir/checkpoint"

    for n_modes in $POD_MODES; do
        run_pod_case "$case_file" "$case_name" "$case_dir/pod_n${n_modes}" "$n_modes"
    done

done

cd "$CURRENT_DIR" || exit 1
