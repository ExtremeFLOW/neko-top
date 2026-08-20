#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 5 ]]; then
    echo "Usage: run_valgrind_check.sh <binary> <case_file> <supp_file> <prepare_script> <checker_script>" >&2
    exit 2
fi

binary="$1"
case_file="$2"
supp_file="$3"
prepare_script="$4"
checker_script="$5"
case_name="$(basename "$case_file")"

if ! command -v valgrind >/dev/null 2>&1; then
    echo "Skipping valgrind regression test: valgrind not found in PATH."
    exit 77
fi

if ! command -v python3 >/dev/null 2>&1; then
    echo "python3 is required for valgrind regression checks." >&2
    exit 2
fi

if [[ ! -x "$binary" ]]; then
    echo "Executable not found or not executable: $binary" >&2
    exit 2
fi

for file in "$case_file" "$supp_file" "$prepare_script" "$checker_script"; do
    if [[ ! -f "$file" ]]; then
        echo "Required file not found: $file" >&2
        exit 2
    fi
done

# Stage all required inputs in the current working directory.
cp -f "$case_file" "$case_name"
cp -f "$supp_file" valgrind_runtime.supp
cp -f "$prepare_script" prepare.sh
cp -f "$checker_script" check_valgrind.py

# Remove stale artifacts from previous runs in this build directory.
rm -f box.nmsh prepare.log valgrind_regression.log valgrind_stdout.log valgrind_stderr.log
rm -f optimization_data.csv neko_valgrind.log
rm -rf checkpoints

mesh_nx="${VALGRIND_MESH_NX:-16}"
mesh_ny="${VALGRIND_MESH_NY:-4}"
mesh_nz="${VALGRIND_MESH_NZ:-4}"

# Generate a deterministic mesh local to this test run.
bash ./prepare.sh -x"${mesh_nx}" -y"${mesh_ny}" -z"${mesh_nz}" > prepare.log

# Execute the case with valgrind and parse leak summary.
if ! valgrind \
    --leak-check=full \
    --show-leak-kinds=all \
    --undef-value-errors=no \
    --num-callers=20 \
    --suppressions=valgrind_runtime.supp \
    --log-file=valgrind_regression.log \
    "$binary" "$case_name" > valgrind_stdout.log 2> valgrind_stderr.log; then
    echo "Valgrind execution failed. Tail of stdout/stderr:" >&2
    tail -n 40 valgrind_stdout.log >&2 || true
    tail -n 40 valgrind_stderr.log >&2 || true
    exit 1
fi

python3 ./check_valgrind.py \
    --log valgrind_regression.log \
    --max-definite-bytes "${VALGRIND_MAX_DEFINITE_BYTES:-0}" \
    --max-indirect-bytes "${VALGRIND_MAX_INDIRECT_BYTES:-0}" \
    --max-possible-bytes "${VALGRIND_MAX_POSSIBLE_BYTES:-0}" \
    --max-reachable-bytes "${VALGRIND_MAX_REACHABLE_BYTES:-460000}"

echo "Valgrind regression test passed."
