#!/usr/bin/env bash
# Run one rugby-ball case and check its optimization history.
#
# Usage:
#   run_rugby_ball_check.sh <binary> <case_file> <mesh_n> <ranks> \
#       <check command ...>
#
# The case runs in the current directory on an <mesh_n> x <mesh_n> mesh
# from prepare.sh. The check command is then called with the path of the
# optimization_data.csv the run wrote appended as its last argument.
#
# NEKO_TOP_RUGBY_BALL_RANKS, when set, overrides <ranks>.
set -euo pipefail

if [[ $# -lt 5 ]]; then
    echo "Usage: run_rugby_ball_check.sh <binary> <case_file> <mesh_n>" \
        "<ranks> <check command ...>" >&2
    exit 2
fi

binary="$1"
case_file="$2"
mesh_n="$3"
ranks="${NEKO_TOP_RUGBY_BALL_RANKS:-$4}"
shift 4
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ ! -x "$binary" ]]; then
    echo "Executable not found or not executable: $binary" >&2
    exit 2
fi
if [[ ! -f "$case_file" ]]; then
    echo "Required case file not found: $case_file" >&2
    exit 2
fi

case_name="$(basename "$case_file")"
cp -f "$case_file" "$case_name"

# The mesh depends only on mesh_n, so regenerate it when that changes.
if [[ ! -f box.nmsh || "$(cat mesh_n.txt 2>/dev/null)" != "$mesh_n" ]]; then
    rm -f box.nmsh
    bash "${script_dir}/prepare.sh" -N"$mesh_n"
    echo "$mesh_n" > mesh_n.txt
fi

# The optimizer appends to an existing optimization_data.csv, so remove
# the one an earlier run left behind: the check must only ever see the
# output of this run.
rm -f optimization_data.csv steady_state_data.csv

export NEKO_LOG_FILE="neko_${case_name%.*}.log"
echo "Running ${case_name} on ${ranks} ranks, ${mesh_n}x${mesh_n} mesh"

# A non-zero exit fails the test. That includes the optimizer's final
# check that every constraint is satisfied, so the case must end on a
# feasible design.
mpirun -n "$ranks" "$binary" "$case_name"

"$@" optimization_data.csv
