#!/usr/bin/bash

# ============================================================================ #
# Define the help function

function help() {
    echo -e "run.sh case"
    echo -e "  Generate a mesh of size [0,2]x[0,1]."
    echo -e "  The input arguments are the number of cells in y direction and"
    echo -e "  the number of processors to run with"
    echo -e ""
    echo -e "  If no input arguments are provided, the default mesh size is 10"
    echo -e "  and run with 4 processors"
    echo -e ""
    echo -e "  Example usage:"
    echo -e "    run.sh -N32 -np=4"
    echo -e ""
    echo -e " Options:"
    echo -e "  -h, --help  Show this help message and exit."
    echo -e "  -N#         Number of cells in the x and y direction (e.g. -N8)."
    echo -e "  -np=#       MPI ranks for mpirun -n (e.g. -np=4)."
    echo -e ""
    echo -e " Environment:"
    echo -e "  NEKO_TOP_BUILD_DIR  Build directory containing reg_mma_bin"
    echo -e "                      (default: <repo root>/build)."
    echo -e ""
    exit 0
}

CURRENT_DIR=$(pwd)
WORKING_DIR=$(cd "$(dirname "$0")" && pwd)
ROOT_DIR=$(realpath "${WORKING_DIR}/../../../")

# The binary is built into the build tree (its RUNTIME_OUTPUT_DIRECTORY),
# not next to this script, so it is not simply "./reg_mma_bin" any more.
# Default to the conventional build directory created by `./setup.sh -t`,
# overridable for a differently-named build directory. Resolved here, before
# the cd below, so a relative NEKO_TOP_BUILD_DIR is taken relative to the
# caller's own directory rather than silently reinterpreted relative to this
# script's directory.
BUILD_DIR=$(realpath -m "${NEKO_TOP_BUILD_DIR:-${ROOT_DIR}/build}")
REG_MMA_BIN="${BUILD_DIR}/tests/regression/mma/reg_mma_bin"

cd "$WORKING_DIR" || exit 1

# Handle options
NP=1
for arg in "$@"; do
    if [ "${arg:0:2}" == "--" ]; then
        case ${arg:2} in
        help) help ;;
        *) echo -e "Invalid option: $arg" >&2 && help ;;
        esac
    elif [ "${arg:0:1}" == "-" ]; then
        case ${arg:1:1} in
        h) help ;;
        n)                             # handle -np=4
           if [ "${arg:0:4}" == "-np=" ]; then
              NP=${arg:4}
           else
              echo -e "Invalid option: $arg" >&2 && help
           fi
           ;;
        *) echo -e "Invalid option: ${arg:1}" >&2 && help ;;
        esac
    fi
done

if [[ ! -x "$REG_MMA_BIN" ]]; then
    echo -e "reg_mma_bin not found or not executable: $REG_MMA_BIN" >&2
    echo -e "Build Neko-TOP's tests first (e.g. ./setup.sh -t), or set" >&2
    echo -e "NEKO_TOP_BUILD_DIR to the build directory used." >&2
    exit 1
fi

# Clean up old files
find . -maxdepth 1 -type f -name "*.csv" -delete
find . -maxdepth 1 -type f -name "*.log" -delete

# ============================================================================ #
# Generate mesh and run case

cases=$(find cases/ -type f -name "*.case")
for case in ${cases[@]}; do
    case=${case#./}
    case_name=$(basename "$case" | cut -f1 -d'.')
    export NEKO_LOG_FILE="neko_${case_name}.log"
    echo "Running ${case} with mpirun -n $NP"
    mpirun -n "$NP" "$REG_MMA_BIN" "${case}" || exit 1
    mv optimization_data.csv optimization_data_${case_name}.csv
done

if [[ -z "$VIRTUAL_ENV" && -f "${ROOT_DIR}/.venv/bin/activate" ]]; then
    source "${ROOT_DIR}/.venv/bin/activate"
fi

python plot_design.py all || exit 1
python check.py || exit 1

echo "MMA regression test passed"
cd "$CURRENT_DIR" || exit 1

# ============================================================================ #
# End of file
# ============================================================================ #
