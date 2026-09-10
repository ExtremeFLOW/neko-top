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
    exit 0
}

CURRENT_DIR=$(pwd)
WORKING_DIR=$(dirname "$0")
cd "$WORKING_DIR" || exit 1

# Handle options
N=8
NP=4
for arg in "$@"; do
    if [ "${arg:0:2}" == "--" ]; then
        case ${arg:2} in
        help) help ;;
        *) echo -e "Invalid option: $arg" >&2 && help ;;
        esac
    elif [ "${arg:0:1}" == "-" ]; then
        case ${arg:1:1} in
        h) help ;;
        N) N=${arg:2} ;;                # e.g. -N8
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

Nx=$((N * 2))
Ny=$N
Nz=1

# Clean up old files
find . -maxdepth 1 -type f -name "steady_state_data*.csv" -delete
find . -maxdepth 1 -type f -name "FD_check_*.csv" -delete
find . -maxdepth 1 -type f -name "*.log" -delete

# ============================================================================ #
# Ensure Neko can be found and set default mesh size

if [ "$NEKO_DIR" ]; then
    PATH=$NEKO_DIR/bin:$PATH
fi

if [[ -z $(which genmeshbox) ]]; then
    echo -e "Neko tool 'genmeshbox' not found." >&2
    echo -e "Please ensure Neko is installed and in your PATH." >&2
    echo -e "Alternatively, set the NEKO_DIR environment variable." >&2
    exit 1
fi

# ============================================================================ #
# Generate mesh and run case

# Compute the depth to retain the aspect ratio of the elements
# The aspect ratio each element is 1:1:1
z=$(python3 -c "print(1.0 / $N)")

echo "Generating mesh with dimensions: $Nx $Ny $Nz"
genmeshbox 0 2 0 1 0 $z $Nx $Ny $Nz .false. .true. .true.

cases=$(find cases/ -type f -name "*.case")
for case in ${cases[@]}; do
    case=${case#./}
    case_name=$(basename "$case")
    export NEKO_LOG_FILE="neko_${case_name%.*}.log"
    export NEKO_TOP_LOG_FILE="nekotop_${case_name%.*}.log"
    echo "Running ${case} with mpirun -n $NP"
    mpirun -n "$NP" ./sensitivity_regression_driver "${case}" || exit 1
    mv steady_state_data.csv steady_state_data_${case_name%.*}.csv
done

python steady_state_plotter.py || exit 1
python FD_check.py || exit 1

# Clean up generated files
find . -maxdepth 1 -type f -name "box.nmsh" -delete
find . -maxdepth 1 -type f -name "*.log" -delete
find . -maxdepth 1 -type f -name "*.nek5000" \
    -exec sh -c 'rm "$1" "${1%.nek5000}".f*' _ {} \;

cd "$CURRENT_DIR" || exit 1

# ============================================================================ #
# End of file
# ============================================================================ #
