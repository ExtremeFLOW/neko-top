#!/usr/bin/bash

# ============================================================================ #
# Define the help function

function help() {
    echo -e "run.sh case"
    echo -e "  Generate a box mesh for the unsteady rugby ball example."
    echo -e "  Lx is the length of the domain in the x direction."
    echo -e "  N sets the base resolution: Nx = N*Lx, Ny = N, Nz = 1."
    echo -e "  The z-extent is chosen so elements are approximately cubic."
    echo -e ""
    echo -e "  If no input arguments are provided, the default N is 10."
    echo -e ""
    echo -e "  Example usage:"
    echo -e "    run.sh -N10 -L2"
    echo -e ""
    echo -e " Options:"
    echo -e "  -h, --help  Show this help message and exit."
    echo -e "  -N#         Base number of cells (Ny = N, Nx = N*Lx)."
    echo -e "  -L#         Domain length in x (sets Nx = N*Lx)."
    echo -e ""
    echo -e "  See Readme for additional details."
    exit 0
}

# Handle options
N=10
Lx=2
for arg in "$@"; do
    if [ "${arg:0:2}" == "--" ]; then
        case ${arg:2} in
        help) help ;;
        Lx=*) Lx=${arg#--Lx=} ;;
        *) echo -e "Invalid option: $arg" >&2 && help ;;
        esac
    elif [ "${arg:0:1}" == "-" ]; then
        case ${arg:1:1} in
        h) help ;;
        N) N=${arg:2} ;;
        L) Lx=${arg:2} ;;
        *) echo -e "Invalid option: ${arg:1}" >&2 && help ;;
        esac
    fi
done
Nx=$((N*Lx))
Ny=$N
Nz=1

# ============================================================================ #
# Ensure Neko can be found and set default mesh size

if [ -z "${NEKO_DIR:-}" ]; then
    MAIN_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
    if [ -d "${MAIN_DIR}/external/neko/bin" ]; then
        export NEKO_DIR="${MAIN_DIR}/external/neko"
    fi
fi

if [ -n "${MAIN_DIR:-}" ] && [ -f "${MAIN_DIR}/build/mpmd_runtime.env" ]; then
    # shellcheck disable=SC1090
    source "${MAIN_DIR}/build/mpmd_runtime.env"
fi

if [ -n "${NEKO_DIR:-}" ]; then
    PATH="${NEKO_DIR}/bin:${PATH}"
fi

if [ -n "${MAIN_DIR:-}" ]; then
    for dep in neko hdf5 json-fortran adios2; do
        dep_lib="${MAIN_DIR}/external/${dep}/lib"
        if [ -d "${dep_lib}" ]; then
            LD_LIBRARY_PATH="${dep_lib}:${LD_LIBRARY_PATH:-}"
        fi
    done
    export LD_LIBRARY_PATH
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
# The aspect ratio each element is  1:1:1
z=$(python3 -c "print(1.0 / $N)")

echo "Generating mesh with dimensions: $Nx $Ny $Nz"
genmeshbox 0 $Lx 0 1 0 $z $Nx $Ny $Nz .false. .true. .true.

# End of file
# ============================================================================ #
