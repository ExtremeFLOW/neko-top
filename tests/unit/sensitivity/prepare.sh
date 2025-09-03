#!/usr/bin/bash

# ============================================================================ #
# Define the help function

function help() {
    echo -e "run.sh case"
    echo -e "  Generate a mesh."
    echo -e "  The input arguments are the number of cells in the x, y, and z"
    echo -e "  directions, respectively."
    echo -e ""
    echo -e "  If no input arguments are provided, the default mesh size is"
    echo -e "  32x8x8."
    echo -e ""
    echo -e "  Example usage:"
    echo -e "    run.sh -x32 -y8 -z8"
    echo -e ""
    echo -e " Options:"
    echo -e "  -h, --help  Show this help message and exit."
    echo -e "  -N#         Number of cells in the x and y direction."
    echo -e ""
    echo -e "  See Readme for additional details."
    exit 0
}

# Handle options
N=10
for arg in "$@"; do
    if [ "${arg:0:2}" == "--" ]; then
        case ${arg:2} in
        help) help ;;
        *) echo -e "Invalid option: $arg" >&2 && help ;;
        esac
    elif [ "${arg:0:1}" == "-" ]; then
        case ${arg:1:1} in
        h) help ;;
        N) N=${arg:2} ;;
        *) echo -e "Invalid option: ${arg:1}" >&2 && help ;;
        esac
    fi
done
Nx=$((N * 2))
Ny=$N
Nz=1

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
# The aspect ratio each element is  1:1:1
z=$(python3 -c "print(1.0 / $N)")

echo "Generating mesh with dimensions: $Nx $Ny $Nz"
genmeshbox 0 2 0 1 0 $z $Nx $Ny $Nz .false. .true. .true.

# End of file
# ============================================================================ #
