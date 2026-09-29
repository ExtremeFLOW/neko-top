#!/bin/bash

function help() {
    echo -e "./prepare.sh"
    echo -e "  Generate the mesh used by the low-Re POD example."
    echo -e ""
    echo -e "  The input arguments are the number of cells in the x, y, and z"
    echo -e "  directions, respectively."
    echo -e ""
    echo -e "  If no input arguments are provided, the default mesh size is"
    echo -e "  48x16x16."
    echo -e ""
    echo -e "  Example usage:"
    echo -e "    ./prepare.sh -x48 -y16 -z16"
    echo -e ""
    echo -e " Options:"
    echo -e "  -h, --help  Show this help message and exit."
    echo -e "  -x#         Number of cells in the x direction."
    echo -e "  -y#         Number of cells in the y direction."
    echo -e "  -z#         Number of cells in the z direction."
    exit 0
}

Nx=48 && Ny=16 && Nz=16
for arg in "$@"; do
    if [ "${arg:0:2}" == "--" ]; then
        case ${arg:2} in
        help) help ;;
        *) echo -e "Invalid option: $arg" >&2 && help ;;
        esac
    elif [ "${arg:0:1}" == "-" ]; then
        case ${arg:1:1} in
        h) help ;;
        x) Nx=${arg:2} ;;
        y) Ny=${arg:2} ;;
        z) Nz=${arg:2} ;;
        *) echo -e "Invalid option: ${arg:1}" >&2 && help ;;
        esac
    fi
done

if [ "${NEKO_DIR:-}" ]; then
    PATH=$NEKO_DIR/bin:$PATH
fi

if [[ -z $(which genmeshbox) ]]; then
    echo -e "Neko tool 'genmeshbox' not found." >&2
    echo -e "Please ensure Neko is installed and in your PATH." >&2
    echo -e "Alternatively, set the NEKO_DIR environment variable." >&2
    exit 1
fi

echo "Generating mesh with dimensions: $Nx $Ny $Nz"
genmeshbox 0 6 0 2 0 2 $Nx $Ny $Nz .false. .false. .false.
