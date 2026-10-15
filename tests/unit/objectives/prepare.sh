#!/usr/bin/bash

# ============================================================================ #
# Generate the box mesh used by the objective time-window tests.
#
# A short channel, deliberately coarse: these tests check that an objective's
# reported value does not depend on the time window it is accumulated over,
# which is a property of the accumulation rather than of the resolution.
# ============================================================================ #

function help() {
    echo -e "prepare.sh"
    echo -e "  Generate the mesh for the objective time-window tests."
    echo -e ""
    echo -e " Options:"
    echo -e "  -h, --help  Show this help message and exit."
    echo -e "  -N#         Number of cells across the channel. Default 4."
    exit 0
}

N=4
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
Nz=$N

# ============================================================================ #
# Ensure Neko can be found

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
# Generate the mesh. Cubic elements, no periodicity, so the six boundary zones
# are the six faces of the box.

echo "Generating mesh with dimensions: $Nx $Ny $Nz"
genmeshbox 0 2 0 1 0 1 $Nx $Ny $Nz .false. .false. .false.

# End of file
# ============================================================================ #
