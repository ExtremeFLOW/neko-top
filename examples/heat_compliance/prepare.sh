#!/usr/bin/env bash

# ============================================================================ #
# Define the help function

function help() {
    echo -e "run.sh [-N#]"
    echo -e "  Generate a Gmsh mesh for a 3x3x3 cube."
    echo -e ""
    echo -e "  N is the number of cells per coarse square in x and y."
    echo -e "  The bottom face is split into 3x3 coarse squares; each"
    echo -e "  is further subdivided into N x N cells."
    echo -e ""
    echo -e "  Bottom center patch has boundary id 2."
    echo -e "  All other external boundaries have boundary id 1."
    echo -e ""
    echo -e "  Example usage:"
    echo -e "    run.sh -N2"
    echo -e ""
    echo -e " Options:"
    echo -e "  -h, --help  Show this help message and exit."
    echo -e "  -N#         Number of cells per coarse square in x,y."
    echo -e ""
    exit 0
}

# Handle options
N=2
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

# ============================================================================ #
# Generate mesh with Gmsh

GEO_FILE="cube_mesh.geo"

echo "Generating mesh with N = $N (bottom has (3N)x(3N) cells)..."
gmsh -0 "${GEO_FILE}" -setnumber N "${N}"

if [ $? -ne 0 ]; then
    echo "Gmsh failed to generate the mesh." >&2
    exit 1
fi

gmsh2nek <<EOF
3
cube_mesh
0
EOF
#
rea2nbin cube_mesh.re2 box.nmsh

echo "Mesh written to box.nmsh"

# Here you can add the command that converts/uses this mesh for Neko, e.g.:
# neko_mesh_converter "${MESH_FILE}" ...
# neko_case_runner ...

# End of file
# ============================================================================ #

