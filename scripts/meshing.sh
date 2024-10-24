# ============================================================================ #
# Functions for meshing using Cubit, Exodus, and Nek5000
# ===========================================================================
#
# This script contains functions to create and convert meshes using Cubit,
# Exodus, and Nek5000. The functions are used in the `mesh.sh` script to
# automate the meshing process.

# ============================================================================ #
# Cubit Journal to neko binary mesh

function jou2nbin() {
    set -e

    find_cubit
    find_exo2nek
    find_rea2nbin

    input_file=$1
    input_name=$(basename ${1%.*})

    # Construct the mesh using Cubit
    $cubit -nojournal -nographics -batch -noecho $1

    if [[ $(find ./ -name "*.exo" | wc -l) -lt 1 ]]; then
        printf >&2 "\n\e[4mError:\e[0m\n"
        printf >&2 "  %-10s %-67s\n" "Cubit:" "Exodus not created: $input_file"
        exit 1
    fi

    # Convert the mesh from Exodus format to Nek5000 format
    (
        echo "$(ls *.exo | wc -l)"
        for file in *.exo; do echo "${file%.*}"; done
        echo "0"
        echo "0"
        echo "$input_name"
    ) | $exo2nek

    if [[ $(find ./ -name "*.re2" | wc -l) -lt 1 ]]; then
        printf >&2 "\n\e[4mError:\e[0m\n"
        printf >&2 "  %-10s %-67s\n" "exo2nek:" "re2 not created: $input_file"
        exit 1
    fi

    # Convert the mesh to Neko mesh format
    $rea2nbin $input_name.re2 $input_name.nmsh

    if [[ $(find ./ -name "*.nmsh" | wc -l) -lt 1 ]]; then
        printf >&2 "\n\e[4mError:\e[0m\n"
        printf >&2 "  %-10s %-67s\n" "rea2nbin:" "re2 not created: $input_file"
        exit 1
    fi
}

# ============================================================================ #
# Exodus to neko binary mesh

function exo2nbin() {
    set -e

    find_exo2nek
    find_rea2nbin

    input_file=$1
    input_name=$(basename ${1%.*})

    ln -sf $input_file $input_name.exo

    # Convert the mesh from Exodus format to Nek5000 format
    (
        echo "1"
        echo "$input_name"
        echo "0"
        echo "0"
        echo "$input_name"
    ) | $exo2nek

    if [[ $(find ./ -name "*.re2" | wc -l) -lt 1 ]]; then
        printf >&2 "\n\e[4mError:\e[0m\n"
        printf >&2 "  %-10s %-67s\n" "exo2nek:" "re2 not created: $input_file"
        exit 1
    fi

    # Convert the mesh to Neko mesh format
    $rea2nbin $input_name.re2 $input_name.nmsh

    if [[ $(find ./ -name "*.nmsh" | wc -l) -lt 1 ]]; then
        printf >&2 "\n\e[4mError:\e[0m\n"
        printf >&2 "  %-10s %-67s\n" "rea2nbin:" "re2 not created: $input_file"
        exit 1
    fi
}

# ============================================================================ #
# Nek5000 mesh to neko binary mesh

function re2nbin() {
    set -e

    find_rea2nbin

    input_file=$1
    input_name=$(basename ${1%.*})

    ln -sf $input_file $input_name.re2

    # Convert the mesh to Neko mesh format
    $rea2nbin $input_name.re2 $input_name.nmsh

    if [[ $(find ./ -name "*.nmsh" | wc -l) -lt 1 ]]; then
        printf >&2 "\n\e[4mError:\e[0m\n"
        printf >&2 "  %-10s %-67s\n" "rea2nbin:" "re2 not created: $input_file"
        exit 1
    fi
}

# ============================================================================ #
# GMsh to neko binary mesh

function geo2nbin() {
    set -e

    find_gmsh
    find_gmsh2nek
    find_rea2nbin

    input_file=$1
    input_name=$(basename ${1%.*})

    ln -sf $input_file $input_name.geo

    gmsh -0 $input_name.geo

    # Grep the name of the mesh file (Save "FILE_NAME.msh";)
    mesh_file=$(grep "Save \"*\"" $input_name.geo | cut -d\" -f2)

    # Convert the mesh from GMsh format to Nek5000 format
    (
        echo "$DIMENSION"
        echo "${mesh_file%.*}"
        echo "0"
    ) | $gmsh2nek

    if [[ $(find ./ -name "*.re2" | wc -l) -lt 1 ]]; then
        printf >&2 "\n\e[4mError:\e[0m\n"
        printf >&2 "  %-10s %-67s\n" "gmsh2nek:" "re2 not created: $input_file"
        exit 1
    fi

    # Convert the mesh to Neko mesh format
    $rea2nbin ${mesh_file%.*}.re2 $input_name.nmsh

    if [[ $(find ./ -name "*.nmsh" | wc -l) -lt 1 ]]; then
        printf >&2 "\n\e[4mError:\e[0m\n"
        printf >&2 "  %-10s %-67s\n" "rea2nbin:" "re2 not created: $input_file"
        exit 1
    fi
}
