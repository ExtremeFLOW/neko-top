#!/usr/bin/bash

# ============================================================================ #
# Generate the unit-square mesh for the Borrvall & Petersson (2003) Sec. 4.4
# "rugby ball" validation case.
#
# Domain: unit square [0,1]x[0,1], NON-periodic in x and y (the four walls
# all carry a prescribed-velocity Dirichlet BC, per the paper), periodic in
# the thin z direction (one element deep, unit aspect ratio, exactly as
# `examples/rugby_ball/prepare.sh` does it for its own box).
#
# genmeshbox zone numbering (see $NEKO_DIR/contrib/genmeshbox/genmeshbox.f90):
#   zone 1 = x=x0 face (left),   zone 2 = x=x1 face (right)
#   zone 3 = y=y0 face (bottom), zone 4 = y=y1 face (top)
#   zone 5 = z=z0 face,          zone 6 = z=z1 face
# When a direction is requested periodic, genmeshbox does not label that
# pair of facets at all (it calls mark_periodic_facet instead of
# mark_labeled_facet), so with x and y non-periodic here, zones 1-4 are the
# four walls of the unit square and are exactly what the case file's
# "zone_indices": [1, 2, 3, 4] boundary condition refers to.
# ============================================================================ #

function help() {
    echo -e "prepare.sh [-N#]"
    echo -e "  Generate the unit-square mesh for the Borrvall & Petersson"
    echo -e "  rugby-ball validation case."
    echo -e ""
    echo -e "  -N#  Number of elements in the x and y direction (both equal)."
    echo -e "       Default: 20. The paper itself used 100x100 U-elements,"
    echo -e "       which is Q1-equivalent; with Neko's much higher-order SEM"
    echo -e "       basis a far coarser N is the expected starting point for"
    echo -e "       a mesh-convergence study."
    echo -e "       Use a multiple of 10: the design domain's edges at 0.05 and"
    echo -e "       0.95 then fall where the GLL-node mask gives an area of"
    echo -e "       exactly 0.81 (measured at polynomial order 5). Other N"
    echo -e "       (e.g. 14) shift the area, and with it the solid amount the"
    echo -e "       volume limit allows, by a few percent."
    echo -e ""
    echo -e " Options:"
    echo -e "  -h, --help  Show this help message and exit."
    exit 0
}

# Default mesh. Keep N a multiple of 10 so that the design-domain mask covers
# exactly the [0.05, 0.95]^2 area (see help).
N=20
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
Nx=$N && Ny=$N && Nz=1

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
# Generate mesh and run case

# Compute the depth to retain the aspect ratio of each element as 1:1:1.
z=$(python3 -c "print(1.0 / $N)")

echo "Generating mesh: $Nx $Ny $Nz (non-periodic x,y; periodic z)"
genmeshbox 0 1 0 1 0 $z $Nx $Ny $Nz .false. .false. .true.

# End of file
# ============================================================================ #
