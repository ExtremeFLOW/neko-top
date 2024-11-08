#!/usr/bin/bash

# Tim, I am in AWE with how beutiful your run scripts are,
# one day I will learn to write them so cleanly... that day is not today...
# For an explanation of the mesh parameters, take a look at:
# https://media.enccs.se/2021/08/Nek5000_hands_on_1.pdf

# This would have been for periodic BCs
# But I think we're better with 'on' conditions

### Immersed meshed
#gmsh -0 immersed_M1.geo
## This is to rebuild the mesh
#gmsh2nek << EOF
#2
#brink_cyl
#1
#4 3
#0 30.0 0
#EOF
##
#rea2nbin brink_cyl.re2 immersed_M1.nmsh

[ -z "$NEKO_DIR" ] && NEKO_DIR=../../external/neko

if [ ! -f $NEKO_DIR/contrib/gmsh2nek/gmsh2nek ]; then
    echo "gmsh2nek not found in $NEKO_DIR/contrib/gmsh2nek/gmsh2nek"
    echo "Please set the NEKO_DIR environment variable to the path of the Neko installation"
    exit 1
fi

if [ ! -f $NEKO_DIR/bin/rea2nbin ]; then
    echo "rea2nbin not found in $NEKO_DIR/bin/rea2nbin"
    echo "Please set the NEKO_DIR environment variable to the path of the Neko installation"
    exit 1
fi

export PATH=$NEKO_DIR/bin:$NEKO_DIR/contrib/gmsh2nek:$PATH

## Immersed meshed
gmsh -0 immersed_M1.geo
# This is to rebuild the mesh
gmsh2nek <<EOF
2
brink_cyl
1
4 3
0 30.0 0
EOF
#
rea2nbin brink_cyl.re2 immersed_M1.nmsh

gmsh -0 immersed_M2.geo
# This is to rebuild the mesh
gmsh2nek <<EOF
2
brink_cyl
1
4 3
0 30.0 0
EOF
#
rea2nbin brink_cyl.re2 immersed_M2.nmsh

gmsh -0 immersed_M3.geo
# This is to rebuild the mesh
gmsh2nek <<EOF
2
brink_cyl
1
4 3
0 30.0 0
EOF
#
rea2nbin brink_cyl.re2 immersed_M3.nmsh

gmsh -0 immersed_M4.geo
# This is to rebuild the mesh
gmsh2nek <<EOF
2
brink_cyl
1
4 3
0 30.0 0
EOF
#
rea2nbin brink_cyl.re2 immersed_M4.nmsh

## MESHED meshed

gmsh -0 meshed_M2.geo
# This is to rebuild the mesh
gmsh2nek <<EOF
2
brink_cyl
1
5 4
0 30.0 0
EOF
#
rea2nbin brink_cyl.re2 meshed_M2.nmsh

gmsh -0 meshed_M3.geo
# This is to rebuild the mesh
gmsh2nek <<EOF
2
brink_cyl
1
5 4
0 30.0 0
EOF
#
rea2nbin brink_cyl.re2 meshed_M3.nmsh

gmsh -0 meshed_M4.geo
# This is to rebuild the mesh
gmsh2nek <<EOF
2
brink_cyl
1
5 4
0 30.0 0
EOF
#
rea2nbin brink_cyl.re2 meshed_M4.nmsh

mv *.nmsh ../data_local/brinkman_parameters/
rm *.re2 *.geo_unrolled *.msh
