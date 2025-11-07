## MESHED meshed

gmsh -0 newton.geo
# This is to rebuild the mesh
gmsh2nek <<EOF
2
newton
0
EOF
#
rea2nbin newton.re2 newton.nmsh
