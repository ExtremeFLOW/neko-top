## MESHED meshed

gmsh -0 ext_cyl.geo
# This is to rebuild the mesh
gmsh2nek <<EOF
2
brink_cyl
0
EOF
#
rea2nbin brink_cyl.re2 ext_cyl.nmsh
