## MESHED meshed

gmsh -0 newton.geo
# This is to rebuild the mesh
gmsh2nek <<EOF
3
newton
1
6 7
0 0 4
EOF
#
rea2nbin newton.re2 newton.nmsh
python3 extrude_baseflow.py baseflow0.f00000 baseflow3d0.f00000 --length 4 --layers 4
