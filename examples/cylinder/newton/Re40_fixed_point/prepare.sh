#!/usr/bin/env bash
set -euo pipefail

gmsh -0 newton.geo
gmsh2nek <<EOF
2
newton
0
EOF
rea2nbin newton.re2 newton.nmsh
