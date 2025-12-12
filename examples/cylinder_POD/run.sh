# Resolve ADIOS2_PATH when launched via the top-level run.sh (env already has
# ADIOS2_DIR/MAIN_DIR from scripts/dependencies.sh). Prefer an explicit
# ADIOS2_PATH, then ADIOS2_DIR, then the local external checkout.
if [ -z "$ADIOS2_PATH" ]; then
    if [ -n "$ADIOS2_DIR" ]; then
        ADIOS2_PATH="$ADIOS2_DIR"
    elif [ -n "$MAIN_DIR" ] && [ -d "$MAIN_DIR/external/adios2" ]; then
        ADIOS2_PATH="$MAIN_DIR/external/adios2"
    fi
fi

if [ -n "$ADIOS2_PATH" ]; then
    pyver=$(python3 -c 'import sys; print(f\"{sys.version_info.major}.{sys.version_info.minor}\")')
    export PYTHONPATH="${PYTHONPATH:+$PYTHONPATH:}$ADIOS2_PATH/lib/python${pyver}/site-packages"
    export LD_LIBRARY_PATH="${ADIOS2_PATH}/lib:${ADIOS2_PATH}/lib64:${LD_LIBRARY_PATH}"
    echo "Using ADIOS2_PATH=${ADIOS2_PATH}"
else
    echo "Warning: ADIOS2_PATH is not set; in-situ Python may fail to import ADIOS2." >&2
fi

gmsh -0 newton.geo
# This is to rebuild the mesh
gmsh2nek <<EOF
2
newton
0
EOF
#
rea2nbin newton.re2 ext_cyl.nmsh


mpirun -n 10 ./neko cylinder_POD.case : -n 4 python3 insitu_task.py
