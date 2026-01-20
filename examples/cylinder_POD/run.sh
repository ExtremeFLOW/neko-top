#!/bin/bash

# Ensure 'conda' is available
if [ -f "/scratch/nobis/miniconda3/etc/profile.d/conda.sh" ]; then
    . "/scratch/nobis/miniconda3/etc/profile.d/conda.sh"
else
    eval "$(/scratch/nobis/miniconda3/bin/conda shell.bash hook)"
fi

# Activate the pySEMTools env (this is where mpi4py was rebuilt)
conda activate /scratch/nobis/POSTDOC/PYSEMTOOLS/miniconda3/.conda/envs/pySEMTools

# Resolve ADIOS2_PATH
if [ -z "${ADIOS2_PATH:-}" ]; then
    if [ -n "${ADIOS2_DIR:-}" ]; then
        ADIOS2_PATH="$ADIOS2_DIR"
    elif [ -n "${MAIN_DIR:-}" ] && [ -d "$MAIN_DIR/external/adios2" ]; then
        ADIOS2_PATH="$MAIN_DIR/external/adios2"
    fi
fi

if [ -n "${ADIOS2_PATH:-}" ]; then
    pyver=$(python3 -c 'import sys; print("{}.{}".format(sys.version_info.major, sys.version_info.minor))')
    export PYTHONPATH="${ADIOS2_PATH}/lib/python${pyver}/site-packages"
    export LD_LIBRARY_PATH="${ADIOS2_PATH}/lib:${ADIOS2_PATH}/lib64:${LD_LIBRARY_PATH:-}"
    echo "Using ADIOS2_PATH=${ADIOS2_PATH}"
else
    echo "Warning: ADIOS2_PATH is not set; in-situ Python may fail to import ADIOS2." >&2
fi


# Optional sanity prints
echo "Using Python:      $(which python3)"
echo "Using mpirun:      $(which mpirun)"
echo "CONDA_PREFIX:      ${CONDA_PREFIX}"
echo "PYTHONPATH (start):${PYTHONPATH:-<unset>}"
echo "LD_LIBRARY_PATH:   ${LD_LIBRARY_PATH:-<unset>}"
echo "---------------------------------------------------------"

CASE_FILE=${1:-cylinder_POD.case}

# Mesh / preprocessing
gmsh -0 newton.geo

gmsh2nek <<EOF
2
newton
0
EOF

rea2nbin newton.re2 ext_cyl.nmsh

# Coupled run: Neko + Python in-situ
mpirun -n 4 python3 insitu_task.py "$CASE_FILE" > python.log &
sleep 3
mpirun -n 4 ./neko "$CASE_FILE" > neko.log
