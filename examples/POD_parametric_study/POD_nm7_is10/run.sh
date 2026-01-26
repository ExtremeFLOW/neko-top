#!/bin/bash
set -euo pipefail

script_dir="$(cd "$(dirname "$0")" && pwd)"
root_dir="${script_dir}"
while [ ! -d "${root_dir}/sources" ] && [ "${root_dir}" != "/" ]; do
    root_dir="$(dirname "${root_dir}")"
done
if [ ! -d "${root_dir}/sources" ]; then
    echo "Error: could not locate repo root with sources/ from ${script_dir}" >&2
    exit 1
fi

PY_SCRIPT=${2:-${root_dir}/sources/state_recovery/POD_state_recover/pod_state_recover.py}

if [ -n "${NEKO_BIN:-}" ]; then
    NEKO_EXE="${NEKO_BIN}"
elif [ -x "${script_dir}/neko" ]; then
    NEKO_EXE="${script_dir}/neko"
else
    example_dir="${script_dir}"
    if [[ "${script_dir}" == "${root_dir}/logs/"* ]]; then
        example_dir="${script_dir/${root_dir}\/logs/${root_dir}\/examples}"
    fi
    example_parent="$(dirname "${example_dir}")"
    if [ -x "${example_parent}/neko" ]; then
        NEKO_EXE="${example_parent}/neko"
    else
        echo "Error: no Neko executable found (set NEKO_BIN or provide ./neko)." >&2
        exit 1
    fi
fi

if [ -x "${script_dir}/neko" ] && [ -z "${NEKO_BIN:-}" ] && \
   [ "${script_dir}/neko" -ot "${root_dir}/sources/problem/problem.f90" ]; then
    echo "Warning: ./neko is older than sources/problem/problem.f90; rebuild to pick up fixes." >&2
fi

NEKO_RANKS=${NEKO_RANKS:-10}
PY_RANKS=${PY_RANKS:-4}
PY_STARTUP_DELAY=${PY_STARTUP_DELAY:-3}

if [ -x "${script_dir}/prepare.sh" ]; then
    (cd "${script_dir}" && ./prepare.sh)
fi

# conda
if [ -f "/scratch/nobis/miniconda3/etc/profile.d/conda.sh" ]; then
    . "/scratch/nobis/miniconda3/etc/profile.d/conda.sh"
else
    eval "$(/scratch/nobis/miniconda3/bin/conda shell.bash hook)"
fi

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
    export PYTHONPATH="${ADIOS2_PATH}/lib/python${pyver}/site-packages:${PYTHONPATH:-}"
    export LD_LIBRARY_PATH="${ADIOS2_PATH}/lib:${ADIOS2_PATH}/lib64:${LD_LIBRARY_PATH:-}"
    echo "Using ADIOS2_PATH=${ADIOS2_PATH}"
else
    echo "Warning: ADIOS2_PATH is not set; in-situ Python may fail to import ADIOS2." >&2
fi

echo "Using Python:      $(which python3)"
echo "Using mpirun:      $(which mpirun)"
echo "CONDA_PREFIX:      ${CONDA_PREFIX}"
echo "PYTHONPATH (start):${PYTHONPATH:-<unset>}"
echo "LD_LIBRARY_PATH:   ${LD_LIBRARY_PATH:-<unset>}"
echo "---------------------------------------------------------"

if [ $# -ge 1 ]; then
    case_path="$1"
else
    case_files=( "${script_dir}"/*.case )
    if [ ${#case_files[@]} -eq 0 ]; then
        echo "Error: no .case files found in ${script_dir}" >&2
        exit 1
    elif [ ${#case_files[@]} -gt 1 ]; then
        echo "Error: multiple .case files found in ${script_dir}; pass one explicitly." >&2
        exit 1
    fi
    case_path="${case_files[0]}"
fi

case_name="$(basename "${case_path}")"
base="${case_name%.case}"
python_log="${script_dir}/python_${base}.log"
neko_log="${script_dir}/neko_${base}.log"

if grep -q '"type"[[:space:]]*:[[:space:]]*"pod"' "${case_path}"; then
    mpirun --tag-output -n "${PY_RANKS}" python3 "${PY_SCRIPT}" "${case_path}" > "${python_log}" 2>&1 &
    PY_PID=$!
    sleep "${PY_STARTUP_DELAY}"
fi

mpirun --tag-output -n "${NEKO_RANKS}" "${NEKO_EXE}" "${case_path}" > "${neko_log}" 2>&1
