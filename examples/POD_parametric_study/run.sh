#!/bin/bash
set -euo pipefail

script_dir="$(cd "$(dirname "$0")" && pwd)"
root_dir="${script_dir}"
while [ ! -d "${root_dir}/sources" ] && [ "${root_dir}" != "/" ]; do
    root_dir="$(dirname "${root_dir}")"
done
if [ ! -d "${root_dir}/sources" ]; then
    echo "Error: could not locate repo root with sources/ from " \
         "${script_dir}" >&2
    exit 1
fi

cd "${script_dir}"

PY_SCRIPT=${2:-${root_dir}/sources/state_recovery/POD_state_recover/\
pod_state_recover.py}

is_parametric_root=false
if [ -x "${script_dir}/generate_cases.sh" ] && \
   [ -f "${script_dir}/checkpoint.case" ] && \
   [ -f "${script_dir}/pod.case" ]; then
    is_parametric_root=true
fi

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
        echo "Error: no Neko executable found " \
             "(set NEKO_BIN or provide ./neko)." >&2
        exit 1
    fi
fi

if [ -x "${script_dir}/neko" ] && [ -z "${NEKO_BIN:-}" ] && \
   [ "${script_dir}/neko" -ot "${root_dir}/sources/problem/problem.f90" ]; then
    echo "Warning: ./neko is older than sources/problem/problem.f90;" \
         " rebuild to pick up fixes." >&2
fi

NEKO_RANKS=${NEKO_RANKS:-10}
PY_RANKS=${PY_RANKS:-4}

if { [ "$is_parametric_root" = false ] || [ $# -ge 1 ]; } && \
   [ -x "${script_dir}/prepare.sh" ]; then
    (cd "${script_dir}" && ./prepare.sh)
fi

# conda
if [ -f "/scratch/nobis/miniconda3/etc/profile.d/conda.sh" ]; then
    . "/scratch/nobis/miniconda3/etc/profile.d/conda.sh"
else
    eval "$(/scratch/nobis/miniconda3/bin/conda shell.bash hook)"
fi

conda activate \
    /scratch/nobis/POSTDOC/PYSEMTOOLS/miniconda3/.conda/envs/pySEMTools

# Resolve ADIOS2_PATH
if [ -z "${ADIOS2_PATH:-}" ]; then
    if [ -n "${ADIOS2_DIR:-}" ]; then
        ADIOS2_PATH="$ADIOS2_DIR"
    elif [ -d "${root_dir}/external/adios2" ]; then
        ADIOS2_PATH="${root_dir}/external/adios2"
    elif [ -n "${MAIN_DIR:-}" ] && [ -d "$MAIN_DIR/external/adios2" ]; then
        ADIOS2_PATH="$MAIN_DIR/external/adios2"
    fi
fi

if [ -n "${ADIOS2_PATH:-}" ]; then
    pyver=$(python3 -c \
        'import sys; print("{}.{}".format(sys.version_info.major, \
sys.version_info.minor))')
    export PYTHONPATH="${ADIOS2_PATH}/lib/python${pyver}/site-packages:\
${PYTHONPATH:-}"
    export LD_LIBRARY_PATH="${ADIOS2_PATH}/lib:${ADIOS2_PATH}/lib64:\
${LD_LIBRARY_PATH:-}"
    echo "Using ADIOS2_PATH=${ADIOS2_PATH}"
else
    echo "Warning: ADIOS2_PATH is not set; in-situ Python may fail to " \
         "import ADIOS2." >&2
fi

echo "Using Python:      $(which python3)"
echo "Using mpirun:      $(which mpirun)"
echo "CONDA_PREFIX:      ${CONDA_PREFIX}"
echo "PYTHONPATH (start):${PYTHONPATH:-<unset>}"
echo "LD_LIBRARY_PATH:   ${LD_LIBRARY_PATH:-<unset>}"
echo "---------------------------------------------------------"

if [ "$is_parametric_root" = true ] && [ $# -lt 1 ]; then
    study_dirs=()

    ./generate_cases.sh

    if [ -d "${script_dir}/checkpoint" ]; then
        study_dirs+=("${script_dir}/checkpoint")
    fi

    while IFS= read -r study_dir; do
        study_dirs+=("${study_dir}")
    done < <(
        find "${script_dir}" -mindepth 1 -maxdepth 1 -type d \
            -name 'POD_nm*_is*' | sort
    )

    if [ ${#study_dirs[@]} -eq 0 ]; then
        echo "Error: no generated study cases found in ${script_dir}" >&2
        exit 1
    fi

    for study_dir in "${study_dirs[@]}"; do
        echo "========================================================="
        echo "Running study case: $(basename "${study_dir}")"
        echo "========================================================="
        "${study_dir}/run.sh" "${study_dir}/case.case" "${PY_SCRIPT}"
    done
    exit 0
fi

if [ $# -ge 1 ]; then
    case_path="$1"
else
    case_files=( *.case )
    if [ ${#case_files[@]} -eq 0 ]; then
        echo "Error: no .case files found in ${script_dir}" >&2
        exit 1
    elif [ ${#case_files[@]} -gt 1 ]; then
        echo "Error: multiple .case files found in ${script_dir}; pass one " \
             "explicitly." >&2
        exit 1
    fi
    case_path="${case_files[0]}"
fi

case_name="$(basename "${case_path}")"
base="${case_name%.case}"

if grep -q '"type"[[:space:]]*:[[:space:]]*"pod"' "${case_path}"; then
    log_file=${LOG_FILE:-"mpmd_${base}.log"}
    mpirun_cmd=(
        mpirun
        --tag-output
        -n "${PY_RANKS}"
        python3 "${PY_SCRIPT}" "${case_path}"
        :
        -n "${NEKO_RANKS}"
        "${NEKO_EXE}" "${case_path}"
    )

    echo "Launching shared MPI job:"
    printf '  %q' "${mpirun_cmd[@]}"
    printf '\n'
    echo "Combined output:   ${log_file}"

    "${mpirun_cmd[@]}" > "${log_file}" 2>&1
else
    log_file=${LOG_FILE:-"neko_${base}.log"}
    mpirun_cmd=(
        mpirun
        --tag-output
        -n "${NEKO_RANKS}"
        "${NEKO_EXE}" "${case_path}"
    )

    echo "Launching MPI job:"
    printf '  %q' "${mpirun_cmd[@]}"
    printf '\n'
    echo "Output:            ${log_file}"

    "${mpirun_cmd[@]}" > "${log_file}" 2>&1
fi
