#!/usr/bin/bash

function help() {
    echo -e "run.sh [options]"
    echo -e "Run rugby_ball for sem_map_option = 0..6."
    echo -e ""
    echo -e "Options:"
    echo -e "  -h, --help      Show this help message and exit."
    echo -e "  -np=#           MPI ranks (default: 1)."
    echo -e "  --max-iters=#   Override optimization.solver.max_iterations."
    echo -e "  --clean-heavy   Remove transient field output files."
    echo -e "                  (default keeps heavy files,"
    echo -e "                   including design0.f*)"
    echo -e ""
    exit 0
}

cleanup_transients() {
    local run_dir=$1

    find "${run_dir}" -maxdepth 1 -type f -name "*.f*" -delete
    find "${run_dir}" -maxdepth 1 -type f -name "*.nek5000" -delete
    find "${run_dir}" -maxdepth 1 -type f -name "checkpoint*.chkp" -delete
}

CURRENT_DIR=$(pwd)
WORKING_DIR=$(dirname "$0")
WORKING_DIR=$(realpath "${WORKING_DIR}")
ROOT_DIR=$(realpath "${WORKING_DIR}/../../../")

NP=1
MAX_ITERS=""
KEEP_HEAVY=1
for arg in "$@"; do
    if [ "${arg:0:2}" == "--" ]; then
        case ${arg:2} in
        help) help ;;
        max-iters=*) MAX_ITERS=${arg:11} ;;
        keep-heavy) KEEP_HEAVY=1 ;;
        clean-heavy) KEEP_HEAVY=0 ;;
        *) echo -e "Invalid option: $arg" >&2 && help ;;
        esac
    elif [ "${arg:0:1}" == "-" ]; then
        case ${arg:1:1} in
        h) help ;;
        n)
           if [ "${arg:0:4}" == "-np=" ]; then
              NP=${arg:4}
           else
              echo -e "Invalid option: $arg" >&2 && help
           fi
           ;;
        *) echo -e "Invalid option: ${arg:1}" >&2 && help ;;
        esac
    fi
done

BIN="${WORKING_DIR}/reg_sem_options_bin"
CASE_SRC="${ROOT_DIR}/logs/rugby_ball/rugby_ball.case"
MESH_SRC="${ROOT_DIR}/logs/rugby_ball/box.nmsh"
RESULTS_DIR="${WORKING_DIR}/results"
SUMMARY_CSV="${RESULTS_DIR}/summary.csv"
OPTIONS=(0 1 2 3 4 5 6)
FAILURES=0
FAILED_OPTS=()

if [ ! -x "${BIN}" ]; then
    echo "Missing executable: ${BIN}" >&2
    echo "Build target 'reg_sem_options_bin' first." >&2
    exit 1
fi
if [ ! -f "${CASE_SRC}" ]; then
    echo "Missing source case: ${CASE_SRC}" >&2
    exit 1
fi
if [ ! -f "${MESH_SRC}" ]; then
    echo "Missing source mesh: ${MESH_SRC}" >&2
    exit 1
fi

mkdir -p "${RESULTS_DIR}"
echo "option,iteration,objective,constraint_1,run_directory,status,error" > \
    "${SUMMARY_CSV}"

for opt in "${OPTIONS[@]}"; do
    RUN_DIR="${RESULTS_DIR}/sem_map_option_${opt}"
    RUN_STATUS="ok"
    RUN_ERROR=""
    CSV_FIELDS=",,"

    rm -rf "${RUN_DIR}"
    mkdir -p "${RUN_DIR}"

    cp "${CASE_SRC}" "${RUN_DIR}/rugby_ball.case"
    cp "${MESH_SRC}" "${RUN_DIR}/box.nmsh"
    ln -sfn "${ROOT_DIR}/data" "${RUN_DIR}/data"
    ln -sfn "${ROOT_DIR}/data_local" "${RUN_DIR}/data_local"

    python3 - "${RUN_DIR}/rugby_ball.case" "${opt}" "${MAX_ITERS}" <<'PY'
import json
import sys

case_path = sys.argv[1]
sem_map_option = int(sys.argv[2])
max_iters = sys.argv[3].strip()

with open(case_path, "r", encoding="utf-8") as f:
    data = json.load(f)

data["optimization"]["design"]["sem_map_option"] = sem_map_option
if max_iters:
    data["optimization"]["solver"]["max_iterations"] = int(max_iters)

with open(case_path, "w", encoding="utf-8") as f:
    json.dump(data, f, indent=4)
    f.write("\n")
PY
    if [ "$?" -ne 0 ]; then
        RUN_STATUS="case_patch_failed"
        RUN_ERROR="failed_to_update_case_json"
    fi

    if [ "${RUN_STATUS}" = "ok" ]; then
        export NEKO_LOG_FILE="neko_sem_map_option_${opt}.log"
        echo "Running sem_map_option=${opt} with mpirun -n ${NP}"
        (
            cd "${RUN_DIR}" || exit 1
            mpirun -n "${NP}" "${BIN}" "rugby_ball.case" > run.log 2>&1
        )
        if [ "$?" -ne 0 ]; then
            RUN_STATUS="run_failed"
            RUN_ERROR="mpirun_failed"
        fi
    fi

    if [ "${RUN_STATUS}" = "ok" ]; then
        if [ ! -f "${RUN_DIR}/optimization_data.csv" ]; then
            RUN_STATUS="missing_optimization_data"
            RUN_ERROR="optimization_data_missing"
        fi
    fi

    if [ "${RUN_STATUS}" = "ok" ]; then
        CSV_FIELDS=$(python3 - "${RUN_DIR}/optimization_data.csv" <<'PY'
import csv
import sys

path = sys.argv[1]
with open(path, "r", encoding="utf-8", newline="") as f:
    rows = list(csv.DictReader(f))

if not rows:
    print(",,")
else:
    last = rows[-1]
    iteration = last.get("iteration", "")
    objective = last.get("objective", "")
    constraint = last.get("constraint_1", "")
    if not constraint:
        for key in last.keys():
            if key.startswith("constraint"):
                constraint = last[key]
                break
    print(f"{iteration},{objective},{constraint}")
PY
)
        if [ "$?" -ne 0 ]; then
            RUN_STATUS="summary_parse_failed"
            RUN_ERROR="failed_to_parse_optimization_data"
            CSV_FIELDS=",,"
        fi
    fi

    if [ "${RUN_STATUS}" != "ok" ]; then
        FAILURES=$((FAILURES + 1))
        FAILED_OPTS+=("${opt}:${RUN_STATUS}")
        echo "Warning: sem_map_option=${opt} failed (${RUN_STATUS})" >&2
    fi

    echo "${opt},${CSV_FIELDS},${RUN_DIR},${RUN_STATUS},${RUN_ERROR}" >> \
        "${SUMMARY_CSV}"

    if [ "${KEEP_HEAVY}" -eq 0 ]; then
        cleanup_transients "${RUN_DIR}"
    fi
done

echo "Completed SEM option regression runs."
echo "Summary: ${SUMMARY_CSV}"
if [ "${FAILURES}" -gt 0 ]; then
    echo "Failures (${FAILURES}): ${FAILED_OPTS[*]}" >&2
fi

cd "${CURRENT_DIR}" || exit 1
