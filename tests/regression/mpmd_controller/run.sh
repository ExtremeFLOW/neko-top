#!/bin/bash
set -euo pipefail

case_path=$1
python_peer=$2
neko_driver=$3
log_file=$4

repo_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../../.." && pwd)

# shellcheck disable=SC1091
source "${repo_root}/scripts/mpmd_run_helpers.sh"

mpmd_validate_mpi4py
if ! mpmd_launch_shared "${case_path}" "${python_peer}" "${neko_driver}" \
    1 1 "${log_file}"
then
    cat "${log_file}" >&2 || true
    exit 1
fi

if ! grep -q "MPMD controller regression passed" "${log_file}"; then
    echo "Error: Python MPMD controller peer did not report success." >&2
    cat "${log_file}" >&2 || true
    exit 1
fi
