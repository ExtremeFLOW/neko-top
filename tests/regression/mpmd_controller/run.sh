#!/bin/bash
set -euo pipefail

case_path=$1
python_peer=$2
neko_driver=$3
log_file=$4

if [ -z "${NEKO_TOP_RUN_MPMD_REGRESSION:-}" ]; then
    echo "Skipping MPMD controller regression: set" \
        "NEKO_TOP_RUN_MPMD_REGRESSION=1 to opt in."
    exit 77
fi

repo_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../../.." && pwd)

# shellcheck disable=SC1091
source "${repo_root}/scripts/mpmd_run_helpers.sh"

mpmd_ensure_runtime "${repo_root}"
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
