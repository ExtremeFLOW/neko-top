#!/bin/bash
set -euo pipefail

script_dir="$(cd "$(dirname "$0")" && pwd)"

# default_modes=(4 8 16)
# default_streams=(1 5 50)
default_modes=(4)
default_streams=(5)

copy_support_files() {
    local out_dir="$1"

    cp "${script_dir}/prepare.sh" "${out_dir}/prepare.sh"
    cp "${script_dir}/run.sh" "${out_dir}/run.sh"
    chmod 755 "${out_dir}/prepare.sh" "${out_dir}/run.sh"
}

generate_checkpoint_case() {
    local out_dir

    out_dir="${script_dir}/checkpoint"
    mkdir -p "${out_dir}"

    cp "${script_dir}/checkpoint.case" "${out_dir}/case.case"
    copy_support_files "${out_dir}"
}

generate_pod_case() {
    local n_modes="$1"
    local i_stream="$2"
    local out_dir

    out_dir="${script_dir}/POD_nm${n_modes}_is${i_stream}"
    mkdir -p "${out_dir}"

    cp "${script_dir}/pod.case" "${out_dir}/case.case"
    sed -i -E \
        "s/(\"n_modes\":[[:space:]]*)[0-9]+/\\1${n_modes}/" \
        "${out_dir}/case.case"
    sed -i -E \
        "s/(\"i_stream\":[[:space:]]*)[0-9]+/\\1${i_stream}/" \
        "${out_dir}/case.case"

    copy_support_files "${out_dir}"
}

main() {
    local n_modes
    local i_stream

    generate_checkpoint_case

    for n_modes in "${default_modes[@]}"; do
        for i_stream in "${default_streams[@]}"; do
            generate_pod_case "${n_modes}" "${i_stream}"
        done
    done
}

main "$@"
