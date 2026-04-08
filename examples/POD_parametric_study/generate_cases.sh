#!/bin/bash
set -euo pipefail

script_dir="$(cd "$(dirname "$0")" && pwd)"
default_dir="${script_dir}/default_case"
cases_dir="${script_dir}/cases"

default_modes=(4 8 16)
default_streams=(1 5 50)
default_memories=(100)
# default_modes=(4)
# default_streams=(5)

copy_support_files() {
    local out_dir="$1"

    cp "${default_dir}/prepare.sh" "${out_dir}/prepare.sh"
    cp "${default_dir}/run.sh.template" "${out_dir}/run.sh"
    chmod 755 "${out_dir}/prepare.sh" "${out_dir}/run.sh"
}

set_state_recovery_type() {
    local casefile="$1"
    local recover_type="$2"

    sed -i \
        "s/\"type\": \"pod\"/\"type\": \"${recover_type}\"/" \
        "${casefile}"
}

generate_checkpoint_case() {
    local n_memory="$1"
    local out_dir

    out_dir="${cases_dir}/checkpoint_nmem${n_memory}"
    mkdir -p "${out_dir}"

    cp "${default_dir}/case.template" "${out_dir}/case.case"
    set_state_recovery_type "${out_dir}/case.case" "checkpoint"
    sed -i -E \
        "s/(\"n_memory\":[[:space:]]*)[0-9]+/\\1${n_memory}/" \
        "${out_dir}/case.case"
    copy_support_files "${out_dir}"
}

generate_pod_case() {
    local n_modes="$1"
    local i_stream="$2"
    local out_dir

    out_dir="${cases_dir}/POD_nm${n_modes}_is${i_stream}"
    mkdir -p "${out_dir}"

    cp "${default_dir}/case.template" "${out_dir}/case.case"
    set_state_recovery_type "${out_dir}/case.case" "pod"
    sed -i -E \
        "s/(\"n_modes\":[[:space:]]*)[0-9]+/\\1${n_modes}/" \
        "${out_dir}/case.case"
    sed -i -E \
        "s/(\"i_stream\":[[:space:]]*)[0-9]+/\\1${i_stream}/" \
        "${out_dir}/case.case"

    copy_support_files "${out_dir}"
}

main() {
    local n_memory
    local n_modes
    local i_stream

    rm -rf "${cases_dir}"

    for n_memory in "${default_memories[@]}"; do
        generate_checkpoint_case "${n_memory}"
    done

    for n_modes in "${default_modes[@]}"; do
        for i_stream in "${default_streams[@]}"; do
            generate_pod_case "${n_modes}" "${i_stream}"
        done
    done
}

main "$@"
