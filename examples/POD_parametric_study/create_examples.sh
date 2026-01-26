#!/bin/bash
set -euo pipefail

example_path="$(cd "$(dirname "$0")" && pwd)"
template_file="${example_path}/POD.template"
checkpoint_template="${example_path}/checkpoint.template"
run_template="${example_path}/run_template.sh"
prepare_script="${example_path}/prepare.sh"

if [ ! -f "${template_file}" ]; then
    echo "Missing template: ${template_file}" >&2
    exit 1
fi
if [ ! -f "${checkpoint_template}" ]; then
    echo "Missing template: ${checkpoint_template}" >&2
    exit 1
fi
if [ ! -f "${run_template}" ]; then
    echo "Missing run template: ${run_template}" >&2
    exit 1
fi
# Sweep parameters (edit as needed).
N_MODES=(3 5 7)
I_STREAM=(5 10 50)

find "${example_path}" -maxdepth 1 -type f -name "*.case" -delete

for n_modes in "${N_MODES[@]}"; do
    for i_stream in "${I_STREAM[@]}"; do
        case_dir="${example_path}/POD_nm${n_modes}_is${i_stream}"
        case_file="${case_dir}/case.case"
        mkdir -p "${case_dir}"
        cp "${template_file}" "${case_file}"
        sed -i "s|__N_MODES__|${n_modes}|g" "${case_file}"
        sed -i "s|__I_STREAM__|${i_stream}|g" "${case_file}"
        cp "${run_template}" "${case_dir}/run.sh"
        chmod +x "${case_dir}/run.sh"
        if [ -f "${prepare_script}" ]; then
            cp "${prepare_script}" "${case_dir}/prepare.sh"
            chmod +x "${case_dir}/prepare.sh"
        fi
    done
done

checkpoint_dir="${example_path}/checkpoint"
checkpoint_case="${checkpoint_dir}/case.case"
mkdir -p "${checkpoint_dir}"
cp "${checkpoint_template}" "${checkpoint_case}"
cp "${run_template}" "${checkpoint_dir}/run.sh"
chmod +x "${checkpoint_dir}/run.sh"
if [ -f "${prepare_script}" ]; then
    cp "${prepare_script}" "${checkpoint_dir}/prepare.sh"
    chmod +x "${checkpoint_dir}/prepare.sh"
fi

case_count=$(find "${example_path}" -type f -name "case.case" | wc -l)
echo "Generated ${case_count} case files."
