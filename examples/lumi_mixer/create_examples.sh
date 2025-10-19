#!/bin/bash

# Global parameters
if [ -z $SLURM_JOB_NAME ]; then
    export MAIN_DIR=$(realpath $(dirname $0)/../..)
else
    export MAIN_DIR=$(realpath $(dirname ./))
fi

if [ "$(basename $MAIN_DIR)" != "neko-top" ]; then
    echo "Invalid MAIN_DIR: $MAIN_DIR" >&2
    exit 1
fi

[ "$1" == "mesh" ] && MESH="true" || MESH="false"

export mesh_pattern="mixer"

export hpc_path="${MAIN_DIR}/scripts/jobscripts"
export data_path="${MAIN_DIR}/data_local/static_mixer"
export example_path="${MAIN_DIR}/examples/lumi_mixer"

export case_path="${example_path}/cases"
export template_path="${example_path}/templates"
export experiment_path="${example_path}/experiments"

[ ! -d "${data_path}" ] && mkdir -p "${data_path}"
[ ! -d "${case_path}" ] && mkdir -p "${case_path}"
[ ! -d "${experiment_path}" ] && mkdir -p "${experiment_path}"

function create_case() {

    if [ $# -lt 6 ]; then
        echo "ERROR: Not enough input arguments" >&2
        exit 1
    fi

    local experiment="$1"     # Experiment type (weak_scaling)
    local Nx="$2"             # Mesh dimension in x
    local Ny="$3"             # Mesh dimension in y
    local Nz="$4"             # Mesh dimension in z
    local cluster="$5"        # HPC cluster (LUMI-G)
    local nodes="$6"          # Number of nodes

    [ $# -ge 7 ] && local N_memory=$7 || local N_memory=""
    [ $# -ge 8 ] && local keep_files=$8 || local keep_files=false

    if [ "${cluster}" == "LUMI-G" ]; then
        Np="8"
    fi

    if [ -z "$N_memory" ]; then
        N_memory=$(grep -oP '(?<="n_memory": )[^,]*' "${template_path}/case.template")
    fi

    # Define a unique case name
    case_name="${cluster,,}_nodes_${nodes}_mesh_${Nx}x${Ny}x${Nz}"
    case_name="${case_name}_n_memory_${N_memory}"
    [ "${keep_files}" == "true" ] && case_name="${case_name}_keep_files"

    # Set file names
    local mesh_file="${data_path}/${mesh_pattern}_${Nx}x${Ny}x${Nz}.nmsh"
    local case_file="${case_path}/${case_name}.case"
    local job_path="${hpc_path}/${cluster}/lumi_mixer/cases"

    # Create the mesh if it does not exist
    if [[ ! -f "${mesh_file}" && "$MESH" == "true" ]]; then
        ./mesh.sh -b 0 4 0 1 0 1 $Nx $Ny $Nz \
           -o "${data_path}" -f ${mesh_pattern}_${Nx}x${Ny}x${Nz}.nmsh
    fi

    # Create the cases from the templates
    cp "${template_path}/case.template" "${case_file}"
    sed -i "s|\"mesh_file\": .*|\"mesh_file\": \"${mesh_file}\",|g" "${case_file}"
    sed -i "s|\"n_memory\": .*|\"n_memory\": ${N_memory},|g" "${case_file}"
    sed -i "s|\"keep_checkpoints\": .*|\"keep_checkpoints\": ${keep_files},|g" "${case_file}"

    [ ! -d "${job_path}" ] && mkdir -p ${job_path}
    [ ! -f ${job_path}/.gitignore ] && echo "*.sh" > ${job_path}/.gitignore

    # Create the jobscript
    cp "${template_path}/${cluster}.sh" "${job_path}/${case_name}.sh"
    sed -i "s|--nodes=.*|--nodes=${nodes}|g" "${job_path}/${case_name}.sh"

    # Create the experiment entry
    experiment_file="${experiment_path}/${experiment}.csv"
    if [ ! -f ${experiment_file} ]; then
        header="case_name, cluster, Nx, Ny, Nz, nodes, Np"
        [ -n "${N_memory}" ] && header+=", N_memory"

        echo "$header" >> ${experiment_file}
    fi

    if ! grep -q "^${case_name}\," "${experiment_file}"; then
        # Determine the data line
        data_line="${case_name}, ${cluster}, ${Nx}, ${Ny}, ${Nz}, ${nodes}, ${Np}"
        [ -n "${N_memory}" ] && data_line+=", ${N_memory}"
        printf "${data_line}\n" >> "${experiment_file}"
    fi
}

# Global settings
cluster="LUMI-G"

# Clean up old cases
find "${case_path}" -type f -name "*.case" -delete
find "${experiment_path}" -type f -name "*.csv" -delete
find "${hpc_path}/${cluster}/lumi_mixer/cases" -type f -name "*.sh" -delete

# Check max capacity of a single node
experiment="single_node_capacity"
create_case ${experiment}  16   4   4 "${cluster}" 1
create_case ${experiment}  16   4   4 "${cluster}" 1
create_case ${experiment}  32   8   8 "${cluster}" 1
create_case ${experiment}  64   8   8 "${cluster}" 1
create_case ${experiment}  64  16   8 "${cluster}" 1
create_case ${experiment}  64  16  16 "${cluster}" 1
create_case ${experiment} 128  16  16 "${cluster}" 1
create_case ${experiment} 128  32  16 "${cluster}" 1
create_case ${experiment} 128  32  32 "${cluster}" 1

# Investigation of the effect of memory checkpoints
experiment="memory_checkpoints"
create_case ${experiment} 128  16  16 "${cluster}" 1 50
create_case ${experiment} 128  16  16 "${cluster}" 1 100
create_case ${experiment} 128  16  16 "${cluster}" 1 200
create_case ${experiment} 128  16  16 "${cluster}" 1 400

create_case ${experiment} 128  32  16 "${cluster}" 1 50
create_case ${experiment} 128  32  16 "${cluster}" 1 100
create_case ${experiment} 128  32  16 "${cluster}" 1 200

create_case ${experiment} 128  32  32 "${cluster}" 1 50
create_case ${experiment} 128  32  32 "${cluster}" 1 100

# Investigate the scaling on a couple of nodes
experiment="weak_scaling_400"
create_case ${experiment} 128  16  16 "${cluster}"   1 400
create_case ${experiment}  64  32  32 "${cluster}"   2 400
create_case ${experiment} 128  32  32 "${cluster}"   4 400
create_case ${experiment} 256  32  32 "${cluster}"   8 400
create_case ${experiment} 128  64  64 "${cluster}"  16 400

experiment="weak_scaling_200"
create_case ${experiment}  64  32  32 "${cluster}"  1 200
create_case ${experiment} 128  32  32 "${cluster}"  2 200
create_case ${experiment} 256  32  32 "${cluster}"  4 200
create_case ${experiment} 128  64  64 "${cluster}"  8 200
create_case ${experiment} 256  64  64 "${cluster}" 16 200

experiment="weak_scaling_100"
create_case ${experiment} 128  32  32 "${cluster}"  1 100
create_case ${experiment} 256  32  32 "${cluster}"  2 100
create_case ${experiment} 128  64  64 "${cluster}"  4 100
create_case ${experiment} 256  64  64 "${cluster}"  8 100
create_case ${experiment} 512  64  64 "${cluster}" 16 100

experiment="weak_scaling"
create_case ${experiment}  128   32   32 "${cluster}"   1 100
create_case ${experiment}  256   32   32 "${cluster}"   2 100
create_case ${experiment}  128   64   64 "${cluster}"   4 100
create_case ${experiment}  256   64   64 "${cluster}"   8 100
# create_case ${experiment}  512   64   64 "${cluster}"  16 100 # Dies on IO

# experiment="strong_scaling_small"
# create_case ${experiment} 128  32  32 "${cluster}" 1
# create_case ${experiment} 128  32  32 "${cluster}" 2
# create_case ${experiment} 128  32  32 "${cluster}" 4
# create_case ${experiment} 128  32  32 "${cluster}" 8
# create_case ${experiment} 128  32  32 "${cluster}" 16
# create_case ${experiment} 128  32  32 "${cluster}" 32
# create_case ${experiment} 128  32  32 "${cluster}" 64

# experiment="strong_scaling_large"
# create_case ${experiment} 1024 256 256 "${cluster}" 16
# create_case ${experiment} 1024 256 256 "${cluster}" 32
# create_case ${experiment} 1024 256 256 "${cluster}" 64
# create_case ${experiment} 1024 256 256 "${cluster}" 128
# create_case ${experiment} 1024 256 256 "${cluster}" 256
# create_case ${experiment} 1024 256 256 "${cluster}" 512

experiment="storage_test"
create_case ${experiment}  16   4   4 "${cluster}" 1 100 true
create_case ${experiment}  16   8   4 "${cluster}" 1 100 true
create_case ${experiment}  16   8   8 "${cluster}" 1 100 true
create_case ${experiment}  32   8   8 "${cluster}" 1 100 true

create_case ${experiment}  16   4   4 "${cluster}" 1 100 true
create_case ${experiment}  16   8   4 "${cluster}" 2 100 true
create_case ${experiment}  16   8   8 "${cluster}" 4 100 true

experiment="run_time_test"
# Two node runs with different mesh sizes and memory settings
# Small
create_case ${experiment}  32   8   8 "${cluster}" 2 100
create_case ${experiment}  32   8   8 "${cluster}" 2 200
create_case ${experiment}  32   8   8 "${cluster}" 2 400
create_case ${experiment}  32   8   8 "${cluster}" 2 800

# Medium
create_case ${experiment}  64  16  16 "${cluster}" 2 100
create_case ${experiment}  64  16  16 "${cluster}" 2 200
create_case ${experiment}  64  16  16 "${cluster}" 2 400
create_case ${experiment}  64  16  16 "${cluster}" 2 800

# Large
create_case ${experiment} 128  32  32 "${cluster}" 2 100
create_case ${experiment} 128  32  32 "${cluster}" 2 200
create_case ${experiment} 128  32  32 "${cluster}" 2 400
create_case ${experiment} 128  32  32 "${cluster}" 2 800

# Four node runs with different mesh sizes and memory settings
# Small
create_case ${experiment}  32   8   8 "${cluster}" 4 100
create_case ${experiment}  32   8   8 "${cluster}" 4 200
create_case ${experiment}  32   8   8 "${cluster}" 4 400
create_case ${experiment}  32   8   8 "${cluster}" 4 800

# Medium
create_case ${experiment}  64  16  16 "${cluster}" 4 100
create_case ${experiment}  64  16  16 "${cluster}" 4 200
create_case ${experiment}  64  16  16 "${cluster}" 4 400
create_case ${experiment}  64  16  16 "${cluster}" 4 800

# Large
create_case ${experiment} 128  32  32 "${cluster}" 4 100
create_case ${experiment} 128  32  32 "${cluster}" 4 200
create_case ${experiment} 128  32  32 "${cluster}" 4 400
create_case ${experiment} 128  32  32 "${cluster}" 4 800

# Extra large
create_case ${experiment} 256  32  32 "${cluster}" 4 100
create_case ${experiment} 256  32  32 "${cluster}" 4 200
create_case ${experiment} 256  32  32 "${cluster}" 4 400
create_case ${experiment} 256  32  32 "${cluster}" 4 800

create_case ${experiment} 128  64  64 "${cluster}" 4 100
create_case ${experiment} 128  64  64 "${cluster}" 4 200
create_case ${experiment} 128  64  64 "${cluster}" 4 400
create_case ${experiment} 128  64  64 "${cluster}" 4 800

# Eight node runs with different mesh sizes and memory settings
# Medium
create_case ${experiment}  64  16  16 "${cluster}" 8 100
create_case ${experiment}  64  16  16 "${cluster}" 8 200
create_case ${experiment}  64  16  16 "${cluster}" 8 400
create_case ${experiment}  64  16  16 "${cluster}" 8 800

# Large
create_case ${experiment} 128  32  32 "${cluster}" 8 100
create_case ${experiment} 128  32  32 "${cluster}" 8 200
create_case ${experiment} 128  32  32 "${cluster}" 8 400
create_case ${experiment} 128  32  32 "${cluster}" 8 800

# Extra large
create_case ${experiment} 256  32  32 "${cluster}" 8 100
create_case ${experiment} 256  32  32 "${cluster}" 8 200
create_case ${experiment} 256  32  32 "${cluster}" 8 400
create_case ${experiment} 256  32  32 "${cluster}" 8 800

create_case ${experiment} 128  64  64 "${cluster}" 8 100
create_case ${experiment} 128  64  64 "${cluster}" 8 200
create_case ${experiment} 128  64  64 "${cluster}" 8 400
create_case ${experiment} 128  64  64 "${cluster}" 8 800

create_case ${experiment} 256  64  64 "${cluster}" 8 100
create_case ${experiment} 256  64  64 "${cluster}" 8 200
create_case ${experiment} 256  64  64 "${cluster}" 8 400
create_case ${experiment} 256  64  64 "${cluster}" 8 800
