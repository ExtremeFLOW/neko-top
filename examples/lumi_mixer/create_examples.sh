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

export mesh_pattern="mixer"

export hpc_path="${MAIN_DIR}/scripts/jobscripts"
export data_path="${MAIN_DIR}/data_local/static_mixer"
export example_path="${MAIN_DIR}/examples/lumi_mixer"

export case_path="${example_path}/cases"
export template_path="${example_path}/templates"
export experiment_path="${example_path}/experiments"

[ ! -d ${data_path} ] && mkdir -p ${data_path}
[ ! -d ${case_path} ] && mkdir -p ${case_path}
[ ! -d ${experiment_path} ] && mkdir -p ${experiment_path}

function create_case() {

    if [ $# -lt 6 ]; then
        echo "ERROR: Not enough input arguments" >&2
        exit 1
    fi

    local experiment=$1     # Experiment type (weak_scaling)
    local Nx=$2             # Mesh dimension in x
    local Ny=$3             # Mesh dimension in y
    local Nz=$4             # Mesh dimension in z
    local cluster=$5        # HPC cluster (LUMI-G)
    local nodes=$6          # Number of nodes

    [ $# -ge 7 ] && local N_memory=$7 || local N_memory=""

    if [ ${cluster} == "LUMI-G" ]; then
        Np=$((nodes * 8))
    fi

    # Define a unique case name
    local case_name="${cluster,,}_nodes_${nodes}_mesh_${Nx}x${Ny}x${Nz}"
    [ -n "${N_memory}" ] && case_name="${case_name}_n_memory_${N_memory}"

    # Set file names
    local mesh_file="${data_path}/${mesh_pattern}_${Nx}x${Ny}x${Nz}.nmsh"
    local case_file="${case_path}/${case_name}.case"
    local job_path="${hpc_path}/${cluster}/lumi_mixer/cases"

    # Create the mesh if it does not exist
    if [ ! -f ${mesh_file} ]; then
        ./mesh.sh -b 0 4 0 1 0 1 $Nx $Ny $Nz \
           -o ${data_path} -f ${mesh_pattern}_${Nx}x${Ny}x${Nz}.nmsh
    fi

    # Create the cases from the templates
    cp ${template_path}/case.template ${case_file}
    sed -i "s|\"mesh_file\": .*|\"mesh_file\": \"${mesh_file}\",|g" ${case_file}

    if [ -n "${N_memory}" ]; then
        sed -i "s|\"n_memory\": .*|\"n_memory\": ${N_memory},|g" ${case_file}
    fi

    [ ! -d ${job_path} ] && mkdir -p ${job_path}
    [ ! -f ${job_path}/.gitignore ] && echo "*.sh" > ${job_path}/.gitignore

    # Create the jobscript
    cp ${template_path}/${cluster}.sh ${job_path}/${case_name}.sh
    sed -i "s|--nodes=.*|--nodes=${nodes}|g" ${job_path}/${case_name}.sh

    # Create the experiment entry
    if [ ! -f ${experiment_path}/${experiment}.txt ]; then
        echo "# Experiment: ${experiment}" > ${experiment_path}/${experiment}.txt
        echo "# Format: case_name cluster Nx Ny Nz Np" >> ${experiment_path}/${experiment}.txt
    fi

    if ! grep -q "^${case_name} " ${experiment_path}/${experiment}.txt; then
        echo "${case_name} ${cluster} ${Nx} ${Ny} ${Nz} ${Np}" >> ${experiment_path}/${experiment}.txt
    fi
}

find ${case_path} -type f -name "*.case" -delete

experiment="preliminary"
cluster="LUMI-G"
create_case $experiment   8   2   2 ${cluster} 1
create_case $experiment  16   4   4 ${cluster} 1
create_case $experiment  32   8   8 ${cluster} 1
create_case $experiment  64   8   8 ${cluster} 1
create_case $experiment  32  16   8 ${cluster} 1
create_case $experiment  32   8  16 ${cluster} 1
create_case $experiment  32   8   8 ${cluster} 1 10
create_case $experiment  32   8   8 ${cluster} 1 100
create_case $experiment  32   8   8 ${cluster} 1 1000
create_case $experiment  32   8   8 ${cluster} 1 10000
create_case $experiment  32   8   8 ${cluster} 1 100000

# experiment="weak_scaling"
# cluster="LUMI-G"
# create_case $experiment 128  32  32 ${cluster} 1
# create_case $experiment 256  32  32 ${cluster} 2
# create_case $experiment 256  64  32 ${cluster} 4
# create_case $experiment 256  64  64 ${cluster} 8
# create_case $experiment 512  64  64 ${cluster} 16
# create_case $experiment 512 128  64 ${cluster} 32
# create_case $experiment 512 128 128 ${cluster} 64

# experiment="strong_scaling_small"
# cluster="LUMI-G"
# create_case $experiment 128  32  32 ${cluster} 1
# create_case $experiment 128  32  32 ${cluster} 2
# create_case $experiment 128  32  32 ${cluster} 4
# create_case $experiment 128  32  32 ${cluster} 8
# create_case $experiment 128  32  32 ${cluster} 16
# create_case $experiment 128  32  32 ${cluster} 32
# create_case $experiment 128  32  32 ${cluster} 64

# experiment="strong_scaling_large"
# cluster="LUMI-G"
# create_case $experiment 1024 256 256 ${cluster} 16
# create_case $experiment 1024 256 256 ${cluster} 32
# create_case $experiment 1024 256 256 ${cluster} 64
# create_case $experiment 1024 256 256 ${cluster} 128
# create_case $experiment 1024 256 256 ${cluster} 256
