#!/bin/bash

# Global parameters
export ROOT_DIR=$(realpath $(dirname $0))
export MAIN_DIR=$(realpath $ROOT_DIR/../..)

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
    local experiment=$1     # Experiment type (weak_scaling)
    local Nx=$2             # Mesh dimension in x
    local Ny=$3             # Mesh dimension in y
    local Nz=$4             # Mesh dimension in z
    local cluster=$5        # HPC cluster (LUMI-G)
    local nodes=$6          # Number of nodes

    if [ ${cluster} == "LUMI-G" ]; then
        Np=$((nodes * 8))
    fi

    # Set file names
    local case_name="${cluster,,}_nodes_${nodes}_mesh_${Nx}x${Ny}x${Nz}"
    local mesh_file="${data_path}/${mesh_pattern}_${Nx}x${Ny}x${Nz}.nmsh"
    local case_file="${case_path}/${case_name}.case"
    local job_path="${hpc_path}/${cluster}/lumi_mixer/cases"

    # Create the mesh if it does not exist
    if [ ! -f ${mesh_file} ]; then
        ./mesh.sh -b 0 4 0 1 0 1 $Nx $Ny $Nz -p $Np \
           -o ${data_path} -f ${mesh_pattern}_${Nx}x${Ny}x${Nz}.nmsh
    fi

    # Create the cases from the templates
    cp ${template_path}/case.template ${case_file}
    sed -i "s|\"mesh_file\": .*|\"mesh_file\": \"${mesh_file}\",|g" ${case_file}

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
create_case $experiment  64  16  16 ${cluster} 1
create_case $experiment 128  32  32 ${cluster} 1

# Too big even with just 5 saves to memory
# create_case $experiment 256  64  64 ${cluster} 1
