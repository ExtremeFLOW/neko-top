#!/bin/bash

# Global parameters
export ROOT_DIR=$(realpath $(dirname $0))
export MAIN_DIR=$(realpath $ROOT_DIR/../..)

export mesh_pattern="mixer"
export data_path="${MAIN_DIR}/data_local/static_mixer"
export example_path="${MAIN_DIR}/examples/lumi_mixer/cases"
export template_path="${MAIN_DIR}/examples/lumi_mixer/templates"
export hpc_path="${MAIN_DIR}/scripts/jobscripts"

[ ! -d ${data_path} ] && mkdir -p ${data_path}
[ ! -d ${example_path} ] && mkdir -p ${example_path}

function create_case() {
    local n=$1              # Number of cells in x direction
    local cluster=$2        # Name of the cluster
    local nodes=$3          # Number of nodes

    # Compute mesh dimensions
    local Nx=$((n / 4 * 4))
    local Ny=$((n / 4))
    local Nz=$((n / 4))

    if [ ${cluster} == "LUMI-G" ]; then
        Np=$((nodes * 8))
    fi

    if [[ Nx -ne n || Ny -lt 1 || Nz -lt 1 ]]; then
        echo "Error: Invalid mesh dimensions: ${n}."
        echo "       Must be a multiple of 4 and at least 4."
        exit 1
    fi

    # Set file names
    local case_name="nodes_${nodes}_mesh_${n}"
    local mesh_file="${data_path}/${mesh_pattern}_${Nx}x${Ny}x${Nz}_${Np}.nmsh"
    local case_file="${example_path}/${case_name}.case"
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

}

N=(8)
cluster="LUMI-G"
nodes=(2 1)
for n in ${N[@]}; do
    for nodes_ in ${nodes[@]}; do
        create_case ${n} ${cluster} ${nodes_}
    done
done
