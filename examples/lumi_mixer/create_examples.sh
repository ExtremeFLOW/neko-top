#!/bin/bash

# Global parameters
export MAIN_DIR=$(realpath $(dirname $0)/../..)
export mesh_pattern="mixer"
export data_path=$(realpath "${MAIN_DIR}/data_local/static_mixer")
export example_path=$(realpath "${MAIN_DIR}/examples/static_mixer/cases")
export template_path=$(realpath "${MAIN_DIR}/examples/static_mixer/templates")

[ ! -d ${data_path} ] && mkdir -p ${data_path}
[ ! -d ${example_path} ] && mkdir -p ${example_path}

function create_case() {
    local case_template=$1  # Template file
    local n=$2              # Number of cells in x direction

    # Compute mesh dimensions
    local Nx=$((n / 4 * 4))
    local Ny=$((n / 4))
    local Nz=$((n / 4))

    if [[ Nx -ne n || Ny -lt 1 || Nz -lt 1 ]]; then
        echo "Error: Invalid mesh dimensions: ${n}."
        echo "       Must be a multiple of 4 and at least 4."
        exit 1
    fi

    # Set file names
    local mesh_file="${data_path}/${mesh_pattern}_${Nx}x${Ny}x${Nz}.nmsh"
    local case_path="${example_path}/${case_template}"
    local case_file="${case_path}/${n}.case"

    # Create the mesh if it does not exist
    if [ ! -f ${mesh_file} ]; then
        ./mesh.sh -b 0 4 0 1 0 1 $Nx $Ny $Nz \
            -o ${data_path} -f ${mesh_pattern}_${Nx}x${Ny}x${Nz}.nmsh
    fi

    # Copy over the CMakeFile
    [ ! -d ${case_path} ] && mkdir -p ${case_path}
    if [ ! -f ${case_path}/CMakeLists.txt ]; then
        cp ${template_path}/CMakeLists.${case_template} ${case_path}/CMakeLists.txt
    fi

    # Copy over the mixer.f90 file
    if [ ! -f ${case_path}/mixer.f90 ]; then
        cp ${template_path}/mixer.f90 ${case_path}/mixer.f90
    fi


    # Create the cases from the templates
    cp ${template_path}/case.${case_template} ${case_file}
    sed -i "s|\"mesh_file\": .*|\"mesh_file\": \"${mesh_file}\",|g" ${case_file}

    if [ ${case_template} == "petsc" ]; then
        sed -i "s|\"cache_file\": .*|\"cache_file\": \"${data_path}/petsc/cache_${n}_\",|g" ${case_file}
    fi
}

N=(4)

for n in ${N[@]}; do
    create_case petsc ${n}
    create_case neko_top ${n}
done
