#!/bin/bash

mkdir -p data_local/static_mixer

mesh_pattern=mixer
data_path=data_local/static_mixer
example_path=examples/static_mixer

N=(8 16 32 64 128 256 512 1024 2048 4096)

for n in ${N[@]}; do
    Nx=$n
    Ny=$((n / 4))
    Nz=$((n / 4))
    output_file="${data_path}/${mesh_pattern}_${Nx}x${Ny}x${Nz}.nmsh"

    if [ ! -f ${output_file} ]; then
        ./mesh.sh -b 0 4 0 1 0 1 $Nx $Ny $Nz
         \
            -o ${data_path} -f ${mesh_pattern}_${Nx}x${Ny}x${Nz}.nmsh
    fi

    # Create the cases from the templates
    cp ${example_path}/petsc/case.template ${example_path}/petsc/${n}.case
    sed -i "s|\"mesh_file\": .*|\"mesh_file\": \"${output_file}\",|g" ${example_path}/petsc/${n}.case
    sed -i "s|\"cache_file\": .*|\"cache_file\": \"${data_path}/petsc/cache_${n}_\",|g" ${example_path}/petsc/${n}.case

    cp ${example_path}/neko_top/case.template ${example_path}/neko_top/${n}.case
    sed -i "s|\"mesh_file\": .*|\"mesh_file\": \"${output_file}\",|g" ${example_path}/neko_top/${n}.case
done
