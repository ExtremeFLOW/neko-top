#!/bin/bash

mkdir -p data_local/static_mixer

mesh_pattern=mixer
data_path=data_local/static_mixer
example_path=examples/static_mixer

find $example_path/petsc -name "*.case" -type f -delete

N=(8 16 32 64 128 256 512 1024 2048 4096)
Pe=(2000)
Re=(1 500)

for pe in ${Pe[@]}; do
    for re in ${Re[@]}; do
        for n in ${N[@]}; do
            Nx=$n
            Ny=$((n / 4))
            Nz=$((n / 4))
            output_file="${data_path}/${mesh_pattern}_${Nx}x${Ny}x${Nz}.nmsh"
            case_file="${example_path}/petsc/${n}_re_${re}_pe_${pe}.case"
            cache_file="${data_path}/petsc/cache_${n}_0.nek5000"

            if [[ "$1" == "mesh" && ! -f ${output_file} ]]; then
                ./mesh.sh -b 0 4 0 1 0 1 $Nx $Ny $Nz \
                    -o ${data_path} -f ${mesh_pattern}_${Nx}x${Ny}x${Nz}.nmsh
            fi

            # Create the cases from the templates
            cp ${example_path}/petsc/case.template ${case_file}
            sed -i "s|\"Re\": .*|\"Re\": ${re},|g" ${case_file}
            sed -i "s|\"Pe\": .*|\"Pe\": ${pe},|g" ${case_file}
            sed -i "s|\"mesh_file\": .*|\"mesh_file\": \"${output_file}\",|g" ${case_file}
            sed -i "s|\"cache_file\": .*|\"cache_file\": \"${data_path}/petsc/cache_${n}_\",|g" ${case_file}

            if [ "$re" == 1 ]; then
                sed -i "s|\"limits\": .*|\"limits\": [0.0, 5000.0],|g" ${case_file}
            else
                sed -i "s|\"limits\": .*|\"limits\": [0.0, 1000.0],|g" ${case_file}
            fi

            if [[ "$1" == "mesh" && ! -f "${cache_file}" ]]; then
                sed -i "s|\"end_time\": .*|\"end_time\": 1e-3,|g" ${case_file}
                ./run.sh -c static_mixer/petsc/${n}_re_${re}_pe_${pe}.case
                sed -i "s|\"end_time\": .*|\"end_time\": 30.0,|g" ${case_file}
            fi

            if [ ! -f "${cache_file}" ]; then
                rm -fr "${case_file}"
            fi
        done
    done
done
