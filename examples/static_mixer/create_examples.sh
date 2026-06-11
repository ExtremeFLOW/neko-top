#!/bin/bash

example_list=(petsc petsc_phmg)
N=(64 128)
Re=(1000 1500 3000)

example_path=examples/static_mixer
data_path=data_local/static_mixer
mesh_pattern=mixer
job_path=scripts/jobscripts/LUMI-G/static_mixer

mkdir -p data_local/static_mixer

for example in "${example_list[@]}"; do

    mkdir -p ${job_path}/$example
    echo ".gitignore" > ${job_path}/$example/.gitignore

    find ${example_path}/${example} -name "*.case" -type f -delete


    # LUMI specific GPU settings
    N_per_gpu=8192
    N_gpu=8

    for n in ${N[@]}; do
        for re in ${Re[@]}; do
            pe=$re
            Nx=$n
            Ny=$((n / 4))
            Nz=$((n / 4))
            nodes=$((Nx * Ny * Nz / N_per_gpu / N_gpu))
            [ $nodes -lt 1 ] && nodes=1 # Ensure at least 1 node is requested

            # Calculate the pde filter radius required: Nx / 4.0
            element_size=$(printf "%.2e" $(echo "scale=2; 4.0 / $Nx" | bc -l))

            output_file="${data_path}/${mesh_pattern}_${Nx}x${Ny}x${Nz}.nmsh"
            case_file="${example_path}/${example}/${n}_re_${re}_pe_${pe}.case"
            cache_file="${data_path}/petsc/cache_${n}_0.nek5000"
            job_file="${job_path}/${example}/${n}_re_${re}_pe_${pe}.sh"

            if [[ "$1" == "mesh" && ! -f ${output_file} ]]; then
                ./mesh.sh -b 0 4 0 1 0 1 $Nx $Ny $Nz \
                    -o ${data_path} -f ${mesh_pattern}_${Nx}x${Ny}x${Nz}.nmsh
            fi

            # Create the cases from the templates
            cp ${example_path}/${example}/case.template ${case_file}
            sed -i "s|\"Re\": .*|\"Re\": ${re}.0,|g" ${case_file}
            sed -i "s|\"Pe\": .*|\"Pe\": ${pe}.0,|g" ${case_file}
            sed -i "s|\"mesh_file\": .*|\"mesh_file\": \"${output_file}\",|g" ${case_file}
            sed -i "s|\"cache_file\": .*|\"cache_file\": \"${cache_file%0.nek5000}\",|g" ${case_file}
            sed -i "s|\"radius\": .*|\"radius\": ${element_size}|g" ${case_file}
            sed -i "s|\"element_size\": .*|\"element_size\": ${element_size},|g" ${case_file}

            # Rule of thumb, Brinkman penalization should be 1e6 / Re.
            brinkman=$(printf "%.2e" $(echo "scale=2; 1000000 / $re" | bc -l))
            sed -i "s|\"limits\": .*|\"limits\": [0.0, ${brinkman}],|g" ${case_file}

            cp ${example_path}/${example}/job.template ${job_file}
            echo $(basename ${job_file}) >> ${job_path}/${example}/.gitignore
            sed -i "s|#SBATCH --nodes=.*|#SBATCH --nodes=${nodes}|g" ${job_file}

            if [[ "$1" == "mesh" && ! -f "${cache_file}" ]]; then
                end_time=$(grep -oP '"end_time": \K[0-9.eE+-]+' ${case_file})
                sed -i "s|\"end_time\": .*|\"end_time\": 1e-3,|g" ${case_file}
                ./run.sh -c static_mixer/${example}/${n}_re_${re}_pe_${pe}.case
                sed -i "s|\"end_time\": .*|\"end_time\": ${end_time},|g" ${case_file}
            fi

            if [ ! -f "${cache_file}" ]; then
                rm -fr "${case_file}" "${job_file}"
            fi
        done
    done
done