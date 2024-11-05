#!/bin/bash
# This script generates the cases for the Brinkman parameters study.
# The script is organized in the following way:
# 1) Define the function make_a_case, which generates the cases for a given
#    experiment.
# 2) Define the experiments, which are a set of cases with different parameters.
# 3) Call make_a_case for each experiment.

export ROOT_FOLDER=$(realpath $(dirname $0))

make_a_case() {
    # 1) Name of the experiment
    experiment_name=$1

    # inputs:
    # 2) Boundary method
    local -n _method_list=$2
    # 3) mesh_list
    local -n _mesh_list=$3
    # 4) Re_list
    local -n _Re_list=$4
    # 5) chi_list
    local -n _chi_list=$5
    # 6) implicit_list
    local -n _implicit_list=$6
    # 7) radius_list
    local -n _radius_list=$7

    # Define the location of the case
    folder=$(realpath $ROOT_FOLDER/$experiment_name)

    rm -r $folder 2>/dev/null
    mkdir -p $folder

    for method in "${_method_list[@]}"; do
        for mesh in "${_mesh_list[@]}"; do
            for Re in "${_Re_list[@]}"; do
                for chi in "${_chi_list[@]}"; do
                    for implicit in "${_implicit_list[@]}"; do
                        for radius in "${_radius_list[@]}"; do

                            # Build the name of the experiment
                            name=""

                            # If there is only one value in the list, don't
                            # include it in the name
                            [ ${#_method_list[@]} -gt 1 ] && name+=${method}_
                            [ ${#_mesh_list[@]} -gt 1 ] && name+=mesh_${mesh}_
                            [ ${#_Re_list[@]} -gt 1 ] && name+=re_${Re}_

                            # Case specific parameters
                            if [ "$method" == "brinkman" ]; then
                                [ ${#_chi_list[@]} -gt 1 ] && name+=chi_${chi}_
                                [ ${#_implicit_list[@]} -gt 1 ] && name+=implicit_${implicit}_
                                [ ${#_radius_list[@]} -gt 1 ] && name+=radius_${radius//./-}_
                            fi
                            name=${name%_}
                            echo "$experiment_name: $name"

                            # Create directory and copy the default files
                            [ -d $folder/$name ] && continue
                            mkdir -p $folder/$name
                            cp -t $folder/$name $ROOT_FOLDER/default_case/cylinder.f90

                            casefile=$folder/$name/cylinder.case
                            if [ $method == "brinkman" ]; then
                                cp $ROOT_FOLDER/default_case/brinkman.template $casefile
                            elif [ $method == "idw" ]; then
                                cp $ROOT_FOLDER/default_case/idw.template $casefile
                            elif [ $method == "meshed" ]; then
                                cp $ROOT_FOLDER/default_case/meshed.template $casefile
                            else
                                echo "Method not recognized"
                                exit 1
                            fi

                            # Set the timestep based on which mesh wass chosen

                            case $mesh in
                            2) dt_mesh="2.50E-03" ;;
                            3) dt_mesh="2.00E-03" ;;
                            4) dt_mesh="1.13E-03" ;;
                            *) dt_mesh="0" ;;
                            esac

                            if [ $method == "meshed" ]; then
                                # Use the meshed timestep directly
                                dt=$dt_mesh
                            else
                                # Compute the expected timestep, in exponential notation
                                dt=$(echo "scale=10; 50.0 / ($Re*$chi)" | bc)
                                # Compute min of the two
                                dt_mesh=$(printf "%f" $dt_mesh)
                                [ $(echo "$dt > $dt_mesh" | bc) -eq 1 ] && dt=$dt_mesh
                                dt=$(printf "%.2e" $dt)
                            fi

                            # Locate the pattern and replace it
                            if [ $method == "meshed" ]; then
                                mesh_pattern='"mesh_file": "data_local/brinkman_parameters/meshed_M2.nmsh"'
                                mesh_replacement='"mesh_file": "data_local/brinkman_parameters/meshed_M'$mesh'.nmsh"'
                            else
                                mesh_pattern='"mesh_file": "data_local/brinkman_parameters/immersed_MX.nmsh"'
                                mesh_replacement='"mesh_file": "data_local/brinkman_parameters/immersed_M'$mesh'.nmsh"'
                            fi

                            re_pattern='"Re": 200.0'
                            re_replacement='"Re": '$Re

                            chi_pattern='"limits": \[ 0.0, 100.0 \]'
                            chi_replacement='"limits": \[ 0.0, '$chi' \]'

                            implicit_pattern='"implicit": true'
                            implicit_replacement='"implicit": '$implicit

                            radius_pattern='"radius": 0.05'
                            radius_replacement='"radius": '$radius

                            timestep_pattern='"timestep": 2.5e-3'
                            timestep_replacement='"timestep": '$dt

                            sed -i "s#$mesh_pattern#$mesh_replacement#" $casefile
                            sed -i "s#$re_pattern#$re_replacement#" $casefile
                            sed -i "s#$chi_pattern#$chi_replacement#" $casefile
                            sed -i "s#$implicit_pattern#$implicit_replacement#" $casefile
                            sed -i "s#$radius_pattern#$radius_replacement#" $casefile
                            sed -i "s#$timestep_pattern#$timestep_replacement#" $casefile
                        done
                    done
                done
            done
        done
    done
}

# CASES
# ---------------------------------------------------------------------------- #
case_name="Implementation"

method_list=("brinkman")
mesh_list=("2")
Re_list=("200")
chi_list=("1" "100" "1000")
implicit_list=("true" "false")
radius_list=("0.05")

make_a_case $case_name method_list mesh_list Re_list chi_list implicit_list \
    radius_list

# ---------------------------------------------------------------------------- #
case_name="Filter_radius"

method_list=("brinkman")
mesh_list=("2")
Re_list=("200")
chi_list=("1000")
implicit_list=("true")
radius_list=("0" "0.01" "0.05" "0.1")

make_a_case $case_name method_list mesh_list Re_list chi_list implicit_list \
    radius_list

# ---------------------------------------------------------------------------- #
case_name="Re_study"

method_list=("brinkman" "meshed" "idw")
mesh_list=("2")
Re_list=("200" "400" "1000" "2000" "3900")
chi_list=("1000")
implicit_list=("true")
radius_list=("0.1")

make_a_case $case_name method_list mesh_list Re_list chi_list implicit_list \
    radius_list

# ---------------------------------------------------------------------------- #
case_name="Mesh_study"

method_list=("brinkman" "meshed" "idw")
mesh_list=("2" "3" "4")
Re_list=("200" "1000")
chi_list=("1000")
implicit_list=("true")
radius_list=("0.1")

make_a_case $case_name method_list mesh_list Re_list chi_list implicit_list \
    radius_list
