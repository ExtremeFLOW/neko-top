#!/bin/bash
# omg Tim you're going to read this eventually and this it's such a stupid way of doing this hahahaha

# Mesh          4
# dt            10
# Re            25
# Chi           69
# implcit       72
# radius        76

make_a_case() {
    # inputs:
    # 1) Boundary method
    local -n _method_list=$1
    # 2) mesh_list
    local -n _mesh_list=$2
    # 3) Re_list
    local -n _Re_list=$3
    # 4) chi_list
    local -n _chi_list=$4
    # 5) implicit_list
    local -n _implicit_list=$5
    # 6) radius_list
    local -n _radius_list=$6
    # 7) name_list
    experiment_name=$7

    # Define the location of the case
    root_folder=$(realpath $(dirname $0))
    folder=$(realpath $root_folder/$experiment_name)

    rm -r $folder 2>/dev/null
    mkdir -p $folder

    for method in "${_method_list[@]}"; do
        for mesh in "${_mesh_list[@]}"; do
            for Re in "${_Re_list[@]}"; do
                for chi in "${_chi_list[@]}"; do
                    for implicit in "${_implicit_list[@]}"; do
                        for radius in "${_radius_list[@]}"; do

                            # Build the name of the experiment
                            name=$experiment_name/

                            # If there is only one value in the list, don't
                            # include it in the name
                            [ ${#_method_list[@]} -gt 1 ] && name+=${method}_
                            [ ${#_mesh_list[@]} -gt 1 ] && name+=mesh_${mesh}_
                            [ ${#_Re_list[@]} -gt 1 ] && name+=re_${Re}_
                            [ ${#_chi_list[@]} -gt 1 ] && name+=chi_${chi}_
                            [ ${#_implicit_list[@]} -gt 1 ] && name+=implicit_${implicit}_
                            [ ${#_radius_list[@]} -gt 1 ] && name+=radius_${radius//./-}_
                            name=${name%_}
                            echo $name

                            # Create directory and copy the default files
                            mkdir -p $folder/$name
                            cp -t $folder/$name $root_folder/default_case/cylinder.f90

                            casefile=$folder/$name/cylinder.case
                            if [ $method == "brinkman" ]; then
                                cp $root_folder/default_case/brinkman.template $casefile
                            elif [ $method == "idw" ]; then
                                cp $root_folder/default_case/idw.template $casefile
                            elif [ $method == "meshed" ]; then
                                cp $root_folder/default_case/meshed.template $casefile
                            else
                                echo "Method not recognized"
                                exit 1
                            fi

                            # now all the replacements

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

                            chi_pattern='"limits": \[ 0.0, 100.0'
                            chi_replacement='"limits": \[ 0.0, '$chi

                            implicit_pattern='"implicit": true'
                            implicit_replacement='"implicit": '$implicit

                            radius_pattern='"radius": 0.05'
                            radius_replacement='"radius": '$radius

                            sed -i "s#$mesh_pattern#$mesh_replacement#" $casefile
                            sed -i "s#$re_pattern#$re_replacement#" $casefile
                            sed -i "s#$chi_pattern#$chi_replacement#" $casefile
                            sed -i "s#$implicit_pattern#$implicit_replacement#" $casefile
                            sed -i "s#$radius_pattern#$radius_replacement#" $casefile
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

make_a_case method_list mesh_list Re_list chi_list implicit_list radius_list $case_name

# ---------------------------------------------------------------------------- #

case_name="Filter_radius"

method_list=("brinkman")
mesh_list=("2")
Re_list=("200")
chi_list=("1000")
implicit_list=("true")
radius_list=("0" "0.01" "0.05" "0.1")

make_a_case method_list mesh_list Re_list chi_list implicit_list radius_list $case_name

# ---------------------------------------------------------------------------- #
case_name="Re_study"

method_list=("brinkman" "meshed" "idw")
mesh_list=("2")
Re_list=("200" "400" "1000")
chi_list=("1000")
implicit_list=("true")
radius_list=("0.1")

make_a_case method_list mesh_list Re_list chi_list implicit_list radius_list $case_name
