#!/bin/bash
# omg Tim you're going to read this eventually and this it's such a stupid way of doing this hahahaha

# Mesh		4
# dt			10
# Re			25
# Chi			69
# implcit	72
# radius		76

make_a_case() {
    # inputs:
    # 1) mesh_list
    local -n _mesh_list=$1
    # 2) Re_list
    local -n _Re_list=$2
    # 3) chi_list
    local -n _chi_list=$3
    # 4) implicit_list
    local -n _implicit_list=$4
    # 5) radius_list
    local -n _radius_list=$5
    # 6) name_list
    big_name=$6

    rm -r $big_name 2>/dev/null
    mkdir -p $big_name

    for mesh in "${_mesh_list[@]}"; do
        for Re in "${_Re_list[@]}"; do
            for chi in "${_chi_list[@]}"; do
                for implicit in "${_implicit_list[@]}"; do
                    for radius in "${_radius_list[@]}"; do
                        name=$big_name/mesh${mesh}_Re${Re}_chi${chi}_radius${radius}_$implicit
                        echo $name
                        cp -r default_case $name
                        mv $name/case.template $name/cylinder.case
                        # now all the replacements

                        new_line='4s#.*#"mesh_file":"data_local/brinkman_parameters/immersed_M'$mesh'.nmsh",#'
                        sed -i $new_line $name/cylinder.case
                        new_line='25s#.*#"Re":'$Re',#'
                        sed -i $new_line $name/cylinder.case
                        new_line='69s#.*#'$chi',#'
                        sed -i $new_line $name/cylinder.case
                        new_line='72s#.*#"implicit":'$implicit',#'
                        sed -i $new_line $name/cylinder.case

                        # radius is also strange if we want no filtering
                        # (not exactly the same as filter radius = 0)
                        if [ $radius = "0" ]; then
                            # note filter type is line 75
                            new_line='75s#.*#"type":"none",#'
                        else
                            new_line='76s#.*#"radius":'$radius',#'
                        fi
                        sed -i $new_line $name/cylinder.case

                        # time step is a bit strange
                        case $mesh in
                        "2")
                            dt="2.5e-3"
                            ;;
                        "3")
                            dt="1.6e-3"
                            ;;
                        "4")
                            dt="1.25e-3"
                            ;;
                        esac
                        new_line='10s#.*#"timestep":'$dt',#'
                        sed -i $new_line $name/cylinder.case

                    done
                done
            done
        done
    done
}

# CASES
# ----------------------------------------------#

# Implementation
mesh_list=("2")
Re_list=("200")
chi_list=("1" "100" "1000")
implicit_list=("true" "false")
radius_list=("0.05")

case_name="Implementation"
make_a_case mesh_list Re_list chi_list implicit_list radius_list $case_name

# I still don't know exactly what the remaining cases should be..

# Can we do the filter radius??
mesh_list=("2")
Re_list=("200")
chi_list=("1000")
implicit_list=("true")
radius_list=("0" "0.01" "0.05" "0.1")

case_name="Filter_radius"
make_a_case mesh_list Re_list chi_list implicit_list radius_list $case_name

# I still don't know exactly what the remaining cases should be..

# omg Tim you're going to read this eventually and this it's such a stupid way of doing this hahahaha

# Mesh		4
# dt			10
# Re			25
# Chi			69
# implcit	72
# radius		76

# mesh study
declare -a mesh_list=("2" "3" "4")
declare -a Re_list=("200")

big_name="Mesh_study"
rm -r $big_name 2>/dev/null
mkdir $big_name

for mesh in "${mesh_list[@]}"; do
    for Re in "${Re_list[@]}"; do
        name=$big_name/mesh${mesh}_Re${Re}
        cp -r default_case_meshed $name
        mv $name/case.template $name/cylinder.case
        # now all the replacements

        new_line='4s#.*#"mesh_file":"data_local/brinkman_parameters/meshed_M'$mesh'.nmsh",#'
        sed -i $new_line $name/cylinder.case
        new_line='25s#.*#"Re":'$Re',#'
        sed -i $new_line $name/cylinder.case
        # time step is a bit strange
        case $mesh in
        "2")
            dt="2.5e-3"
            ;;
        "3")
            dt="1.6e-3"
            ;;
        "4")
            dt="1.25e-3"
            ;;
        esac
        new_line='10s#.*#"timestep":'$dt',#'
        sed -i $new_line $name/cylinder.case
    done
done

# Reynolds study
declare -a mesh_list=("3")
declare -a Re_list=("200" "1000" "2000" "3900")

big_name="Re_study"
rm -r $big_name 2>/dev/null
mkdir $big_name

for mesh in "${mesh_list[@]}"; do
    for Re in "${Re_list[@]}"; do
        name=$big_name/mesh${mesh}_Re${Re}
        cp -r default_case_meshed $name
        mv $name/case.template $name/cylinder.case
        # now all the replacements

        new_line='4s#.*#"mesh_file":"data_local/brinkman_parameters/meshed_M'$mesh'.nmsh",#'
        sed -i $new_line $name/cylinder.case
        new_line='25s#.*#"Re":'$Re',#'
        sed -i $new_line $name/cylinder.case
        # time step is a bit strange
        case $Re in
        "200") dt="2.50e-3" ;;
        "1000") dt="0.50e-3" ;;
        "2000") dt="0.25e-3" ;;
        "3900") dt="0.01e-3" ;;
        esac
        new_line='10s#.*#"timestep":'$dt',#'
        sed -i $new_line $name/cylinder.case
    done
done
