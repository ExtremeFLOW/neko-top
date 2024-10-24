#!/bin/bash
# ============================================================================ #
# Define the help function
function help() {
    echo -e "From Journal files to Binary Neko meshes."

    echo -e "\n\e[4mUsage:\e[0m"
    echo -e "  mesh.sh [options] [example]"

    echo -e "\n\e[4mDescription:\e[0m"
    echo -e "  This script automates the creation and convertions of meshes."
    echo -e "  Cubit is used to create meshes from the input files."
    echo -e "  The meshes are then converted to Nek5000 format using exo2nek."
    echo -e "  Finally, the meshes are converted to Neko format using rea2nbin."
    echo -e ""
    echo -e "  The user should specify the pattern used to find input files"
    echo -e "  in the INPUT_PATH folder, which defaults to data. Completed"
    echo -e "  meshes are stored in the OUTPUT_PATH folder, which defaults to"
    echo -e "  data_local."
    echo -e ""
    echo -e "  This code depends on the following external libraries:"
    echo -e "    - Cubit"
    echo -e "    - Nek5000"
    echo -e "    - Neko"
    echo -e "    - json-fortran"
    echo -e ""

    printf "\e[4mOptions:\e[0m\n"
    printf "  -%-1s, --%-10s %-60s\n" "h" "help" "Print help."
    printf "  -%-1s, --%-10s %-60s\n" "a" "all" "Run all input files available."
    printf "  -%-1s, --%-10s %-60s\n" "k" "keep" "Keep logs and temporaries."
    printf "  -%-1s, --%-10s %-60s\n" "r" "remesh" "Do complete remesh."
    printf "  -%-1s, --%-10s %-60s\n" "d" "dimension" "Dimension of mesh file."

    exit 0
}
if [ $# -lt 1 ]; then help; fi

# Assign default values to the options
ALL=false    # Run all meshing
KEEP=false   # Keep logs and temporaries
REMESH=false # Do complete remesh
DIMENSION=2  # Dimension of GMSH file

# List possible options
OPTIONS=help,all,keep,remesh,dimension:
OPT=h,a,k,r,d:

# Parse the inputs for options
PARSED=$(getopt --options=$OPT --longoptions=$OPTIONS --name "$0" -- "$@")
eval set -- "$PARSED"

# Loop through the options and set the variables
while true; do
    case "$1" in
    "-h" | "--help") help && exit ;;                   # Print help
    "-a" | "--all") ALL=true && shift ;;               # Run all file_list available
    "-k" | "--keep") KEEP=true && shift ;;             # Keep logs and temporaries
    "-r" | "--remesh") REMESH=true && shift ;;         # Do complete remesh
    "-d" | "--dimension") DIMENSION="$2" && shift 2 ;; # Dimension of GMSH file

    # End of options
    "--") shift && break ;;
    esac
done
export ALL KEEP REMESH DIMENSION

# ============================================================================ #
# User defined inputs.

# Define all needed folders relative to the project folder. (without trailing /)
CURRENT_DIR=$(pwd)
MAIN_DIR=$(dirname $(realpath $0))

# Set the path to the input files and the output meshes
[ -z $INPUT_PATH ] && INPUT_PATH="$MAIN_DIR/data"         # Input files
[ -z $OUTPUT_PATH ] && OUTPUT_PATH="$MAIN_DIR/data_local" # Meshes

# Set the path to the external dependencies
[ -z "$NEK5000_DIR" ] && NEK5000_DIR="$MAIN_DIR/external/Nek5000"
[ -z "$NEKO_DIR" ] && NEKO_DIR="$MAIN_DIR/external/neko"
[ -z "$JSON_FORTRAN_DIR" ] && JSON_FORTRAN_DIR="$MAIN_DIR/external/json-fortran"

# Force the paths to be absolute
export INPUT_PATH=$(realpath $INPUT_PATH)
export OUTPUT_PATH=$(realpath $OUTPUT_PATH)
export NEK5000_DIR=$(realpath $NEK5000_DIR)
export NEKO_DIR=$(realpath $NEKO_DIR)
export JSON_FORTRAN_DIR=$(realpath $JSON_FORTRAN_DIR)

# ============================================================================ #
# Ensure executables are available

source $MAIN_DIR/scripts/dependencies.sh
source $MAIN_DIR/scripts/meshing.sh

# ============================================================================ #
# Loop through the inputs and extract the file_list

SUPPORTED_TYPES=(".jou", ".exo", ".re2", ".geo")

file_list=""
for input in $@; do
    [[ $ALL == "true" ]] && break
    input_name="$(basename $input)"
    input_dir=$(realpath $INPUT_PATH/$(dirname $input))

    tmp_list=$(find $input_dir -type f -name "$input_name")

    # Ignore invalid inputs
    if [ -z "${tmp_list[@]}" ]; then
        printf '  %-10s %-67s\n' "Not Found:" "$input"
        continue
    fi

    # Extract the file_list from the input
    for type in ${SUPPORTED_TYPES[@]}; do
        for file in $tmp_list; do
            [[ $file == *$type ]] && file_list+="$file "
        done
    done
done

if [ "$ALL" == "true" ]; then
    file_list=""

    tmp_list=$(find $INPUT_PATH -type f)

    # Extract the file_list from the input
    for type in ${SUPPORTED_TYPES[@]}; do
        for file in $tmp_list; do
            [[ $file == *$type ]] && file_list+="$file "
        done
    done
fi
[ -z "$file_list" ] && exit 0

# ============================================================================ #
# Run the file_list

full_start=$(date +%s.%N)

mkdir -p $OUTPUT_PATH
printf "\n\e[4mQueueing file_list.\e[0m\n"

for input_file in $file_list; do
    input_name=$(basename ${input_file%.*})
    input_type=$(basename ${input_file##*.})
    input_dir=$(dirname ${input_file#$INPUT_PATH/})

    if [ -f "$OUTPUT_PATH/$input_dir/$input_name.nmsh" ]; then
        if [ $REMESH == "true" ]; then
            printf '  %-11s' "Remeshing:"
            rm -f $OUTPUT_PATH/${input_name%.*}.nmsh
        else
            printf '  %-11s' "Skipping:"
            continue
        fi
    else
        printf '  %-11s' "Meshing:"
    fi
    printf '%-67s\n' "$input_dir/$input_name.$input_type"

    tmp="$OUTPUT_PATH/$input_dir/$input_name.tmp"
    mkdir -p $tmp
    cd $tmp

    case $input_type in
    "jou") jou2nbin $input_file 1>${input_name%.*}.log 2>error.log ;;
    "exo") exo2nbin $input_file 1>${input_name%.*}.log 2>error.log ;;
    "re2") re2nbin $input_file 1>${input_name%.*}.log 2>error.log ;;
    "geo") geo2nbin $input_file 1>${input_name%.*}.log 2>error.log ;;
    esac

    cp *.nmsh -ft $OUTPUT_PATH/$input_dir

    # Clean up the temporary files
    [ $KEEP == "false" ] && rm -fr $tmp
done

# ============================================================================ #
# Print the results

for input_file in $file_list; do
    input_name=$(basename $input_file)
    input_dir=$(dirname $input_file)

    if [ -f "$INPUT_PATH/${input_file%.*}.log" ]; then
        printf '  %-10s %-67s\n' "Error:" "Mesh not created: $input_file"
    fi

done

printf "\n"
# # EOF # #
