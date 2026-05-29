#!/bin/bash
# ============================================================================ #
# Define the help function
function help() {
    echo -e "From Journal files to Binary Neko meshes."

    echo -e "\n\e[4mUsage:\e[0m"
    echo -e "  mesh.sh [options] [example]"

    echo -e "\n\e[4mDescription:\e[0m"
    echo -e "  This script automates the creation and conversions of meshes."
    echo -e "  Cubit is used to create meshes from the input files."
    echo -e "  The meshes are then converted to Nek5000 format using exo2nek."
    echo -e "  Finally, the meshes are converted to Neko format using rea2nbin."
    echo -e ""
    echo -e "  The user should specify the pattern used to find input files"
    echo -e "  in the INPUT_PATH folder, which defaults to data. Completed"
    echo -e "  meshes are stored in the OUTPUT_PATH folder, which defaults to"
    echo -e "  data_local."
    echo -e ""
    echo -e "  The script additionally support the creation of box meshes"
    echo -e "  directly from the Neko genmeshbox command."
    echo -e "  This is done by specifying the -b option followed by:"
    echo -e "  - Dimensions of the box: x0,x1,y0,y1,z0,z1 (Float)"
    echo -e "  - Number of elements in each direction: nx,ny,nz (Integer)"
    echo -e "  - Periodic boundary conditions: px,py,pz (Boolean 0/1 or false/true)"
    echo -e "  Example: mesh.sh -b 0.0 3.0 0.0 2.0 0.0 1.0 30 20 10 0 0 1"
    echo -e ""
    echo -e "Please note, periodic boundary conditions are not supported"
    echo -e "directly."
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
    printf "  -%-1s, --%-10s %-60s\n" "i" "input" "Input path for the mesh files."
    printf "  -%-1s, --%-10s %-60s\n" "o" "output" "Output path for the meshes."
    printf "  -%-1s, --%-10s %-60s\n" "f" "file" "Output file for the mesh."
    printf "  -%-1s, --%-10s %-60s\n" "d" "dimension" "Dimension of mesh file."
    printf "  -%-1s, --%-10s %-60s\n" "b" "box" "Create a box mesh."

    exit 0
}
if [ $# -lt 1 ]; then help; fi

# Assign default values to the options
ALL=false      # Run all meshing
KEEP=false     # Keep logs and temporaries
REMESH=false   # Do complete remesh
DIMENSION=2    # Dimension of GMSH file
PREPART=0      # Pre-partition the mesh
BOX=false      # Create a box mesh
OUTPUT_FILE="" # Output file for the mesh
OUTPUT_PATH="" # Path to the output meshes
INPUT_PATH=""  # Path to the input files

# List possible options
OPTIONS=help,all,keep,remesh,file:,dimension:,box,input:,output:,prepart:
OPT=h,a,k,r,f:,d:,b,i:,o:,p:

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
    "-b" | "--box") BOX=true && shift ;;               # Create a box mesh

    "-i" | "--input") INPUT_PATH="$2" && shift 2 ;;   # Input path for the mesh files
    "-o" | "--output") OUTPUT_PATH="$2" && shift 2 ;; # Output path for the meshes
    "-f" | "--file") OUTPUT_FILE="$2" && shift 2 ;;   # Output file for the mesh

    "-p" | "--prepart") PREPART="$2" && shift 2 ;;     # Prepart the mesh

    # End of options
    "--") shift && break ;;
    esac
done
export ALL KEEP REMESH DIMENSION

# ============================================================================ #
# User defined inputs.

# Define all needed folders relative to the project folder. (without trailing /)
export CURRENT_DIR=$(pwd)
export MAIN_DIR=$(dirname $(realpath $0))
export EXTERNAL_DIR="$MAIN_DIR/external"

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

[ -f $MAIN_DIR/prepare.env ] && source $MAIN_DIR/prepare.env
source $MAIN_DIR/scripts/dependencies.sh
source $MAIN_DIR/scripts/meshing.sh

# ============================================================================ #
# Create a box mesh if requested

if [ "$BOX" == "true" ]; then
    # Create the box mesh using Neko genmeshbox
    find_neko $NEKO_DIR

    [ -z "$OUTPUT_FILE" ] && OUTPUT_FILE="box.nmsh"

    if [[ ! -f "$OUTPUT_PATH/$OUTPUT_FILE" || $REMESH == "true" ]]; then

        mkdir -p $OUTPUT_PATH/box_mesh.tmp
        cd $OUTPUT_PATH/box_mesh.tmp

        printf "\n\e[4mCreating box mesh.\e[0m\n"
        echo "Finding Neko in $NEKO_DIR"

        if ! command -v genmeshbox &>/dev/null; then
            echo "Error: genmeshbox command not found."
            echo "Please ensure Neko is installed and the path is set correctly."
            exit 1
        fi
        echo "genmeshbox $@"
        genmeshbox $@ 1>box_mesh.log 2>error.log

        if [ $? -ne 0 ]; then
            echo "Error: Failed to create box mesh."
            exit 1
        fi

        cp box.nmsh $OUTPUT_PATH/$OUTPUT_FILE
        printf '  %-11s %-67s\n' "Box Mesh:" "Created in $OUTPUT_FILE"
    else
        printf '  %-11s %-67s\n' "Box Mesh:" "Already exists, skipping."
    fi

    if [[ $PREPART -gt 1 && ! -f $OUTPUT_PATH/${OUTPUT_FILE%.*}_$PREPART.nmsh ]]; then
        prepart $OUTPUT_PATH/$OUTPUT_FILE $PREPART
    fi

    cd $CURRENT_DIR

    # Clean up the temporary files
    [ $KEEP == "false" ] && rm -fr $OUTPUT_PATH/box_mesh.tmp
    exit 0
fi

# ============================================================================ #
# Loop through the inputs and extract the file_list

SUPPORTED_TYPES=(".jou" ".e" ".exo" ".rea" ".re2" ".geo")

file_list=""
for input in $@; do
    [[ $ALL == "true" ]] && break

    if [ -f "$input" ]; then
        input="$(realpath $input)"
        if [ -z "${input#"$INPUT_PATH"/}" ]; then
            printf '  %-10s %-67s\n' "Not found in INPUT_PATH:" "$input"
            continue
        fi
        file_list+="$input"
        continue
    fi

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
            printf '  %-10s %-67s\n' "Skipping:" \
                "$input_dir/$input_name.$input_type"
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
    "e" | "exo") exo2nbin $input_file 1>${input_name%.*}.log 2>error.log ;;
    "rea" | "re2") re2nbin $input_file 1>${input_name%.*}.log 2>error.log ;;
    "geo") geo2nbin $input_file 1>${input_name%.*}.log 2>error.log ;;
    "msh") msh2nbin $input_file 1>${input_name%.*}.log 2>error.log ;;
    esac

    if [ $? -ne 0 ]; then
        printf '  %-10s %-67s\n' "Error:" "Mesh not created: $input_file"
        cd $CURRENT_DIR
        continue
    fi

    if [ -n "$OUTPUT_FILE" ]; then
        mkdir -p $(dirname $OUTPUT_PATH/$OUTPUT_FILE)
        cp -f $input_name.nmsh $OUTPUT_PATH/$OUTPUT_FILE
    else
        find -name "*.nmsh" -exec cp -t $OUTPUT_PATH/$input_dir {} \;
    fi
    cd $CURRENT_DIR

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
