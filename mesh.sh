#!/bin/bash
# ============================================================================ #
# Define the help function
function help() {
    echo -e "From Journal files to Binary Neko meshes."

    echo -e "\n\e[4mUsage:\e[0m"
    echo -e "  mesh.sh [options] [example]"

    echo -e "\n\e[4mDescription:\e[0m"
    echo -e "  This script automates the creation and convertions of meshes."
    echo -e "  Cubit is used to create meshes from the journal files."
    echo -e "  The meshes are then converted to Nek5000 format using exo2nek."
    echo -e "  Finally, the meshes are converted to Neko format using rea2nbin."
    echo -e ""
    echo -e "  The user should specify the pattern used to find journal files"
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
    printf "  -%-1s, --%-10s %-60s\n" "a" "all" "Run all journals available."
    printf "  -%-1s, --%-10s %-60s\n" "k" "keep" "Keep logs and temporaries."
    printf "  -%-1s, --%-10s %-60s\n" "d" "delete" "Delete previous runs."
    printf "  -%-1s, --%-10s %-60s\n" "r" "remesh" "Do complete remesh."

    exit 0
}

if [ $# -lt 1 ]; then help; fi

# ============================================================================ #
# User defined inputs.

# Define all needed folders relative to the project folder. (without trailing /)
CURRENT_DIR=$(pwd)
MAIN_DIR=$(dirname $(realpath $0))
INPUT_PATH="$MAIN_DIR/data"        # Journal files
OUTPUT_PATH="$MAIN_DIR/data_local" # Meshes

# Set the path to the external dependencies
if [ -z "$NEK5000_DIR" ]; then
    NEK5000_DIR="$MAIN_DIR/external/Nek5000"
fi
if [ -z "$NEKO_DIR" ]; then
    NEKO_DIR="$MAIN_DIR/external/neko"
fi
if [ -z "$JSON_FORTRAN_DIR" ]; then
    JSON_FORTRAN_DIR="$MAIN_DIR/external/json-fortran"
fi

# If the user specified a relative path, then make it absolute
if [[ ! "$NEK5000_DIR" = /* ]]; then
    NEK5000_DIR=$(realpath $NEK5000_DIR)
fi
if [[ ! "$NEKO_DIR" = /* ]]; then
    NEKO_DIR=$(realpath $NEKO_DIR)
fi
if [[ ! "$JSON_FORTRAN_DIR" = /* ]]; then
    JSON_FORTRAN_DIR=$(realpath $JSON_FORTRAN_DIR)
fi

# ============================================================================ #
# Ensure executables are available

# Check if Cubit are available
if [ ! -z "$(which cubit)" ]; then
    cubit=$(which cubit)
elif [ ! -z "$(which coreform_cubit)" ]; then
    cubit=$(which coreform_cubit)
elif [ ! -z "$(which trelis)" ]; then
    cubit=$(which trelis)
else
    echo "Cubit not found. Please ensure it is installed and available in the PATH."
    exit 1
fi
cubit="$cubit -nojournal -nographics -batch -noecho"

source $MAIN_DIR/scripts/dependencies.sh
find_nek5000 $NEK5000_DIR

# Check if exo2nek is available
if [ ! -z "$(which exo2nek)" ]; then
    exo2nek=$(which exo2nek)
elif [ -f "$NEK5000_DIR/bin/exo2nek" ]; then
    exo2nek="$NEK5000_DIR/bin/exo2nek"
elif [ -f "$NEK5000_DIR/tools/maketools" ]; then
    cd $NEK5000_DIR/tools
    ./maketools exo2nek
    cd $CURRENT_DIR
    exo2nek="$NEK5000_DIR/bin/exo2nek"
else
    echo "exo2nek not found. Please ensure it is installed and available in the PATH."
    exit 1
fi

# Check if rea2nbin is available
if [ ! -z "$(which rea2nbin)" ]; then
    rea2nbin=$(which rea2nbin)
elif [ -f "$NEKO_DIR/bin/rea2nbin" ]; then
    rea2nbin="$NEKO_DIR/bin/rea2nbin"
else
    echo "rea2nbin not found. Please ensure it is installed and available in the PATH."
    exit 1
fi

# Set the LD_LIBRARY_PATH to ensure the libraries are found
if [ -d "$JSON_FORTRAN_DIR/lib" ]; then
    export LD_LIBRARY_PATH="$JSON_FORTRAN_DIR/lib:$LD_LIBRARY_PATH"
elif [ -d "$JSON_FORTRAN_DIR/lib64" ]; then
    export LD_LIBRARY_PATH="$JSON_FORTRAN_DIR/lib64:$LD_LIBRARY_PATH"
fi

# ============================================================================ #
# Handle inputs and settings

journals=""
for in in $@; do

    # Opions and flags
    if [[ ${in:0:2} == "--" ]]; then
        case "${in:2}" in
        "help") help ;;          # Print help
        "all") ALL=true ;;       # Run all journals available
        "keep") KEEP=true ;;     # Clean logs
        "delete") DELETE=true ;; # Delete previous runs
        "remesh") REMESH=true ;; # Do complete remesh
        \?) printf '  %-10s %-67s\n' "Invalid option:" "$in" && exit 1 ;;
        esac

    elif [[ ${in:0:1} == "-" ]]; then
        for ((i = 1; i < ${#in}; i++)); do
            case "${in:$i:1}" in
            "h") help ;;        # Print help
            "a") ALL=true ;;    # Run all journals available
            "k") KEEP=true ;;   # Clean logs
            "d") DELETE=true ;; # Delete previous runs
            "r") REMESH=true ;; # Do complete remesh
            \?) printf '  %-10s %-67s\n' "Invalid option:" "${in:$i:1}" && exit 1 ;;
            esac
        done

    # Ignore invalid inputs
    elif [[ -z $(ls $INPUT_PATH/$in 2>>/dev/null) ]]; then
        printf '  %-10s %-67s\n' "Not Found:" "$in"

    # Extract the journals from the input
    elif [[ ! $ALL ]]; then
        for journal in $(find $INPUT_PATH/$in -name "*.jou"); do
            journals+="$journal "
        done
    fi
done

if [ $ALL ]; then
    journals=""
    for journal in $(find $INPUT_PATH -name "*.jou" 2>>/dev/null); do
        journals+="$journal "
    done
fi

if [ ! "$journals" ]; then
    exit
fi

if [ $DELETE ]; then
    printf 'Do you wish to delete ALL datasets? [Yes No]\n'
    read -p '> ' yn
    case $yn in
    [Yy]*) echo "Removing..." && rm -fr $OUTPUT_PATH && echo "Results removed" ;;
    *) echo "Results not removed" ;;
    esac
    printf 'Logs have been cleaned.\n'
    rm -fr $LPATH
fi

# ============================================================================ #
# Function to run meshing

function mesh() {
    journal_name=$(basename $1)

    rm -fr ${journal_name%.*}.log error.log

    # Construct the mesh using Cubit
    $cubit $1 1>>${journal_name%.*}.log 2>>error.log

    if [[ ! -s error.log && $(find ./ -name "*.exo" | wc -l) -lt 1 ]]; then
        printf "\n\e[4mError:\e[0m\n" 2>>error.log
        printf "  %-10s %-67s\n" "Cubit:" "Exodus not created: $journal" 2>>error.log
        exit 1
    fi

    # Convert the mesh from Exodus format to Nek5000 format
    (
        echo "$(ls *.exo 2>>error.log | wc -l)"
        for file in *.exo; do echo "${file%.*}"; done
        echo "0"
        echo "0"
        echo "${journal_name%.*}"
    ) | $exo2nek 1>>${journal_name%.*}.log 2>>error.log

    if [[ ! -s error.log && $(find ./ -name "*.re2" | wc -l) -lt 1 ]]; then
        printf "\n\e[4mError:\e[0m\n" 2>>error.log
        printf "  %-10s %-67s\n" "exo2nek:" "re2 not created: $journal" 2>>error.log
        exit 1
    fi

    # Convert the mesh to Neko mesh format
    for file in *.re2; do
        $rea2nbin $file ${file%.*}.nmsh 1>>${journal_name%.*}.log 2>>error.log
    done

    if [[ ! -s error.log && $(find ./ -name "*.nmsh" | wc -l) -lt 1 ]]; then
        printf "\n\e[4mError:\e[0m\n" 2>>error.log
        printf "  %-10s %-67s\n" "rea2nbin:" "re2 not created: $journal" 2>>error.log
        exit 1
    fi
}

# ============================================================================ #
# Run the journals

full_start=$(date +%s.%N)

mkdir -p $OUTPUT_PATH $LPATH
printf "\n\e[4mQueueing journals.\e[0m\n"
for journal in $journals; do
    journal_name=${journal#$INPUT_PATH/}
    journal_dir=$(dirname $journal_name)

    if [ -f "$OUTPUT_PATH/${journal_name%.*}.nmsh" ]; then
        if [ $REMESH ]; then
            printf '  %-10s %-67s\n' "Remeshing:" "$journal_name"
            rm -f $OUTPUT_PATH/${journal_name%.*}.nmsh
        else
            printf '  %-10s %-67s\n' "Skipping:" "$journal_name"
            continue
        fi
    else
        printf '  %-10s %-67s\n' "Meshing:" "$journal_name"
    fi

    mkdir -p $OUTPUT_PATH/$journal_dir/tmp
    cd $OUTPUT_PATH/$journal_dir/tmp

    mesh $journal

    cp ./*.nmsh -ft $OUTPUT_PATH/$journal_dir
    if [ ! $KEEP ]; then
        rm -fr $OUTPUT_PATH/$journal_dir/tmp
    fi
    cd $CURRENT_DIR
done

# ============================================================================ #
# Print the results

for journal in $journals; do
    journal_name=$(basename $journal)
    journal_dir=$(dirname $journal)

    if [ -f "$INPUT_PATH/${journal%.*}.log" ]; then
        printf '  %-10s %-67s\n' "Error:" "Mesh not created: $journal"
    fi

done

printf "\n"
# # EOF # #
