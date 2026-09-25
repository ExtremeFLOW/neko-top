#!/bin/bash
# ============================================================================ #
# Define the help function
function help() {
    echo -e "Construction of a scaling example based on an existing example."

    echo -e "\n\e[4mUsage:\e[0m"
    echo -e "  $0 <scaling_example> <base_example> <num_processes...>"

    echo -e "\n\e[4mDescription:\e[0m"
    echo -e "  This script creates a scaling example based on an existing example."
    echo -e "  It copies the files from the base example to the scaling example"
    echo -e "  and creates a job script for the scaling example."
    echo -e "  The job script is created based on the job script of the base example."
    echo -e "  The number of processes is set in the job script."

    printf "\n\e[4mOptions:\e[0m\n"
    printf "  %-20s %s\n" "-h, --help" "Show this help message and exit"
    printf "  %-20s %s\n" "-c, --cluster" "Cluster to use (default: empty)"
    printf "  %-20s %s\n" "--clean" "Clean the output directory"

}

# ============================================================================ #
# Input handling

# List possible options
OPTIONS=help,cluster:,clean
OPT=h,c:

# Parse the inputs for options
PARSED=$(getopt --options=$OPT --longoptions=$OPTIONS --name "$0" -- "$@")
eval set -- "$PARSED"

# Default values and global variables
MAIN_DIR=$(dirname "$(dirname "$(realpath "$0")")")
CURRENT_DIR=$(pwd)
CLUSTER=""
CLEAN=false

# Loop through the options and set the variables
while true; do
    case "$1" in
    "-h" | "--help") help && exit ;;             # Print help
    "-c" | "--cluster") CLUSTER=$2 && shift 2 ;; # Clean logs

    # Long options only
    "--clean") CLEAN=true && shift 1 ;; # Clean logs

    # End of options
    "--") shift && break ;;
    esac
done

if [ $# -lt 3 ]; then
    echo "Usage: $0 <scaling_example> <base_example> <num_processes...>" >&2
    echo "Example: $0 passive_scalar passive_scalar 1 2 4 8" >&2
    echo "Call with -h or --help for more information." >&2
    exit 1
fi

# Extract the required inputs
SCALING_EXAMPLE=$1
BASE_EXAMPLE=$2
# We assume the rest of the arguments are the number of processes
shift 2
if [ $# -eq 1 ]; then
    N=($1)
else
    N=("$@")
fi

# ============================================================================ #
# Setup the parameters for the scaling example

EXAMPLE_DIR="$MAIN_DIR/examples"
SCRIPT_DIR="$MAIN_DIR/scripts/jobscripts"

CASE_FILE=$(find "$EXAMPLE_DIR/$BASE_EXAMPLE" -type f \( -name "*.case" -o -name "*.json" \) | head -n 1)
if [ -z "$CASE_FILE" ]; then
    echo "No case file found in $EXAMPLE_DIR/$BASE_EXAMPLE"
    exit 1
fi

# Create the directory for the scaling example
if [[ -d "$EXAMPLE_DIR/$SCALING_EXAMPLE" && $CLEAN == true ]]; then
    rm -rf "$EXAMPLE_DIR/$SCALING_EXAMPLE"
elif [[ -d "$EXAMPLE_DIR/$SCALING_EXAMPLE" ]]; then
    printf "Directory $EXAMPLE_DIR/$SCALING_EXAMPLE already exists.\n" >&2
    printf "Use --clean to clean it.\n" >&2
    exit 1
fi
mkdir -p "$EXAMPLE_DIR/$SCALING_EXAMPLE"

# Copy extra files to the example directory
find "$EXAMPLE_DIR/$BASE_EXAMPLE" -type f \
    \( -name "CMakeLists.txt" \
    -o -name "*.f90" \
    -o -name "*.sh" \
    -o -name "*.nmsh" \) \
    -exec cp {} $EXAMPLE_DIR/$SCALING_EXAMPLE \;

# Define a gitignore file in the scaling example directory
cat <<EOF >$EXAMPLE_DIR/$SCALING_EXAMPLE/.gitignore
# Ignore all files
*
# Ignore all directories
*/
EOF

# Set parameters for the cluster specific job script
if [ -n "$CLUSTER" ]; then
    # Set the scheduler and tasks per node based on the cluster
    if [ "$CLUSTER" == "LUMI-G" ]; then
        SCHEDULER=SLURM
        TASKS_PER_NODE=8
    else
        printf "ERROR: Unknown cluster $CLUSTER" >&2
        exit 1
    fi

    # Clean up the old job scripts
    rm -fr "$SCRIPT_DIR/$CLUSTER/$SCALING_EXAMPLE"
    mkdir -p "$SCRIPT_DIR/$CLUSTER/$SCALING_EXAMPLE"

    # Locate the base job script
    if [ -f "$SCRIPT_DIR/$CLUSTER/$BASE_EXAMPLE/$(basename $CASE_FILE).sh" ]; then
        JOBSCRIPT="$SCRIPT_DIR/$CLUSTER/$BASE_EXAMPLE/$(basename $CASE_FILE).sh"
    elif [ -f "$SCRIPT_DIR/$CLUSTER/$BASE_EXAMPLE/default.sh" ]; then
        JOBSCRIPT="$SCRIPT_DIR/$CLUSTER/$BASE_EXAMPLE/default.sh"
    elif [ -f "$SCRIPT_DIR/$CLUSTER/default.sh" ]; then
        JOBSCRIPT="$SCRIPT_DIR/$CLUSTER/default.sh"
    else
        printf "ERROR: No job script found for cluster $CLUSTER\n" >&2
        exit 1
    fi

    cp $EXAMPLE_DIR/$SCALING_EXAMPLE/.gitignore $SCRIPT_DIR/$CLUSTER/$SCALING_EXAMPLE
fi

# Loop over the number of processes
for n in "${N[@]}"; do
    if [[ ! $n =~ ^[0-9]+$ ]]; then
        printf "ERROR: Number of processes must be a positive integer\n" >&2
        printf "Got: '$n'\n" >&2
        exit 1
    fi

    # Create a case file for the number of processes
    cp $CASE_FILE $EXAMPLE_DIR/$SCALING_EXAMPLE/$n.case

    # Create the job script
    if [ -f "$JOBSCRIPT" ]; then
        JOBSCRIPT_N="$SCRIPT_DIR/$CLUSTER/$SCALING_EXAMPLE/$n.sh"
        cp "$JOBSCRIPT" "$JOBSCRIPT_N"

        # Assign the number of processes to the job script
        if [ $SCHEDULER == "SLURM" ]; then
            # Compute the number of nodes needed
            N_NODES=$((n / TASKS_PER_NODE))
            if [ $((n % TASKS_PER_NODE)) -ne 0 ]; then
                N_NODES=$((N_NODES + 1))
            fi

            sed -i "s/--nodes=[0-9]*/--nodes=$N_NODES/" "$JOBSCRIPT_N"
            sed -i "s/--ntasks=[0-9]*/--ntasks=$n/" "$JOBSCRIPT_N"
        else
            prinft "ERROR: Unknown scheduler $SCHEDULER" >&2
            exit 1
        fi
    fi
done
