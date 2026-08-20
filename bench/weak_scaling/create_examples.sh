#!/bin/bash
# This script is designed to construct a sequence of experiments which can be
# used to benchmark the code. It will create all necessary files and folders for
# the experiment. The user will specify a tag to keep track of experiments.
#
# Usage: create_examples.sh [options] Tag
# Options:
#   -h, --help        Show this help message and exit
#   -m, --mesh        Create the meshes for the cases (if not already created)

# ============================================================================ #
# Print the help message
function help() {
    echo -e "Usage: $0 [options] Tag"
    echo -e "Options:"
    echo -e "\t-h, --help        Show this help message and exit"
    echo -e "\t-m, --mesh        Create the meshes for the cases (if not already created)"
    echo -e "\t-c, --cluster     HPC cluster to create cases for (LUMI-G)"
    echo -e ""
    echo -e "This script is designed to construct a sequence of experiments which can be"
    echo -e "used to benchmark the code. It will create all necessary files and folders for"
    echo -e "the experiment. The user will specify a tag to keep track of experiments."
}

# ============================================================================ #
# Set general options and directories

export MAIN_DIR=$(realpath $(dirname $0)/../../)
export BENCH_DIR=$(realpath $(dirname $0))
if [ "$(basename $MAIN_DIR)" != "neko-top" ]; then
    echo "Unable to determine MAIN_DIR" >&2
    exit 1
fi

# Assign default values to the options
MESH="false"
CLUSTER="local"

# List possible options
OPTIONS=help,mesh,cluster:
OPT=h,m,c:

# Parse the inputs for options
PARSED=$(getopt --options=$OPT --longoptions=$OPTIONS --name "$0" -- "$@")
eval set -- "$PARSED"

# Loop through the options and set the variables
while true; do
    case "$1" in
    "-h" | "--help") help && exit ;;                  # Print help
    "-m" | "--mesh") MESH="true" && shift ;;          # Create the meshes
    "-c" | "--cluster") CLUSTER="$2" && shift 2 ;;    # HPC cluster

    # End of options
    "--") shift && break ;;
    esac
done
export MESH CLUSTER

if [ $# -lt 1 ]; then
    echo "ERROR: Tag is required" >&2
    exit 1
fi
TAG="$1"

# Make sure Cluster is not case-sensitive
CLUSTER="${CLUSTER^^}"

# Check if the cluster is supported and load the corresponding modules
if [ "${CLUSTER}" == "LUMI-G" ]; then
    command -v python >/dev/null 2>&1 || module load cray-python
else
    echo "Unsupported cluster: ${CLUSTER}" >&2
    exit 1
fi

# ============================================================================ #
# Set up the directories and paths for the experiment

export mesh_pattern="mixer"

export hpc_path="${MAIN_DIR}/scripts/jobscripts"
export data_path="${MAIN_DIR}/data_local/static_mixer"
export example_path="${MAIN_DIR}/examples/benchmark"
export template_path="${BENCH_DIR}/templates"

export experiment_path="${MAIN_DIR}/results/benchmark/experiments/${TAG}/"

[ ! -d "${data_path}" ] && mkdir -p "${data_path}"
[ ! -d "${example_path}/${TAG}" ] && mkdir -p "${example_path}/${TAG}"
[ ! -d "${experiment_path}" ] && mkdir -p "${experiment_path}"

# Add CMake file with the Tagged experiment path
if [ ! -f "${example_path}/CMakeLists.txt" ]; then
    touch "${example_path}/CMakeLists.txt"
fi
found=$(grep -c "add_subdirectory(${TAG})" "${example_path}/CMakeLists.txt")
if [ $found -eq 0 ]; then
    echo "add_subdirectory(${TAG})" >> "${example_path}/CMakeLists.txt"
fi

function create_case() {

    if [ $# -lt 6 ]; then
        echo "ERROR: Not enough input arguments" >&2
        exit 1
    fi

    local experiment="$1"     # Experiment type (weak_scaling)
    local Nx="$2"             # Mesh dimension in x
    local Ny="$3"             # Mesh dimension in y
    local Nz="$4"             # Mesh dimension in z
    local nodes="$5"          # Number of nodes

    [ $# -ge 6 ] && local N_memory=$6 || local N_memory=""
    [ $# -ge 7 ] && local keep_files=$7 || local keep_files=false
    [ $# -ge 8 ] && local end_time=$8 || local end_time=""

    # Compute element size, largest of 4/Nx, 1/Ny, 1/Nz
    element_size=$(python -c "print(max(4.0/${Nx}, 1.0/${Ny}, 1.0/${Nz}))")

    if [ "${CLUSTER}" == "LUMI-G" ]; then
        Np="8"
    else
        nodes=1
        Np="1"
    fi

    if [ -z "$N_memory" ]; then
        N_memory=$(grep -oP '(?<="n_memory": )[^,]*' "${template_path}/case.template")
    fi

    # Define a unique case name
    case_name="${CLUSTER,,}_nodes_${nodes}_mesh_${Nx}x${Ny}x${Nz}"
    case_name="${case_name}_n_memory_${N_memory}"
    [ "${keep_files}" == "true" ] && case_name="${case_name}_keep_files"
    [ -n "${end_time}" ] && case_name="${case_name}_end_time_${end_time}"

    # Set file names
    local mesh_file="${data_path}/${mesh_pattern}_${Nx}x${Ny}x${Nz}.nmsh"
    local case_file="${example_path}/${TAG}/${case_name}.case"
    local job_path="${hpc_path}/${CLUSTER}/benchmark/${TAG}"
    local experiment_file="${experiment_path}/${experiment}.csv"

    # Create the mesh if it does not exist
    if [[ "$MESH" == "true" ]]; then
        ${MAIN_DIR}/mesh.sh -b 0 4 0 1 0 1 $Nx $Ny $Nz \
           -o "${data_path}" -f ${mesh_pattern}_${Nx}x${Ny}x${Nz}.nmsh
    fi

    # Create the cases from the templates
    cp "${template_path}/case.template" "${case_file}"
    sed -i "s|\"mesh_file\": .*|\"mesh_file\": \"${mesh_file}\",|g" "${case_file}"
    sed -i "s|\"n_memory\": .*|\"n_memory\": ${N_memory},|g" "${case_file}"
    sed -i "s|\"keep_checkpoints\": .*|\"keep_checkpoints\": ${keep_files}|g" "${case_file}"
    sed -i "s|\"r\": .*|\"r\": ${element_size}|g" "${case_file}"
    [ -n "${end_time}" ] && sed -i "s|\"end_time\": .*|\"end_time\": ${end_time},|g" "${case_file}"

    # Create the jobscript
    if [ "${CLUSTER}" != "LOCAL" ]; then
        [ ! -d "${job_path}" ] && mkdir -p ${job_path}
        if [ ! -f ${job_path}/.gitignore ]; then
            echo "*.sh" > ${job_path}/.gitignore
            echo ".gitignore" >> ${job_path}/.gitignore
        fi
        cp "${template_path}/${CLUSTER}.sh" "${job_path}/${case_name}.sh"
        sed -i "s|--nodes=.*|--nodes=${nodes}|g" "${job_path}/${case_name}.sh"
    fi

    # Create the experiment entry
    if [ ! -f ${experiment_file} ]; then
        header="tag, case_name, cluster, Nx, Ny, Nz, nodes, Np"
        [ -n "${N_memory}" ] && header+=", N_memory"

        echo "$header" >> ${experiment_file}
    fi

    if ! grep -q "^${TAG}, ${case_name}\," "${experiment_file}"; then
        # Determine the data line
        data_line="${TAG}, ${case_name}, ${CLUSTER}, ${Nx}, ${Ny}, ${Nz}, ${nodes}, ${Np}"
        [ -n "${N_memory}" ] && data_line+=", ${N_memory}"
        printf "${data_line}\n" >> "${experiment_file}"
    fi
}

# Clean up old cases
if [ -d "${example_path}/${TAG}" ]; then
    find "${example_path}/${TAG}" -type f -name "*.case" -delete
fi
if [ -d "${experiment_path}" ]; then
    find "${experiment_path}" -type f -name "${TAG}_*.csv" -delete
fi
if [ -d "${hpc_path}/${CLUSTER}/benchmark/${TAG}" ]; then
    find "${hpc_path}/${CLUSTER}/benchmark/${TAG}" -type f -name "*.sh" -delete
fi

if [ ! -f "${example_path}/.gitignore" ]; then
    echo "*" > "${example_path}/.gitignore"
fi

# Update template files in case they have changed
rsync -u "${template_path}/mixer.f90" "${template_path}/CMakeLists.txt" \
    "${example_path}/${TAG}/"

experiment="single_node"
create_case ${experiment}  64  32  32  1 100
create_case ${experiment}  64  32  32  1 250
create_case ${experiment}  64  32  32  1 500
# create_case ${experiment}  64  32  32  1 1000 # OOM

create_case ${experiment} 128  32  32  1 100
create_case ${experiment} 128  32  32  1 250
# create_case ${experiment} 128  32  32  1 500  # OOM
# create_case ${experiment} 128  32  32  1 1000 # OOM

experiment="weak_scaling"
create_case ${experiment}  64  32  32  1 500
create_case ${experiment} 128  32  32  2 500
create_case ${experiment} 256  32  32  4 500
create_case ${experiment} 128  64  64  8 500
create_case ${experiment} 256  64  64 16 500

create_case ${experiment} 128  32  32  1 250
create_case ${experiment} 256  32  32  2 250
create_case ${experiment} 128  64  64  4 250
create_case ${experiment} 256  64  64  8 250
create_case ${experiment} 512  64  64 16 250
